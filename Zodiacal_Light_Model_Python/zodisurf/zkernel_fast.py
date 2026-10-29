"""
zkernel_fast.py
---------------
Vectorized, multithreaded line-of-sight integrator for the ZodiSURF model.

Gives the same results as zodisurf.zkernelpy.zkernel (same Simpson grid,
same float32 weights, same physics functions, same order of summing the
components), but evaluates all LOS samples in large NumPy batches instead
of one pixel at a time, and spreads those batches over threads.

All physics functions are imported from zkernelpy, so zkernelpy.py does not
need to be modified.
"""
import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np

# Important all the important physics functions from the original code :)
from zodisurf.zkernelpy_original import (earthsun, scattfunc, thermfunc, zsrcfunc,
                                zcloud, migband, solring, new_isocloud)

__all__ = ["zkernel_fast"]

RMAX = 5.2                                   # same as zkernelpy
_SIMP_H = 0.025                              # simpint default step
_H32 = np.float32(_SIMP_H)
_W_END = np.float32(0.333333) * _H32         # identical to simpint's float32 weights
_W_ODD = np.float32(1.333333) * _H32
_W_EVEN = np.float32(0.666667) * _H32
_DBWAVE = np.array([1.25, 2.2, 3.5, 4.9, 12., 25., 60., 100., 140., 240.])
_CHUNK_POINTS = 131072   # This splits certain calculations into "chunks" to save RAM
                         # I am trusting that this actually helps, per Claude's suggestion

# (density start index, source-function start index) for dust bands 1-4
_BANDS = ((45, 59), (71, 85), (97, 111), (123, 137))
# resonant ring / blob density parameter indices (same as zkernel)
_IRB = np.r_[149:170, 182]


# ----------------------------------------------------------------------
#  Batched Simpson grid (identical to zkernelpy.simpint)
# ----------------------------------------------------------------------
def _simpint_batch(nsteps):
    """
    Simpson abscissae/weights for many lines of sight at once.

    Returns
    -------
    los    : float32 (N,)  distance along each LOS
    wts    : float32 (N,)  Simpson weights
    owner  : int     (N,)  local LOS index of each element
    starts : int     (nlos,) first element of each LOS
    """
    nsteps = np.asarray(nsteps, dtype=np.int64)
    starts = np.cumsum(nsteps) - nsteps
    owner = np.repeat(np.arange(nsteps.size), nsteps)
    k = np.arange(owner.size, dtype=np.int64) - starts[owner]

    los = k.astype(np.float32) * _H32
    wts = np.where((k & 1) == 1, _W_ODD, _W_EVEN).astype(np.float32)
    wts[starts] = _W_END
    wts[starts + nsteps - 1] = _W_END
    return los, wts, owner, starts


# ----------------------------------------------------------------------
#  Scattering: vectorized call, or per-LOS fallback
# ----------------------------------------------------------------------
def _scatt_vector(phase_type, Det, Lambda, los, R, Re_e, SE_e, aScatt, sirr):
    out = scattfunc(phase_type, Det, Lambda, los, R, Re_e, SE_e, aScatt,
                    solar_irr=sirr)
    out = np.asarray(out, dtype=np.float64)
    if out.shape != R.shape:
        out = np.broadcast_to(out, R.shape).copy()
    return out


def _scatt_loop(phase_type, Det, Lambda, los, R, Re_e, SE_e, aScatt, sirr,
                starts, nst):
    """Calls scattfunc once per LOS with scalar Re/SolElong, as the original does."""
    out = np.empty(R.shape, dtype=np.float64)
    for s, n in zip(starts.tolist(), nst.tolist()):
        sl = slice(s, s + n)
        out[sl] = scattfunc(phase_type, Det, Lambda, los[sl], R[sl],
                            Re_e[s], SE_e[s], aScatt, solar_irr=sirr)
    return out


# ----------------------------------------------------------------------
#  Main kernel
# ----------------------------------------------------------------------
def zkernel_fast(data, a, phase_type='kelsall', L2=False, no_colcorr=True,
                 solar_irr=None, new_iso_comp=False, iso_comp_only=False,
                 group=None, n_jobs=1, scatt_mode='auto'):
    """
    Fast equivalent of zkernelpy.zkernel (without losinfo / partials).

    Parameters
    ----------
    data       : structured array from mk_zdata
    a          : (npar,), (nuniq, npar) or (npts, npar) parameter array
    group      : optional (npts,) int array mapping each pixel to a row of `a`.
                 Lets get_zmod build parameters once per unique wavelength
                 instead of once per pixel.
    n_jobs     : number of threads (default: 1; use -1 for all cores)
    scatt_mode : 'auto'   - use a vectorized scattfunc call if it gives the
                            same answer as the per-LOS call, else fall back
                 'vector' - always vectorized
                 'loop'   - call scattfunc once per LOS (safest, slower)

    The scatt mode actually used is stored in zkernel_fast.last_scatt_mode.
    """
    if phase_type not in ('kelsall', 'skysurf'):
        raise ValueError("phase_type must be 'kelsall' or 'skysurf'")
    if scatt_mode not in ('auto', 'vector', 'loop'):
        raise ValueError("scatt_mode must be 'auto', 'vector' or 'loop'")
    if n_jobs == -1 or n_jobs is None:
        n_jobs = os.cpu_count() or 1

    d2r = np.pi / 180.0
    npts = len(data)
    if npts == 0:
        return np.zeros(0)

    lam_all = np.asarray(data['wave_len'], dtype=np.float64)
    day = np.asarray(data['day1990'], dtype=np.float64)
    lon_deg = np.asarray(data['longitude'], dtype=np.float64)
    lat_deg = np.asarray(data['latitude'], dtype=np.float64)

    # ---- parameter rows / grouping ------------------------------------
    a = np.atleast_2d(np.asarray(a, dtype=np.float64))
    if group is None:
        if a.shape[0] == 1:
            group = np.zeros(npts, dtype=np.int64)
        elif a.shape[0] == npts:
            a, group = np.unique(a, axis=0, return_inverse=True)
        else:
            raise ValueError("zpar must have 1 row, npts rows, or be given with `group`.")
    group = np.asarray(group, dtype=np.int64).ravel()
    if group.shape != (npts,):
        raise ValueError("`group` must have one entry per data point.")

    if phase_type == 'skysurf' and not np.any(a[:, 183:189] != 0):
        raise ValueError("SKYSURF mode requires non-zero Hong parameters in zpar[183:189].")

    # Each unique (param row, wavelength[, solar_irr]) runs with scalar parameters
    cols = [group.astype(np.float64), lam_all]
    if solar_irr is not None:
        solar_irr = np.broadcast_to(
            np.asarray(solar_irr, dtype=np.float64).ravel(), (npts,))
        cols.append(solar_irr)
    keys, gid = np.unique(np.column_stack(cols), axis=0, return_inverse=True)
    gid = np.asarray(gid).ravel()

    gparams = []
    for g in range(len(keys)):
        prm = a[int(keys[g, 0])]
        Lambda = keys[g, 1]
        sirr = keys[g, 2] if solar_irr is not None else None
        if phase_type == 'kelsall':
            hit = np.nonzero(_DBWAVE == Lambda)[0]
            Det = int(hit[0]) if hit.size else -1
            aScatt = prm[1:10]
        else:
            Det = -99
            aScatt = prm[183:189]
        gparams.append((prm, Lambda, sirr, Det, aScatt))

    # ---- per-pixel geometry (vectorized) ------------------------------
    SolElong, Earth_Dis, Earth_Lon, Earth_Mean_Lon, L2_Dis = earthsun(day, lon_deg, lat_deg)
    Re = np.broadcast_to(np.asarray(L2_Dis if L2 else Earth_Dis, dtype=np.float64), (npts,))
    SolElong = np.broadcast_to(np.asarray(SolElong, dtype=np.float64), (npts,))
    Earth_Lon = np.broadcast_to(np.asarray(Earth_Lon, dtype=np.float64), (npts,))
    Earth_Mean_Lon = np.broadcast_to(np.asarray(Earth_Mean_Lon, dtype=np.float64), (npts,))

    lat = lat_deg * d2r
    lon = lon_deg * d2r
    COSLAT, SINLAT = np.cos(lat), np.sin(lat)
    COSLON, SINLON = np.cos(lon), np.sin(lon)
    X0 = Re * np.cos(Earth_Lon)
    Y0 = Re * np.sin(Earth_Lon)

    B = 2.0 * (X0 * COSLAT * COSLON + Y0 * COSLAT * SINLON)
    C = Re**2 - RMAX**2
    with np.errstate(divide='ignore', invalid='ignore'):
        Q = -0.5 * B * (1.0 + np.sqrt(B**2 - 4.0 * C) / np.abs(B))
        RANGE = np.maximum(Q, C / Q)
    RANGE = np.where(B == 0.0, np.sqrt(-C), RANGE)       # exact root when B == 0
    if not np.all(np.isfinite(RANGE)):
        raise ValueError("Non-finite line-of-sight range encountered; check inputs.")
    nsteps = np.rint(RANGE / _SIMP_H).astype(np.int64) + 1  # == int(round(b/h)) + 1

    # ---- work units: chunks of pixels with ~_CHUNK_POINTS LOS samples ----
    order = np.argsort(gid, kind='stable')
    bounds = np.cumsum(np.bincount(gid, minlength=len(keys)))[:-1]
    tasks = []
    for g, idx in enumerate(np.split(order, bounds)):
        if idx.size == 0:
            continue
        cum = np.cumsum(nsteps[idx])
        cuts = np.searchsorted(cum, np.arange(_CHUNK_POINTS, cum[-1], _CHUNK_POINTS))
        tasks.extend((g, c) for c in np.split(idx, cuts) if c.size)

    f = np.zeros(npts)
    state = {'mode': None if scatt_mode == 'auto' else scatt_mode}

    def _detect_mode(args, starts, nst):
        """Use vectorized scattfunc only if it reproduces the per-LOS call."""
        phase, Det, Lambda, los, R, Re_e, SE_e, aScatt, sirr = args
        try:
            sv = _scatt_vector(*args)
        except Exception:
            return 'loop'
        k = min(8, len(starts))
        m = int(starts[k - 1] + nst[k - 1])
        sl = _scatt_loop(phase, Det, Lambda, los[:m], R[:m], Re_e[:m], SE_e[:m],
                         aScatt, sirr, starts[:k], nst[:k])
        same = np.allclose(sv[:m], sl, rtol=1e-10, atol=0.0, equal_nan=True)
        return 'vector' if same else 'loop'

    def _run(task):
        g, pix = task
        prm, Lambda, sirr, Det, aScatt = gparams[g]

        nst = nsteps[pix]
        los, wts, own, starts = _simpint_batch(nst)
        p = pix[own]                                   # global pixel of each element

        Sxy = los * COSLAT[p]
        X = X0[p] + Sxy * COSLON[p]
        Y = Y0[p] + Sxy * SINLON[p]
        Z = los * SINLAT[p]
        R = np.sqrt(X**2 + Y**2 + Z**2)

        Re_e = Re[p]
        SE_e = SolElong[p]
        args = (phase_type, Det, Lambda, los, R, Re_e, SE_e, aScatt, sirr)

        mode = state['mode']
        if mode is None:                               # only happens in the first (serial) task
            mode = _detect_mode(args, starts, nst)
            state['mode'] = mode
        if mode == 'vector':
            Scatt = _scatt_vector(*args)
        else:
            Scatt = _scatt_loop(*args, starts, nst)

        Therm = thermfunc(phase_type, Det, Lambda, R, prm[10:14], no_colcorr=no_colcorr)

        # Components are added in the same order as zkernel. Components whose
        # normalisation is zero are skipped (the original adds an exact 0).
        Flux = np.zeros_like(R)
        Src_C = None

        if not iso_comp_only:
            if prm[14] != 0.0:                          # smooth cloud
                Src_C = zsrcfunc(Det, Scatt, Therm, prm[33:45], phase_type=phase_type)
                Flux += Src_C * zcloud(X, Y, Z, R, prm[14:33], FuncIndx=prm[0])
            for d0, s0 in _BANDS:                       # dust bands 1-4
                if prm[d0] != 0.0:
                    Dens, _ = migband(X, Y, Z, R, prm[d0:d0 + 14])
                    Flux += zsrcfunc(Det, Scatt, Therm, prm[s0:s0 + 12],
                                     phase_type=phase_type) * Dens
            aRB = prm[_IRB]                             # ring + blobs
            if aRB[0] != 0.0 or aRB[4] != 0.0 or aRB[10] != 0.0:
                Dens, _ = solring(X, Y, Z, R, Earth_Mean_Lon[p], aRB)
                Flux += zsrcfunc(Det, Scatt, Therm, prm[170:182],
                                 phase_type=phase_type) * Dens

        if new_iso_comp:                                # isotropic component
            if Src_C is None:
                Src_C = zsrcfunc(Det, Scatt, Therm, prm[33:45], phase_type=phase_type)
            Flux += Src_C * new_isocloud(X, Y, Z, R, prm[14:33])

        f[pix] = np.bincount(own, weights=wts * Flux, minlength=pix.size)

    # First task runs serially: decides the scatt mode and warms any lazy
    # caches (e.g. colcorr table, solar spectrum) before threads start.
    _run(tasks[0])

    rest = tasks[1:]
    if n_jobs == 1 or len(rest) <= 1:
        for t in rest:
            _run(t)
    else:
        with ThreadPoolExecutor(max_workers=n_jobs) as ex:
            for _ in ex.map(_run, rest):
                pass

    zkernel_fast.last_scatt_mode = state['mode']
    return f


zkernel_fast.last_scatt_mode = None