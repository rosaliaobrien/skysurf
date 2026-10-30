import json
import os
import warnings

import numpy as np
from scipy.interpolate import interp1d

from zodisurf.zkernelpy_original import get_parnames
from zodisurf.zkernel_fast import zkernel_fast
from zodisurf.skysurf_params import (get_albedo, get_hong_params, get_mult,
                                     get_emiss, put_zpar)
from zodisurf.solar_sp import solar_sp

try:
    from zodisurf.zkernel_numba import zkernel_numba, NUMBA_AVAILABLE
except Exception as e:
    zkernel_numba, NUMBA_AVAILABLE = None, False
    print(f"[get_zmod] numba engine disabled: {type(e).__name__}: {e}")

__all__ = ['get_zmod', 'read_zpars', 'mk_zdata',
           'set_default_engine', 'get_default_engine']

# ----------------------------------------------------------------------
#  Engine selection
#
#  'numpy' (default) : batched NumPy kernel.
#  'numba'           : compiled scalar kernel. ~4x faster single-core.
#                      Falls back to 'numpy' if numba is unavailable or a
#                      safety guard trips.
# ----------------------------------------------------------------------
_ENGINES = ('numpy', 'numba')
_DEFAULT_ENGINE = 'numpy'
_FALLBACK_WARNED = set()

_PHASE_TYPES = ('kelsall', 'skysurf', 'fitting')

# Slot layout (matches put_zpar with det1=0, and zsrcfunc/scattfunc for 'skysurf')
_SRC_SMOOTH = (33,)                    # smooth cloud (also used by the iso component)
_SRC_BANDS = (59, 85, 111, 137)        # dust bands 1-4
_SRC_RING = (170,)                     # resonant ring + blobs
_SRC_ALL = _SRC_SMOOTH + _SRC_BANDS + _SRC_RING
_HONG = slice(183, 189)                # g1, g2, g3, w1, w2, w3
_ALB_OFF = 0
_EM_OFF = 4

_NUMBA_FIT_OK = set()


def set_default_engine(engine):
    """Set the engine used when get_zmod is called without `engine=`."""
    global _DEFAULT_ENGINE
    if engine not in _ENGINES:
        raise ValueError(f"engine must be one of {_ENGINES}, got {engine!r}")
    _DEFAULT_ENGINE = engine


def get_default_engine():
    return _DEFAULT_ENGINE


def read_zpars(filename=None):
    """
    Reads the ZODI model parameters from a JSON file and returns them in the
    exact order expected by the model. Uses a mapping from model parameter
    names to JSON keys. Any missing parameters are filled with 0.0 for albedo,
    1.0 for emissivity, otherwise 0.0.
    """
    if filename is None:
        filename = os.path.join(os.path.dirname(__file__), 'aend_og_wlabels.json')
    param_map = {
        # Smooth cloud
        'FuncIndx': 'density_flag',
        'PF1_C0': 'C01', 'PF1_C1': 'C11', 'PF1_Ka': 'C21',
        'PF2_C0': 'C02', 'PF2_C1': 'C12', 'PF2_Ka': 'C22',
        'PF3_C0': 'C03', 'PF3_C1': 'C13', 'PF3_Ka': 'C23',
        'To1': 'T', 'To2': 'idk2', 'K': 'idk1', 'Delta': 'T_powerlaw',
        'C_No': 'n0', 'C_Alpha': 'alpha', 'C_Beta': 'beta', 'C_Gamma': 'gamma',
        'C_Mu': 'mu', 'C_Mu2': 'idk4', 'C_Mu3': 'idk5', 'C_Mu4': 'idk6', 'C_Mu5': 'idk7',
        'C_Mu6': 'idk8', 'C_Mu7': 'idk9', 'C_Mu8': 'idk10', 'C_Mu9': 'idk11', 'C_Mu10': 'idk12',
        'C_Omega': 'Omega', 'C_Incl': 'i', 'C_Xo': 'X0', 'C_Yo': 'Y0', 'C_Zo': 'Z0',
        'C_A1': 'A1_sc', 'C_A2': 'A2_sc', 'C_A3': 'A3_sc', 'C_A4': 'A4_sc',
        'C_E3': 'E3_sc', 'C_E4': 'E4_sc', 'C_E5': 'E5_sc', 'C_E6': 'E6_sc', 'C_E7': 'E7_sc',
        'C_E8': 'E8_sc', 'C_E9': 'E9_sc', 'C_E10': 'E10_sc',
        # Band 1
        'B1_No': 'nB1', 'B1_Dz': 'dXB1', 'B1_Dr': 'dRB1', 'B1_Ro': 'R0_B1',
        'B1_Vi': 'nuB1', 'B1_Vr': 'idk15', 'B1_Pi': 'pB1', 'B1_Pr': 'idk17', 'B1_P2r': 'idk18',
        'B1_Omega': 'OmegaB1', 'B1_Incl': 'iB1',
        'B1_Xo': 'X0_B1', 'B1_Yo': 'Y0_B1', 'B1_Zo': 'Z0_B1',
        'B1_A1': 'A1_B1', 'B1_A2': 'A2_B1', 'B1_A3': 'A3_B1', 'B1_A4': 'A4_B1',
        'B1_E3': 'E3_B1', 'B1_E4': 'E4_B1', 'B1_E5': 'E5_B1', 'B1_E6': 'E6_B1', 'B1_E7': 'E7_B1',
        'B1_E8': 'E8_B1', 'B1_E9': 'E9_B1', 'B1_E10': 'E10_B1',
        # Band 2
        'B2_No': 'nB2', 'B2_Dz': 'dXB2', 'B2_Dr': 'dRB2', 'B2_Ro': 'idk23',
        'B2_Vi': 'nuB2', 'B2_Vr': 'idk24', 'B2_Pi': 'idk25', 'B2_Pr': 'idk26', 'B2_P2r': 'idk27',
        'B2_Omega': 'OmegaB2', 'B2_Incl': 'iB2',
        'B2_Xo': 'idk28', 'B2_Yo': 'idk29', 'B2_Zo': 'idk30',
        'B2_A1': 'A1_B2', 'B2_A2': 'A2_B2', 'B2_A3': 'A3_B2', 'B2_A4': 'A4_B2',
        'B2_E3': 'E3_B2', 'B2_E4': 'E4_B2', 'B2_E5': 'E5_B2', 'B2_E6': 'E6_B2', 'B2_E7': 'E7_B2',
        'B2_E8': 'E8_B2', 'B2_E9': 'E9_B2', 'B2_E10': 'E10_B2',
        # Band 3 (uses ignored_BZ parameters)
        'B3_No': 'ignored_BZ0', 'B3_Dz': 'ignored_BZ1', 'B3_Dr': 'ignored_BZ2',
        'B3_Ro': 'ignored_BZ3',
        'B3_Vi': 'ignored_BZ4', 'B3_Vr': 'ignored_BZ5', 'B3_Pi': 'ignored_BZ6',
        'B3_Pr': 'ignored_BZ7', 'B3_P2r': 'ignored_BZ8',
        'B3_Omega': 'ignored_BZ9', 'B3_Incl': 'ignored_BZ10',
        'B3_Xo': 'ignored_BZ11', 'B3_Yo': 'ignored_BZ12', 'B3_Zo': 'ignored_BZ13',
        'B3_A1': 'ignored_BZ14', 'B3_A2': 'ignored_BZ15', 'B3_A3': 'ignored_BZ16',
        'B3_A4': 'ignored_BZ17',
        'B3_E3': 'ignored_BZ18', 'B3_E4': 'ignored_BZ19', 'B3_E5': 'ignored_BZ20',
        'B3_E6': 'ignored_BZ21',
        'B3_E7': 'ignored_BZ22', 'B3_E8': 'ignored_BZ23', 'B3_E9': 'ignored_BZ24',
        'B3_E10': 'ignored_BZ25',
        # Band 4 (uses nB3/dXB3 parameters)
        'B4_No': 'nB3', 'B4_Dz': 'dXB3', 'B4_Dr': 'dRB3', 'B4_Ro': 'idk33',
        'B4_Vi': 'nuB3', 'B4_Vr': 'idk34', 'B4_Pi': 'idk35', 'B4_Pr': 'idk36', 'B4_P2r': 'idk37',
        'B4_Omega': 'OmegaB3', 'B4_Incl': 'iB3',
        'B4_Xo': 'idk38', 'B4_Yo': 'idk39', 'B4_Zo': 'idk40',
        'B4_A1': 'A1_B3', 'B4_A2': 'A2_B3', 'B4_A3': 'A3_B3', 'B4_A4': 'A4_B3',
        'B4_E3': 'E3_B3', 'B4_E4': 'E4_B3', 'B4_E5': 'E5_B3', 'B4_E6': 'E6_B3', 'B4_E7': 'E7_B3',
        'B4_E8': 'E8_B3', 'B4_E9': 'E9_B3', 'B4_E10': 'E10_B3',
        # SR/Resonant Ring
        'SR_No': 'nSR', 'SR_R': 'R_SR', 'SR_dR': 'sigma_rSR', 'SR_dZ': 'sigma_zSR',
        # Leading Blob
        'LB_No': 'leading_blob_n', 'LB_R': 'leading_blob_R', 'LB_dR': 'leading_blob_dR',
        'LB_Theta': 'leading_blob_theta', 'LB_dTheta': 'leading_blob_dtheta',
        'LB_dZ': 'leading_blob_dZ',
        # Trailing Blob
        'TB_No': 'trailing_blob_n', 'TB_R': 'R_TB', 'TB_dR': 'sigma_rTB',
        'TB_Theta': 'theta_TB', 'TB_dTheta': 'sigma_thetaTB', 'TB_dZ': 'sigma_zTB',
        # RB Geometry
        'RB_Omega': 'OmegaRB', 'RB_Incl': 'iRB', 'RB_Xo': 'blob_X0',
        'RB_Yo': 'blob_Y0', 'RB_Zo': 'blob_Z0',
        # RB Albedo and Emissivity
        'RB_A1': 'A1_SR', 'RB_A2': 'A2_SR', 'RB_A3': 'A3_SR', 'RB_A4': 'A4_SR',
        'RB_E3': 'E3_SR', 'RB_E4': 'E4_SR', 'RB_E5': 'E5_SR', 'RB_E6': 'E6_SR',
        'RB_E7': 'E7_SR', 'RB_E8': 'E8_SR', 'RB_E9': 'E9_SR', 'RB_E10': 'E10_SR'
    }
    try:
        with open(filename, 'r') as f:
            data = json.load(f)
        parnames = get_parnames()
        zpar = []
        for name in parnames:
            json_key = param_map.get(name, name)
            if name.startswith(('C_A', 'B1_A', 'B2_A', 'B3_A', 'SR_A', 'B4_A', 'RB_A')):
                val = data.get(json_key, 0.0)
            elif name.startswith(('C_E', 'B1_E', 'B2_E', 'B3_E', 'SR_E', 'B4_E', 'RB_E')):
                val = data.get(json_key, 1.0)
            else:
                val = data.get(json_key, 0.0)
            zpar.append(val)
        return np.array(zpar, dtype=float)
    except Exception as e:
        print(f"Error reading data from '{filename}': {e}")
        return None


# Default parameters are read from disk once and then reused
_DEFAULT_ZPAR = None


def _default_zpar():
    global _DEFAULT_ZPAR
    if _DEFAULT_ZPAR is None:
        _DEFAULT_ZPAR = read_zpars()
        if _DEFAULT_ZPAR is None:
            raise RuntimeError("Could not read default zodi parameters "
                               "(aend_og_wlabels.json).")
    return _DEFAULT_ZPAR.copy()        # copy: protects the cache from in-place edits


def mk_zdata(lambda_, day, lon, lat):
    """
    Generate a NumPy structured array for use in zodi model calculations.
    """
    lambda_ = np.atleast_1d(lambda_)
    day = np.atleast_1d(day)
    lon = np.atleast_1d(lon)
    lat = np.atleast_1d(lat)

    npts = max(len(lambda_), len(day), len(lon), len(lat))

    for arr, name in zip([lambda_, day, lon, lat], ['lambda', 'day', 'lon', 'lat']):
        if len(arr) != 1 and len(arr) != npts:
            raise ValueError(f"{name} must be a scalar or match the length of "
                             f"the longest input ({npts}).")

    lambda_ = np.resize(lambda_, npts)
    day = np.resize(day, npts)
    lon = np.resize(lon, npts)
    lat = np.resize(lat, npts)

    dtype = np.dtype([
        ('pixel_no', 'i4'),
        ('day1990', 'f8'),
        ('wave_len', 'f8'),
        ('longitude', 'f8'),
        ('latitude', 'f8'),
        ('photometry', 'f8')
    ])

    zdata = np.zeros(npts, dtype=dtype)
    zdata['day1990'] = day
    zdata['wave_len'] = lambda_
    zdata['longitude'] = lon
    zdata['latitude'] = lat

    return zdata


def _broadcast_fitting(npts, **fits):
    """Scalar, (npts,) for albedo/emissivity; (6,) or (npts, 6) for Hong."""
    out = {}
    for name, val in fits.items():
        if val is None:
            continue
        v = np.asarray(val, dtype=float)
        if name == 'hong':
            if v.shape == (6,):
                v = np.broadcast_to(v, (npts, 6))
            if v.shape != (npts, 6):
                raise ValueError(f"fitting_hong must have shape (6,) or "
                                 f"({npts}, 6); got {v.shape}")
        else:
            if v.size == 1:
                v = np.full(npts, float(v.ravel()[0]))
            if v.shape != (npts,):
                raise ValueError(f"fitting_{name} must be a scalar or shape "
                                 f"({npts},); got {v.shape}")
        if not np.all(np.isfinite(v)):
            raise ValueError(f"fitting_{name} contains non-finite values")
        out[name] = v
    return out


def _apply_fitting(zpar_u, inv, fit):
    """One parameter row per unique (wavelength, fitting values). Returns (rows, group)."""
    npts = inv.size
    cols = [inv.reshape(npts, 1).astype(float)]
    cols += [v.reshape(npts, -1) for v in fit.values()]
    _, first, new_inv = np.unique(np.hstack(cols), axis=0,
                                  return_index=True, return_inverse=True)

    rows = zpar_u[inv[first]].copy()
    for r, p in enumerate(first):
        if 'hong' in fit:
            rows[r, _HONG] = fit['hong'][p]
        if 'albedo' in fit:
            for s0 in _SRC_ALL:
                rows[r, s0 + _ALB_OFF] = fit['albedo'][p]
        for key, blocks in (('em_smooth', _SRC_SMOOTH),
                            ('em_band', _SRC_BANDS),
                            ('em_ring', _SRC_RING)):
            if key in fit:
                for s0 in blocks:
                    rows[r, s0 + _EM_OFF] = fit[key][p]
    return rows, np.asarray(new_inv).ravel()


def get_zmod(lambda_, phase_type, day, lon, lat, L2=False, zpar=None,
             solar_irr=None, new_iso_comp=False, iso_comp_only=False,
             engine=None, n_jobs=1, scatt_mode='auto',
             fitting_hong=None, fitting_albedo=None, fitting_em_band=None,
             fitting_em_smooth=None, fitting_em_ring=None):
    """
    Compute the Kelsall et al. 1998 (ApJ,508,44) ZODI model intensity
    for a given set of LOS, Wavelength, and Time.

    Parameters:
    -----------
    lambda_ : float or array
        Wavelength in microns
    phase_type : str
        "kelsall" - standard Kelsall+1998 phase function
        "skysurf" - O'Brien+2026 phase function / albedo / emissivity
        "fitting" - To be used when performing model fitting.
                    Utilizes SKYSURF physics (Hong phase function), with every
                    scattering/emission value supplied by the caller:
                      * ALL fitting_* arguments and solar_irr are REQUIRED
                      * get_albedo, get_hong_params, get_mult, get_emiss and
                        solar_sp are NOT used
                      * the 3.5 micron limit is not enforced
                      * the output is NOT rounded
                    Cloud geometry and thermal parameters (To1, Delta) still
                    come from zpar.
    day : array_like
        1990 day number(s), where 1.0 = 1 Jan 1990
    lon : array_like
        Ecliptic longitude(s) in degrees
    lat : array_like
        Ecliptic latitude(s) in degrees
    L2 : bool
        If True, then put observer at Earth-Sun L2 (e.g., JWST)
    zpar : array_like, optional
        Array of zodi model params. If not specified, defaults are loaded
    solar_irr : float or array (required for 'fitting')
        Solar irradiance (in MJy/sr) corresponding to bandpass. In 'fitting'
        mode a single value is applied to every pixel.
    new_iso_comp : bool, optional
        If True, include isotropic component (default: False)
    iso_comp_only : bool, optional
        If True, model ONLY isotropic component (default: False)
    engine : {'numpy', 'numba'}, optional
        Compute backend. Default is the module default (see
        set_default_engine), initially 'numpy'.
          'numpy' - batched NumPy; bit-identical to the original code
          'numba' - compiled scalar loops; ~1e-13 relative difference,
                    ~4x faster on a single core. 
    n_jobs : int, optional
        Threads to use (default: 1 core). Set to 1 inside MCMC/MPI pools.
        Set to -1 to use all available cores.
    scatt_mode : str, optional
        'auto' (default), 'vector' or 'loop'; see zkernel_fast

    Fitting parameters (REQUIRED when phase_type='fitting', not allowed otherwise):
    -------------------------------------------------------------------------------
    fitting_hong : array_like
        Hong parameters [g1, g2, g3, w1, w2, w3]; written to zpar[183:189].
        Shape (6,) for all pixels, or (npts, 6) per pixel.
    fitting_albedo : float or array_like
        Albedo applied to every component. Scalar or shape (npts,).
    fitting_em_smooth : float or array_like
        Emissivity of the smooth cloud (also used by the isotropic
        component). Scalar or shape (npts,).
    fitting_em_band : float or array_like
        Emissivity of dust bands 1-4. Scalar or shape (npts,).
    fitting_em_ring : float or array_like
        Emissivity of the resonant ring + blobs. Scalar or shape (npts,).

    Returns:
    --------
    numpy.ndarray
        Zodi model intensity in MJy sr^-1 (rounded to 5 decimals, except
        when phase_type='fitting')
    """

    # Use all cores if n_jobs == -1
    # Use just 1 core if n_jobs is None (but it should never be None)
    if n_jobs == -1 or n_jobs is None:
        n_jobs = os.cpu_count() or 1

    lambda_ = np.atleast_1d(np.asarray(lambda_, dtype=float))  # wavelengths
    day = np.atleast_1d(np.asarray(day, dtype=float))
    lon = np.atleast_1d(np.asarray(lon, dtype=float))
    lat = np.atleast_1d(np.asarray(lat, dtype=float))

    # The wavelength array may be a single number, so you can compute many
    # pointings at the same wavelength
    if lambda_.size == 1 and lambda_.shape != day.shape:
        lambda_ = np.full_like(day, lambda_[0], dtype=float)

    # ---- phase type ---------------------------------------------------
    if phase_type not in _PHASE_TYPES:
        raise ValueError(f"phase_type must be one of {_PHASE_TYPES}, got {phase_type!r}")

    fit_mode = (phase_type == 'fitting')
    kernel_phase = 'skysurf' if fit_mode else phase_type   # what the kernels see

    # Keys match the argument names (fitting_<key>), so error messages
    # can name exactly what is missing.
    fit_args = dict(hong=fitting_hong, albedo=fitting_albedo,
                    em_smooth=fitting_em_smooth, em_band=fitting_em_band,
                    em_ring=fitting_em_ring)

    if fit_mode:
        missing = [f"fitting_{k}" for k, v in fit_args.items() if v is None]
        if solar_irr is None:
            # Without solar_irr the kernels would call solar_sp, which
            # 'fitting' must not use.
            missing.append("solar_irr")
        if missing:
            raise ValueError("phase_type='fitting' requires every fitting value; "
                             "missing: " + ", ".join(missing))
    
    # If any of the fitting params are defined, then phase_type must be "fittin"
    elif any(v is not None for v in fit_args.values()):
        raise ValueError("fitting_* arguments require phase_type='fitting'")

    # ---- solar irradiance ---------------------------------------------
    if solar_irr is not None:
        solar_irr = np.atleast_1d(np.asarray(solar_irr, dtype=float))
        if fit_mode and solar_irr.size == 1 and lambda_.size > 1:
            solar_irr = np.full(lambda_.shape, solar_irr[0])

    # Base shapes that must always match
    if not (lambda_.shape == day.shape == lon.shape == lat.shape):
        raise ValueError(
            f"Shape mismatch: lambda_ {lambda_.shape} must match "
            f"day {day.shape}, lon {lon.shape}, and lat {lat.shape}."
        )

    # Include solar_irr in the shape check only if it's provided
    if solar_irr is not None and solar_irr.shape != lambda_.shape:
        raise ValueError(
            f"Shape mismatch: solar_irr {solar_irr.shape} must match lambda_ {lambda_.shape}."
        )
    if fit_mode and not np.all(np.isfinite(solar_irr)):
        raise ValueError("solar_irr contains non-finite values")

    # At 3.5 micron, thermal emission starts to kick in, and HST really couldn't
    # calibrate that. So ignore these wavelengths -- except when fitting.
    if phase_type == 'skysurf' and np.any(lambda_ > 3.5):
        raise ValueError("This code does not work at lambda > 3.5 micron.")
    if np.any(lon < 0) or np.any(lon > 360):
        raise ValueError("Longitude (lon) must be in the range 0 to 360 degrees.")
    if np.any(lat < -90) or np.any(lat > 90):
        raise ValueError("Latitude (lat) must be in the range -90 to 90 degrees.")

    # ---- engine resolution -------------------------------------------
    if engine is None:
        engine = _DEFAULT_ENGINE
    if engine not in _ENGINES:
        raise ValueError(f"engine must be one of {_ENGINES}, got {engine!r}")

    # ---- parameters ---------------------------------------------------
    if zpar is None:
        zpar = _default_zpar()
    else:
        zpar = np.array(zpar, dtype=float)  # custom input parameters

    # ---- build parameters ONCE per unique wavelength -------------------
    lam_u, inv = np.unique(lambda_, return_inverse=True)
    inv = np.asarray(inv).ravel()

    if phase_type in ('kelsall', 'fitting'):
        zpar_u = np.tile(zpar, (len(lam_u), 1)) # vectorized to have zpar list for all wavelengths or pointings

    elif phase_type == 'skysurf':
        albedo_arr = get_albedo(lam_u)
        hong_params_arr = get_hong_params(lam_u)
        mult_arr = get_mult(lam_u)
        hong_params_arr[:, -3:] = hong_params_arr[:, 3:] * mult_arr
        emiss_arr = get_emiss(lam_u)
        zpar_u = put_zpar(zpar.copy(), 0, 0, 0, albedo_arr, det1=0,
                          hg3=hong_params_arr, E1=emiss_arr)

    zpar_u = np.atleast_2d(np.asarray(zpar_u, dtype=float))
    if zpar_u.shape[0] != len(lam_u):
        raise RuntimeError(f"Expected {len(lam_u)} parameter rows, "
                           f"got {zpar_u.shape[0]}.")

    if fit_mode:
        fit = _broadcast_fitting(lambda_.size, **fit_args)
        zpar_u, inv = _apply_fitting(zpar_u, inv, fit)
        zpar_u[:, 1:10] = 0.0      # unused by the skysurf kernels; matches put_zpar

    data = mk_zdata(lambda_, day, lon, lat)

    kw = dict(phase_type=kernel_phase, L2=L2, no_colcorr=True,
              solar_irr=solar_irr, new_iso_comp=new_iso_comp,
              iso_comp_only=iso_comp_only)

    # Don't round when fitting cause you might want high precision
    def _finish(zodi):
        return zodi if fit_mode else np.round(zodi, 5)

    # ---- 'numba' engine (no fallback: fail loudly) --------------------
    if engine == 'numba':
        if zkernel_numba is None or not NUMBA_AVAILABLE:
            raise RuntimeError(
                "engine='numba' was requested but numba is not available. "
                "Install numba, or call with engine='numpy'.")
        zodi = zkernel_numba(data, zpar_u, group=inv, n_jobs=n_jobs, **kw)
        if fit_mode:
            # Records that numba actually ran (checked by the tests)
            _NUMBA_FIT_OK.add((L2, new_iso_comp, iso_comp_only,
                               tuple(lam_u.tolist())))
        return _finish(zodi)

    # ---- 'numpy' engine ------------------------------------------------
    zodi = zkernel_fast(data, zpar_u, group=inv, n_jobs=n_jobs,
                        scatt_mode=scatt_mode, **kw)
    return _finish(zodi)