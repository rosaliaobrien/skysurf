"""
zkernel_numba.py
----------------
Numba-compiled scalar line-of-sight integrator for ZodiSURF.

Instead of building 40M-element temporary arrays (as the NumPy kernel does),
this walks one point at a time in compiled code, keeping the CPU cache hot.
Memory traffic drops from gigabytes to kilobytes.

SAFETY DESIGN
  * Guards: refuses to run (caller falls back to zkernel_fast) when the
    reconstructed physics cannot be guaranteed -- see _guard_reason().
  * Source-function coefficients are extracted EXACTLY from the user's own
    zsrcfunc by linear probing, so no assumption is made about albedo/
    emissivity selection.
  * First use of any parameter set is verified against zkernel_fast.

Accuracy: ~1e-13 relative vs zkernel_fast (summation order differs).
"""
import os
import numpy as np

from zodisurf.zkernelpy_original import zsrcfunc, scattfunc, thermfunc, solar_sp

try:
    from numba import njit, prange, set_num_threads, get_num_threads
    NUMBA_AVAILABLE = True
except Exception as _numba_err:          # includes NumPy-version rejection
    NUMBA_IMPORT_ERROR = _numba_err                                   # graceful degradation
    NUMBA_AVAILABLE = False
    prange = range
    def njit(*args, **kwargs):
        if len(args) == 1 and callable(args[0]):
            return args[0]
        return lambda f: f
    def set_num_threads(n): pass
    def get_num_threads(): return 1

__all__ = ["zkernel_numba", "NUMBA_AVAILABLE", "guard_reason", "warmup"]

# ---------------------------------------------------------------- constants
RMAX   = 5.2
_SIMP_H = 0.025
H32     = np.float32(_SIMP_H)
W_END   = np.float64(np.float32(0.333333) * H32)        # identical to simpint
W_ODD   = np.float64(np.float32(1.333333) * H32)
W_EVEN  = np.float64(np.float32(0.666667) * H32)

C_LIGHT = 2.9979e14
H_PL    = 6.6262e-34
K_B     = 1.3807e-23
SIGMA   = 5.6697e-12
COEFF   = 15.0 / np.pi**5 * SIGMA
D2R     = np.pi / 180.0
TWOPI   = 2.0 * np.pi

DBWAVE = np.array([1.25, 2.2, 3.5, 4.9, 12., 25., 60., 100., 140., 240.])
SOLFLUX_1AU = np.array([2.3405606e+08, 1.2309874e+08, 64292872., 35733824.,
                        5763843.0, 1327989.4, 230553.73, 82999.336,
                        42346.605, 14409.608])

# --------------------------------------------- packed parameter layout (P)
# Using a flat float64 array keeps the Numba signature small and stable.
P_SIZE = 128
# smooth cloud
C_NO, C_ALPHA, C_BETA, C_GAMMA, C_MU = 0, 1, 2, 3, 4
C_MU3, C_MU4, C_MU5, C_MU6, C_MU7, C_MU8 = 5, 6, 7, 8, 9, 10
C_SINO, C_COSO, C_SINI, C_COSI, C_XO, C_YO, C_ZO = 11, 12, 13, 14, 15, 16, 17
C_FUNCIDX, C_A, C_E = 18, 19, 20
# bands: 16 slots each, b = 0..3
B_BASE, B_STRIDE = 21, 16
B_NO, B_DZ, B_DR, B_RO, B_VI, B_PI, B_PR = 0, 1, 2, 3, 4, 5, 6
B_SINO, B_COSO, B_SINI, B_COSI, B_XO, B_YO, B_ZO, B_A, B_E = 7, 8, 9, 10, 11, 12, 13, 14, 15
# ring + blobs
R_BASE = 85
R_SRNO, R_SRR, R_SRDR, R_SRDZ = 0, 1, 2, 3
R_LBNO, R_LBR, R_LBDR, R_LBTH, R_LBDTH, R_LBDZ = 4, 5, 6, 7, 8, 9
R_TBNO, R_TBR, R_TBDR, R_TBTH, R_TBDTH, R_TBDZ = 10, 11, 12, 13, 14, 15
R_SINO, R_COSO, R_SINI, R_COSI, R_A, R_E = 16, 17, 18, 19, 20, 21
# thermal / scattering / modes
T_TO1, T_DELTA, T_CFACT, T_LAM = 107, 108, 109, 110
S_DO, S_FLUX, S_ISHONG = 111, 112, 113
S_P0 = 114                       # kelsall: C0,C1,ka,v | hong: g1,g2,g3,w1,w2,w3
M_ISO, M_ISOONLY = 120, 121


# ======================================================================
#  Compiled scalar physics
# ======================================================================
@njit(cache=True, inline='always')
def _zplanck(T, lam):
    """Planck function, W/cm^2/sr. Scalar form of zPlanck()."""
    chi = H_PL * C_LIGHT / (lam * K_B * T)
    T4 = T * T * T * T
    if chi < 0.001:
        return COEFF * T4 * chi * chi * chi
    elif chi <= 50.0:
        c4 = chi * chi * chi * chi
        return COEFF * T4 * c4 / (np.exp(chi) - 1.0)
    else:
        c4 = chi * chi * chi * chi
        return COEFF * T4 * c4 * np.exp(-chi)


@njit(cache=True, inline='always')
def _therm(R, P):
    """Thermal emission. Only valid when K == 0 (enforced by guard)."""
    T1 = P[T_TO1] / R**P[T_DELTA]
    return _zplanck(T1, P[T_LAM]) * P[T_CFACT]


@njit(cache=True, inline='always')
def _scatt(los, R, Re, sin_se, cos_se, P):
    """Scattered sunlight. Scalar form of scattfunc()."""
    if P[S_DO] == 0.0:
        return 0.0
    solflux = P[S_FLUX] / (R * R)

    s = Re / R * sin_se
    if s > 1.0:
        s = 1.0
    elif s < -1.0:
        s = -1.0
    pa = np.arcsin(s)
    if los >= Re * cos_se:
        x = np.pi - pa
    else:
        x = pa

    if P[S_ISHONG] != 0.0:
        cx = np.cos(x)
        phi = 0.0
        for i in range(3):
            g = P[S_P0 + i]
            w = P[S_P0 + 3 + i]
            phi += w * (1.0 - g * g) / (1.0 + g * g - 2.0 * g * cx)**1.5
    else:
        phi = (P[S_P0] + P[S_P0 + 1] * x + np.exp(P[S_P0 + 2] * x)) / P[S_P0 + 3]
    return solflux * phi


@njit(cache=True, inline='always')
def _zcloud(x, y, z, R, P):
    """Smooth-cloud density. Scalar form of zcloud()."""
    No = P[C_NO]
    if No == 0.0:
        return 0.0
    Xp = x - P[C_XO]
    Yp = y - P[C_YO]
    Zp = z - P[C_ZO]
    Rc = np.sqrt(Xp * Xp + Yp * Yp + Zp * Zp)
    Zc = P[C_SINO] * P[C_SINI] * Xp - P[C_COSO] * P[C_SINI] * Yp + P[C_COSI] * Zp
    Zeta = abs(Zc / Rc)

    dR = R - 1.0
    m6, m7, m8 = P[C_MU6], P[C_MU7], P[C_MU8]
    AlphaR = P[C_ALPHA] + (m6 * m6 * dR + m7 * m7 * dR * dR + m8 * m8 * dR * dR * dR)
    if AlphaR > 50.0:
        AlphaR = 50.0

    Beta, Gamma, Mu = P[C_BETA], P[C_GAMMA], P[C_MU]
    fi = int(P[C_FUNCIDX])

    if fi == 0:                                        # John Good
        d2 = Rc * Rc - Zc * Zc
        if d2 < 1e-10:
            d2 = 1e-10
        Dc = np.sqrt(d2)
        vert = np.exp(-Beta * abs(Zc / Dc)**Gamma)
        rad = 1.0 / Dc**AlphaR
    elif fi == 1:                                      # Modified Fan
        vert = np.exp(-Beta * Zeta**Gamma)
        rad = 1.0 / Rc**AlphaR
    elif fi == 2:                                      # Widened Modified Fan
        if Zeta < Mu:
            GZR = 0.5 * Zeta * Zeta / Mu
        else:
            GZR = Zeta - 0.5 * Mu
        vert = np.exp(-Beta * GZR**Gamma)
        rad = 1.0 / Rc**AlphaR
    elif fi == 3:                                      # Ellipsoid
        bt = 1.0 + (Beta * Zeta)**2
        vert = 1.0 / bt**Gamma
        rad = 1.0 / Rc**AlphaR
    elif fi == 4:                                      # Sombrero
        c2 = 1.0 - Zeta * Zeta
        if c2 < 1e-10:
            c2 = 1e-10
        cosB = np.sqrt(c2)
        vert = (1.0 + Beta * cosB**Gamma) / (1.0 + Beta)
        rad = 1.0 / Rc**AlphaR
    else:                                              # 5: polynomial z-height
        aZ = abs(z)
        BetaZ = P[C_BETA] + P[C_MU3] * aZ + P[C_MU4] * aZ * aZ + P[C_MU5] * aZ * aZ * aZ
        vert = np.exp(-BetaZ * Zeta)
        rad = 1.0 / Rc**AlphaR

    return No * vert * rad


@njit(cache=True, inline='always')
def _isocloud(R, P):
    """Isotropic component. Scalar form of new_isocloud()."""
    dR = R - 1.0
    m6, m7, m8 = P[C_MU6], P[C_MU7], P[C_MU8]
    AlphaR = -1.46 + (m6 * m6 * dR + m7 * m7 * dR * dR + m8 * m8 * dR * dR * dR)
    if AlphaR > 50.0:
        AlphaR = 50.0
    return 1.356e-9 / R**AlphaR


@njit(cache=True, inline='always')
def _migband(x, y, z, P, o):
    """Dust-band density. Scalar form of migband(). `o` = band base offset."""
    No = P[o + B_NO]
    if No == 0.0:
        return 0.0
    Xp = x - P[o + B_XO]
    Yp = y - P[o + B_YO]
    Zp = z - P[o + B_ZO]
    Zc = P[o + B_SINO] * P[o + B_SINI] * Xp - P[o + B_COSO] * P[o + B_SINI] * Yp \
        + P[o + B_COSI] * Zp
    Rc = np.sqrt(Xp * Xp + Yp * Yp + Zp * Zp)

    ZDz = abs(Zc / Rc) / P[o + B_DZ]
    RDr = Rc / P[o + B_DR]
    ViTerm = 1.0 + (ZDz**P[o + B_PI]) / P[o + B_VI]

    z2 = ZDz * ZDz
    arg2 = z2 * z2 * z2                                # ZDz**6
    if arg2 > 86.0:
        return 0.0
    r2 = RDr * RDr
    r4 = r2 * r2
    r8 = r4 * r4
    arg1 = r8 * r8 * r4                                # RDr**20

    base = No * (P[o + B_RO] / Rc)**P[o + B_PR] * np.exp(-arg2) * ViTerm
    if arg1 <= 86.0:
        return base * (1.0 - np.exp(-arg1))
    return base


@njit(cache=True, inline='always')
def _wrap_pi(d):
    """Replicates  (d % 2pi) then shift into (-pi, pi]."""
    d = d - TWOPI * np.floor(d / TWOPI)                # Python-style mod -> [0, 2pi)
    if d > np.pi:
        d -= TWOPI
    return d


@njit(cache=True, inline='always')
def _solring(x, y, z, R, Theta, P):
    """Resonant ring + leading/trailing blobs. Scalar form of solring()."""
    o = R_BASE
    srn, lbn, tbn = P[o + R_SRNO], P[o + R_LBNO], P[o + R_TBNO]
    if srn == 0.0 and lbn == 0.0 and tbn == 0.0:
        return 0.0

    SR_Z = P[o + R_SINO] * P[o + R_SINI] * x - P[o + R_COSO] * P[o + R_SINI] * y \
        + P[o + R_COSI] * z
    dTheta = np.arctan2(y, x) - Theta
    dens = 0.0

    if srn != 0.0:
        dr = R - P[o + R_SRR]
        arg = dr * dr / (P[o + R_SRDR]**2) + abs(SR_Z) / P[o + R_SRDZ]
        if arg <= 86.0:
            dens += srn * np.exp(-arg)
    if lbn != 0.0:
        d = _wrap_pi(dTheta - P[o + R_LBTH])
        dr = R - P[o + R_LBR]
        arg = dr * dr / (P[o + R_LBDR]**2) + d * d / (P[o + R_LBDTH]**2) \
            + abs(SR_Z) / P[o + R_LBDZ]
        if arg <= 86.0:
            dens += lbn * np.exp(-arg)
    if tbn != 0.0:
        d = _wrap_pi(dTheta - P[o + R_TBTH])
        dr = R - P[o + R_TBR]
        arg = dr * dr / (P[o + R_TBDR]**2) + d * d / (P[o + R_TBDTH]**2) \
            + abs(SR_Z) / P[o + R_TBDZ]
        if arg <= 86.0:
            dens += tbn * np.exp(-arg)
    return dens


# ======================================================================
#  Compiled driver: one thread per block of pixels
# ======================================================================
@njit(parallel=True, cache=True)
def _core(P, nsteps, X0, Y0, Re, COSLAT, SINLAT, COSLON, SINLON,
          SIN_SE, COS_SE, EML, out):
    npts = nsteps.shape[0]
    iso_only = P[M_ISOONLY] != 0.0
    do_iso = P[M_ISO] != 0.0

    for p in prange(npts):
        n = nsteps[p]
        re = Re[p]
        x0 = X0[p]
        y0 = Y0[p]
        cl, sl = COSLAT[p], SINLAT[p]
        clo, slo = COSLON[p], SINLON[p]
        sse, cse = SIN_SE[p], COS_SE[p]
        eml = EML[p]
        acc = 0.0

        for k in range(n):
            los = np.float64(np.float32(k) * H32)
            if k == 0 or k == n - 1:
                w = W_END
            elif k % 2 == 1:
                w = W_ODD
            else:
                w = W_EVEN

            sxy = los * cl
            x = x0 + sxy * clo
            y = y0 + sxy * slo
            z = los * sl
            r = np.sqrt(x * x + y * y + z * z)

            sc = _scatt(los, r, re, sse, cse, P)
            th = _therm(r, P)
            flux = 0.0

            if not iso_only:
                d = _zcloud(x, y, z, r, P)
                if d != 0.0:
                    flux += (P[C_A] * sc + P[C_E] * th) * d
                for b in range(4):
                    o = B_BASE + b * B_STRIDE
                    d = _migband(x, y, z, P, o)
                    if d != 0.0:
                        flux += (P[o + B_A] * sc + P[o + B_E] * th) * d
                d = _solring(x, y, z, r, eml, P)
                if d != 0.0:
                    flux += (P[R_BASE + R_A] * sc + P[R_BASE + R_E] * th) * d

            if do_iso:
                flux += (P[C_A] * sc + P[C_E] * th) * _isocloud(r, P)

            acc += w * flux
        out[p] = acc


# ======================================================================
#  Parameter packing (plain Python, runs once per group)
# ======================================================================
def _src_coeffs(det, a_src, phase_type):
    """
    Exact albedo/emissivity extraction by linear probing.
    Source = Albedo*Scatt + Emiss*(1-Albedo)*Therm  is linear in (Scatt, Therm),
    so two probes recover both coefficients without knowing zsrcfunc internals.
    """
    one = np.ones(1)
    zero = np.zeros(1)
    A = float(np.asarray(zsrcfunc(det, one, zero, a_src, phase_type=phase_type)).ravel()[0])
    E = float(np.asarray(zsrcfunc(det, zero, one, a_src, phase_type=phase_type)).ravel()[0])
    return A, E


def _pack(prm, Lambda, det, phase_type, sirr, new_iso_comp, iso_comp_only):
    P = np.zeros(P_SIZE, dtype=np.float64)

    # --- smooth cloud ---
    c = prm[14:33]
    P[C_NO], P[C_ALPHA], P[C_BETA], P[C_GAMMA], P[C_MU] = c[0], c[1], c[2], c[3], c[4]
    P[C_MU3], P[C_MU4], P[C_MU5] = c[6], c[7], c[8]
    P[C_MU6], P[C_MU7], P[C_MU8] = c[9], c[10], c[11]
    om, inc = c[14] * D2R, c[15] * D2R
    P[C_SINO], P[C_COSO] = np.sin(om), np.cos(om)
    P[C_SINI], P[C_COSI] = np.sin(inc), np.cos(inc)
    P[C_XO], P[C_YO], P[C_ZO] = c[16], c[17], c[18]
    P[C_FUNCIDX] = prm[0]
    P[C_A], P[C_E] = _src_coeffs(det, prm[33:45], phase_type)

    # --- dust bands 1-4 ---
    for b, (d0, s0) in enumerate(((45, 59), (71, 85), (97, 111), (123, 137))):
        o = B_BASE + b * B_STRIDE
        g = prm[d0:d0 + 14]
        P[o + B_NO] = g[0]
        P[o + B_DZ] = g[1] * D2R
        P[o + B_DR], P[o + B_RO], P[o + B_VI] = g[2], g[3], g[4]
        P[o + B_PI], P[o + B_PR] = g[6], g[7]
        om, inc = g[9] * D2R, g[10] * D2R
        P[o + B_SINO], P[o + B_COSO] = np.sin(om), np.cos(om)
        P[o + B_SINI], P[o + B_COSI] = np.sin(inc), np.cos(inc)
        P[o + B_XO], P[o + B_YO], P[o + B_ZO] = g[11], g[12], g[13]
        P[o + B_A], P[o + B_E] = _src_coeffs(det, prm[s0:s0 + 12], phase_type)

    # --- resonant ring + blobs ---
    o = R_BASE
    g = prm[np.r_[149:170, 182]]
    for i in range(16):
        P[o + i] = g[i]
    P[o + R_LBTH] *= D2R
    P[o + R_LBDTH] *= D2R
    P[o + R_TBTH] *= D2R
    P[o + R_TBDTH] *= D2R
    om, inc = g[16] * D2R, g[17] * D2R
    P[o + R_SINO], P[o + R_COSO] = np.sin(om), np.cos(om)
    P[o + R_SINI], P[o + R_COSI] = np.sin(inc), np.cos(inc)
    P[o + R_A], P[o + R_E] = _src_coeffs(det, prm[170:182], phase_type)

    # --- thermal ---
    P[T_TO1], P[T_DELTA] = prm[10], prm[13]
    P[T_CFACT] = Lambda / 3.0e-10
    P[T_LAM] = Lambda

    # --- scattering ---
    P[S_DO] = 1.0 if det <= 2 else 0.0
    if sirr is not None:
        P[S_FLUX] = float(sirr)
    elif phase_type == 'skysurf':
        P[S_FLUX] = float(np.asarray(solar_sp(np.atleast_1d(Lambda))).ravel()[0])
    elif det >= 0:
        P[S_FLUX] = SOLFLUX_1AU[det]
    else:
        P[S_FLUX] = 0.0

    if phase_type == 'skysurf':
        P[S_ISHONG] = 1.0
        P[S_P0:S_P0 + 6] = prm[183:189]
    else:
        P[S_ISHONG] = 0.0
        C0, C1, ka = prm[1 + 3 * det:1 + 3 * det + 3] if det >= 0 else (0.0, 0.0, 0.0)
        P[S_P0], P[S_P0 + 1], P[S_P0 + 2] = C0, C1, ka
        qk = (np.exp(ka * np.pi) + 1.0) / (ka * ka + 1.0)
        P[S_P0 + 3] = (2.0 * np.pi) * (C0 * 2.0 + C1 * np.pi + qk)

    P[M_ISO] = 1.0 if new_iso_comp else 0.0
    P[M_ISOONLY] = 1.0 if iso_comp_only else 0.0
    return P


# ======================================================================
#  Guards + verification
# ======================================================================
def guard_reason(a, phase_type, no_colcorr):
    """Returns None if Numba may be used, else a human-readable reason."""
    if not NUMBA_AVAILABLE:
        return "numba is not installed"
    if not no_colcorr:
        return "no_colcorr=False (colour correction not implemented in Numba)"
    a = np.atleast_2d(a)
    if np.any(a[:, 11] != 0.0):
        return "To2 != 0 (second thermal component source unverified)"
    if np.any(a[:, 12] != 0.0):
        return "K != 0 (second thermal component source unverified)"
    if np.any((a[:, 0] < -0.5) | (a[:, 0] > 5.5)):
        return "FuncIndx outside 0-5"
    if phase_type not in ('kelsall', 'skysurf'):
        return f"phase_type {phase_type!r} not supported"
    return None


_VERIFIED = {}


def _verify(data, a, group, kw, nsample=96, rtol=1e-9):
    """One-time check of the Numba kernel against zkernel_fast on a subset."""
    from zodisurf.zkernel_fast import zkernel_fast
    key = (kw['phase_type'], kw['L2'], kw['new_iso_comp'], kw['iso_comp_only'],
           hash(np.ascontiguousarray(a).tobytes()))
    if key in _VERIFIED:
        return _VERIFIED[key]

    n = len(data)
    idx = np.linspace(0, n - 1, min(nsample, n)).astype(int)
    sub, gsub = data[idx], group[idx]
    try:
        ref = zkernel_fast(sub, a, group=gsub, n_jobs=1, **kw)
        got = zkernel_numba(sub, a, group=gsub, _verified=True, **kw)
        denom = np.maximum(np.abs(ref), 1e-30)
        ok = bool(np.all(np.abs(got - ref) / denom < rtol))
        worst = float(np.max(np.abs(got - ref) / denom))
    except Exception as e:                             # noqa: BLE001
        ok, worst = False, float('nan')
        print(f"[zkernel_numba] verification raised {type(e).__name__}: {e}")
    if not ok:
        print(f"[zkernel_numba] verification FAILED (worst rel diff {worst:.3g}); "
              "falling back to the NumPy kernel.")
    _VERIFIED[key] = ok
    return ok


def warmup(phase_type='skysurf'):
    """Trigger JIT compilation (and populate the on-disk cache) up front."""
    from zodisurf.get_zmod import get_zmod
    get_zmod(1.0 if phase_type == 'skysurf' else 1.25, phase_type,
             200.0, 10.0, 5.0, engine='numba')


# ======================================================================
#  Public kernel
# ======================================================================
def zkernel_numba(data, a, phase_type='kelsall', L2=False, no_colcorr=True,
                  solar_irr=None, new_iso_comp=False, iso_comp_only=False,
                  group=None, n_jobs=None, _verified=False):
    """
    Numba equivalent of zkernel_fast. Raises RuntimeError if a guard trips;
    callers should fall back to zkernel_fast in that case.
    """
    from zodisurf.zkernelpy_original import earthsun

    reason = guard_reason(a, phase_type, no_colcorr)
    if reason is not None:
        raise RuntimeError(f"Numba kernel unavailable: {reason}")

    npts = len(data)
    if npts == 0:
        return np.zeros(0)

    lam_all = np.asarray(data['wave_len'], dtype=np.float64)
    day = np.asarray(data['day1990'], dtype=np.float64)
    lon_deg = np.asarray(data['longitude'], dtype=np.float64)
    lat_deg = np.asarray(data['latitude'], dtype=np.float64)

    a = np.atleast_2d(np.asarray(a, dtype=np.float64))
    if group is None:
        if a.shape[0] == 1:
            group = np.zeros(npts, dtype=np.int64)
        elif a.shape[0] == npts:
            a, group = np.unique(a, axis=0, return_inverse=True)
        else:
            raise ValueError("zpar must have 1 row, npts rows, or be given with `group`.")
    group = np.asarray(group, dtype=np.int64).ravel()

    kw = dict(phase_type=phase_type, L2=L2, no_colcorr=no_colcorr,
              solar_irr=solar_irr, new_iso_comp=new_iso_comp,
              iso_comp_only=iso_comp_only)
    if not _verified and not _verify(data, a, group, kw):
        raise RuntimeError("Numba kernel failed self-verification")

    if n_jobs is not None and NUMBA_AVAILABLE:
        set_num_threads(max(1, int(n_jobs)))

    # ---- per-pixel geometry ----
    SolElong, Earth_Dis, Earth_Lon, Earth_Mean_Lon, L2_Dis = earthsun(day, lon_deg, lat_deg)
    Re = np.ascontiguousarray(np.broadcast_to(
        np.asarray(L2_Dis if L2 else Earth_Dis, dtype=np.float64), (npts,)))
    SolElong = np.broadcast_to(np.asarray(SolElong, dtype=np.float64), (npts,))
    Earth_Lon = np.broadcast_to(np.asarray(Earth_Lon, dtype=np.float64), (npts,))
    EML = np.ascontiguousarray(np.broadcast_to(
        np.asarray(Earth_Mean_Lon, dtype=np.float64), (npts,)))

    lat, lon = lat_deg * D2R, lon_deg * D2R
    COSLAT, SINLAT = np.cos(lat), np.sin(lat)
    COSLON, SINLON = np.cos(lon), np.sin(lon)
    X0 = Re * np.cos(Earth_Lon)
    Y0 = Re * np.sin(Earth_Lon)
    SIN_SE, COS_SE = np.sin(SolElong), np.cos(SolElong)

    B = 2.0 * (X0 * COSLAT * COSLON + Y0 * COSLAT * SINLON)
    C = Re**2 - RMAX**2
    with np.errstate(divide='ignore', invalid='ignore'):
        Q = -0.5 * B * (1.0 + np.sqrt(B**2 - 4.0 * C) / np.abs(B))
        RANGE = np.maximum(Q, C / Q)
    RANGE = np.where(B == 0.0, np.sqrt(-C), RANGE)
    nsteps = (np.rint(RANGE / _SIMP_H).astype(np.int64) + 1)

    # ---- group pixels so every group has scalar parameters ----
    cols = [group.astype(np.float64), lam_all]
    if solar_irr is not None:
        solar_irr = np.broadcast_to(
            np.asarray(solar_irr, dtype=np.float64).ravel(), (npts,))
        cols.append(solar_irr)
    keys, gid = np.unique(np.column_stack(cols), axis=0, return_inverse=True)
    gid = np.asarray(gid).ravel()

    f = np.zeros(npts)
    for g in range(len(keys)):
        pix = np.nonzero(gid == g)[0]
        if pix.size == 0:
            continue
        prm = a[int(keys[g, 0])]
        Lambda = keys[g, 1]
        sirr = keys[g, 2] if solar_irr is not None else None
        if phase_type == 'kelsall':
            hit = np.nonzero(DBWAVE == Lambda)[0]
            det = int(hit[0]) if hit.size else -1
            if det < 0:
                raise RuntimeError(
                    f"kelsall wavelength {Lambda} is not a DIRBE band; "
                    "phase-function indexing is ambiguous")
        else:
            det = -99

        P = _pack(prm, Lambda, det, phase_type, sirr, new_iso_comp, iso_comp_only)
        out = np.zeros(pix.size)
        _core(P,
              np.ascontiguousarray(nsteps[pix]),
              np.ascontiguousarray(X0[pix]), np.ascontiguousarray(Y0[pix]),
              np.ascontiguousarray(Re[pix]),
              np.ascontiguousarray(COSLAT[pix]), np.ascontiguousarray(SINLAT[pix]),
              np.ascontiguousarray(COSLON[pix]), np.ascontiguousarray(SINLON[pix]),
              np.ascontiguousarray(SIN_SE[pix]), np.ascontiguousarray(COS_SE[pix]),
              np.ascontiguousarray(EML[pix]), out)
        f[pix] = out
    return f