"""
zkernel_ref.py
--------------
Per-pixel REFERENCE kernel: the original zkernel loop (scalar Re per pixel,
simpint with 0.025 AU steps, as in the IDL code). Slow; used only for
validating zkernel_fast. Physics functions are imported unchanged from zkernelpy.
"""
import numpy as np

from zodisurf.zkernelpy import (earthsun, simpint, scattfunc, thermfunc, zsrcfunc,
                                zcloud, migband, solring, new_isocloud)

RMAX = 5.2
DBWAVE = [1.25, 2.2, 3.5, 4.9, 12., 25., 60., 100., 140., 240.]


def zkernel_ref(data, a, phase_type='kelsall', L2=False, no_colcorr=True,
                solar_irr=None, new_iso_comp=False, iso_comp_only=False):
    d2r = np.pi / 180.0
    npts = len(data)

    a = np.atleast_2d(np.asarray(a, dtype=np.float64))
    if a.shape[0] == 1 and npts > 1:
        a = np.broadcast_to(a, (npts, a.shape[1]))
    if a.shape[0] != npts:
        raise ValueError(f"zpar has {a.shape[0]} rows, expected {npts}.")

    FuncIndx = a[:, 0]
    if phase_type == 'skysurf':
        aScatt = a[:, 183:189]
        if not np.any(aScatt != 0):
            raise ValueError("SKYSURF mode requires non-zero Hong parameters in zpar[183:189].")
    elif phase_type == 'kelsall':
        aScatt = a[:, 1:10]
    else:
        raise ValueError("phase_type must be 'kelsall' or 'skysurf'")

    aTherm = a[:, 10:14]
    aDens_C, aSrc_C = a[:, 14:33], a[:, 33:45]
    aDens_B1, aSrc_B1 = a[:, 45:59], a[:, 59:71]
    aDens_B2, aSrc_B2 = a[:, 71:85], a[:, 85:97]
    aDens_B3, aSrc_B3 = a[:, 97:111], a[:, 111:123]
    aDens_B4, aSrc_B4 = a[:, 123:137], a[:, 137:149]
    aDens_RB = a[:, np.r_[149:170, 182]]
    aSrc_RB = a[:, 170:182]

    SolElong, Earth_Dis, Earth_Lon, Earth_Mean_Lon, L2_Dis_from_sun = earthsun(
        data['day1990'], data['longitude'], data['latitude'])

    if solar_irr is not None:
        solar_irr = np.broadcast_to(np.asarray(solar_irr, dtype=np.float64).ravel(), (npts,))

    if phase_type == 'skysurf':
        Det = -99
    else:
        wl = np.asarray(data['wave_len'])
        detnum = np.zeros(npts, dtype=int) - 1
        for i in range(10):
            detnum += (wl == DBWAVE[i]) * (i + 1)

    f = np.zeros(npts)

    for ilos in range(npts):
        lat = data[ilos]['latitude'] * d2r
        lon = data[ilos]['longitude'] * d2r

        Re = L2_Dis_from_sun[ilos] if L2 else Earth_Dis[ilos]
        Theta = Earth_Lon[ilos]

        COSTHETA, SINTHETA = np.cos(Theta), np.sin(Theta)
        COSLON, SINLON = np.cos(lon), np.sin(lon)
        COSLAT, SINLAT = np.cos(lat), np.sin(lat)

        X0 = Re * COSTHETA
        Y0 = Re * SINTHETA
        B = 2.0 * (X0 * COSLAT * COSLON + Y0 * COSLAT * SINLON)
        C = Re**2 - RMAX**2
        Q = -0.5 * B * (1.0 + np.sqrt(B**2 - 4.0 * C) / abs(B))
        RANGE = max([Q, C / Q])

        los, gqwts = simpint(RANGE)

        Sxy = los * COSLAT
        X = Re * COSTHETA + Sxy * COSLON
        Y = Re * SINTHETA + Sxy * SINLON
        Z = los * SINLAT
        R = np.sqrt(X**2 + Y**2 + Z**2)

        Lambda = data[ilos]['wave_len']
        if phase_type == 'kelsall':
            Det = detnum[ilos]
        sirr = solar_irr[ilos] if solar_irr is not None else None

        Scatt = scattfunc(phase_type, Det, Lambda, los, R, Re, SolElong[ilos],
                          aScatt[ilos], solar_irr=sirr)
        Therm = thermfunc(phase_type, Det, Lambda, R, aTherm[ilos], no_colcorr=no_colcorr)

        Dens_C = zcloud(X, Y, Z, R, aDens_C[ilos], FuncIndx=FuncIndx[ilos])
        Src_C = zsrcfunc(Det, Scatt, Therm, aSrc_C[ilos], phase_type=phase_type)

        if new_iso_comp:
            Dens_new = new_isocloud(X, Y, Z, R, aDens_C[ilos])
        else:
            Dens_new = np.zeros_like(R)

        Dens_B1, _ = migband(X, Y, Z, R, aDens_B1[ilos])
        Src_B1 = zsrcfunc(Det, Scatt, Therm, aSrc_B1[ilos], phase_type=phase_type)
        Dens_B2, _ = migband(X, Y, Z, R, aDens_B2[ilos])
        Src_B2 = zsrcfunc(Det, Scatt, Therm, aSrc_B2[ilos], phase_type=phase_type)
        Dens_B3, _ = migband(X, Y, Z, R, aDens_B3[ilos])
        Src_B3 = zsrcfunc(Det, Scatt, Therm, aSrc_B3[ilos], phase_type=phase_type)
        Dens_B4, _ = migband(X, Y, Z, R, aDens_B4[ilos])
        Src_B4 = zsrcfunc(Det, Scatt, Therm, aSrc_B4[ilos], phase_type=phase_type)
        Dens_RB, _ = solring(X, Y, Z, R, Earth_Mean_Lon[ilos], aDens_RB[ilos])
        Src_RB = zsrcfunc(Det, Scatt, Therm, aSrc_RB[ilos], phase_type=phase_type)

        if iso_comp_only:
            Dens_C = Dens_B1 = Dens_B2 = Dens_B3 = Dens_B4 = Dens_RB = np.zeros_like(R)

        Flux = (Src_C * Dens_C + Src_B1 * Dens_B1 + Src_B2 * Dens_B2 + Src_B3 * Dens_B3
                + Src_B4 * Dens_B4 + Src_RB * Dens_RB + Src_C * Dens_new)

        f[ilos] = np.sum(gqwts * Flux)

    return f