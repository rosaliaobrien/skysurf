def hong_phase_func(x, a):
    """
    Hong phase function for SKYSURF model.
    Implementation based on O'Brien+2025 and IDL phase_func_fns.pro

    Parameters:
    -----------
    x : array-like
        Scattering angles
    a : array-like
        Hong parameters [g1, g2, g3, w1, w2, w3]

    Returns:
    --------
    phi : np.ndarray
        Phase function values
    """
    x = np.asarray(x, dtype=np.float64)
    a = np.asarray(a, dtype=np.float64)

    nf = len(a) // 2  # Number of components
    phi = np.zeros_like(x, dtype=np.float64)

    for i in range(nf):
        g = a[i]
        w = a[i + nf]
        # Henyey-Greenstein phase function component
        phi += w * (1 - g**2) / (1 + g**2 - 2*g*np.cos(x))**1.5

    return phi

def phasefunc(x, a):
    """
    Trial Phase Function for ZKERNEL (Kelsall model)
    Returns:
        phi (np.ndarray): phase function
        pder (np.ndarray): partial derivatives w.r.t. each parameter
    """
    x = np.asarray(x, dtype=np.float64)

    C0 = a[0]
    C1 = a[1]
    ka = a[2]

    pisq = PI**2
    picu = PI**3
    q0 = 2.0
    q1 = PI
    qk = (np.exp(ka * PI) + 1.0) / (ka**2 + 1.0)
    v = (2.0 * PI) * (C0 * q0 + C1 * q1 + qk)

    u = C0 + C1 * x + np.exp(ka * x)
    phi = u / v

    # Partial derivatives
    npts = len(x)

    pder = np.zeros((npts, 3), dtype=np.float64)

    # ∂phi/∂C0
    pder[:, 0] = (1.0 - phi * 4.0 * PI) / v

    # ∂phi/∂C1
    pder[:, 1] = (x - 2.0 * PI * phi * q1) / v

    # ∂phi/∂ka
    term1 = x * np.exp(ka * x)
    term2_numer = (ka**2 + 1.0) * PI * np.exp(ka * PI) - 2.0 * ka * (np.exp(ka * PI) + 1.0)
    term2_denom = (ka**2 + 1.0)**2
    term2 = 2.0 * PI * phi * term2_numer / term2_denom
    pder[:, 2] = (term1 - term2) / v

    return phi, pder

def scattfunc(phase_type, det, Lambda, LOS, R, Re, SolElong, a, solar_irr=None):
    """
    Compute scattered light contribution of zodiacal dust.

    Parameters:
    - phase_type: str
        'kelsall' or 'skysurf'
    - det: int
        DIRBE detector number (0–9) or -99 if non-DIRBE
    - Lambda: float
        Wavelength in microns
    - LOS: float or np.ndarray
        Line-of-sight distance from Earth (AU)
    - R: float or np.ndarray
        Distance from the Sun to the particle (AU)
    - Re: float
        Distance from Earth to Sun (AU)
    - SolElong: float
        Solar elongation (in radians)
    - a: np.ndarray
        Model parameters (assumed length at least 3 per detector)
    - df_out: np.ndarray or None
        Optional output array to receive partial derivatives [npts, nparams]
    - solar_irr: float or None
        Optional solar irradiance override

    Returns:
    - Scatt: np.ndarray
        Scattered light brightness per LOS element [npts]
    """
    if det <= 2:

        # Solar flux at 1 AU (W/m^2/um) - from IDL code
        SolFlux1AU = np.array([
            2.3405606e+08,  # 1.25 micron
            1.2309874e+08,  # 2.2 micron
            64292872.,      # 3.5 micron
            35733824.,      # 4.9 micron
            5763843.0,      # 12 micron
            1327989.4,      # 25 micron
            230553.73,      # 60 micron
            82999.336,      # 100 micron
            42346.605,      # 140 micron
            14409.608       # 240 micron
        ])

        # Solar flux calculation 
        if solar_irr is not None:
            SolFlux = solar_irr / R**2
        elif (phase_type == 'skysurf'):
            SolFlux = solar_sp(Lambda) / R**2
        elif det >= 0 and phase_type == 'kelsall':
            SolFlux = SolFlux1AU[det] / R**2
        else:
            SolFlux = 0.0

        # print(np.shape(R), np.shape(SolFlux), np.shape(Lambda))
        # print(R[:5])
        # print(SolFlux[:5])

        # Phase angle calculation (match IDL logic)
        phase_ang = np.arcsin(np.clip(Re / R * np.sin(SolElong), -1.0, 1.0))
        itest = (LOS >= Re * np.cos(SolElong))
        scat_ang = np.where(itest, np.pi - phase_ang, phase_ang)

        # Phase function parameter selection - match IDL logic
        if len(a) == 9:
            # Original Kelsall: extract 3 parameters based on detector
            aaa = a[3*det:3*det+3]
        elif len(a) % 2 == 0:
            # Hong functions: use all parameters
            aaa = a
        else:
            raise Exception('Your phase function parameters are messed up!')

        # Phase function - use Hong function for SKYSURF, regular phase function for Kelsall
        if phase_type == 'skysurf' and len(a) % 2 == 0:
            # SKYSURF uses Hong phase function
            phase_func = hong_phase_func(scat_ang, aaa)
        if phase_type == 'kelsall':
            # Use regular Kelsall phase function
            phase_func, _ = phasefunc(scat_ang, aaa)

        phase_func = np.asarray(phase_func)

        # Scattered light calculation using conditional solar flux
        Scatt = SolFlux * phase_func
        
    else:
        Scatt = np.zeros_like(LOS)
    
    return Scatt


def thermfunc(phase_type, det, Lambda, R, a, no_colcorr=False):

    eps = 1e-20

    # Temp, Temp2, idk, Temperature powerlaw
    # To2 and Kappa are zero in the Kelsall model :)
    To1, To2, Kappa, Delta = a[:4]

    # Conversion factor from W/cm^2/sr to MJy/sr
    cfact = Lambda / 3.0e-10

    # First component
    Temp1 = To1 / R**Delta
    Bnu1, dBdT1 = zPlanck(Temp1, Lambda, return_derivative=True)
    Bnu1 *= cfact
    dBdT1 *= cfact

    if (phase_type == 'kelsall') and (det >= 0) and not no_colcorr:
        CCtherm1 = colcorr(det, Temp1, return_derivative=True)
    else:
        CCtherm1 = 1.0, 0.0

    # Second component
    if Kappa != 0.0:
        Temp2 = To2 / R**Delta
        Bnu2, dBdT2 = zPlanck(Temp2, Lambda, return_derivative=True)
        Bnu2 *= cfact
        dBdT2 *= cfact

        if (phase_type == 'kelsall') and (det >= 0) and not no_colcorr:
            CCtherm2, dCdT2 = colcorr(det, Temp2, return_derivative=True)
        else:
            Temp2 = 0.0
            Bnu2 = eps
            dBdT2 = 0.0
            CCtherm2 = (1.0, 0.0)
            dCdT2 = 0.0
    else:
        Temp2 = 0.0
        Bnu2 = eps
        dBdT2 = 0.0
        CCtherm2 = (1.0, 0.0)
        dCdT2 = 0.0

    # Total thermal brightness
    # Note: Kappa is zero
    Therm = Bnu1 * CCtherm1[0] + Kappa * Bnu2 * CCtherm2[0]
    # dBdT = dBdT1 + Kappa * dBdT2

    # df = None
    # if want_partials:
    #     npts = len(R)
    #     npar = len(a)
    #     df = np.zeros((npts, npar), dtype=np.float64)

    #     dT1 = dBdT1 * CCtherm1[0] + Bnu1 * CCtherm1[1]
    #     dT2 = (dBdT2 * CCtherm2[0] + Bnu2 * CCtherm2[1]) * Kappa

    #     df[:, 0] = dT1 / R**Delta                      # To1
    #     df[:, 1] = dT2 / R**Delta                      # To2
    #     df[:, 2] = Bnu2 * CCtherm2[0]                     # Kappa
    #     df[:, 3] = -np.log(R) * (Temp1 * dT1 + Temp2 * dT2)  # Delta

    return Therm

def zsrcfunc(det, Scatt, Therm, a, phase_type='kelsall'):
    """
    Computes the zodiacal dust source function (brightness per particle).

    Parameters:
    - det: int
        DIRBE detector number (0–9) or -1 if non-DIRBE
    - Scatt: float or np.ndarray
        Scattered light contribution
    - Therm: float or np.ndarray
        Thermal emission contribution
    - a: np.ndarray
        Model parameters [Albedo 0–2, shared Albedo (3), Emissivities 4–11]
    - df_out: np.ndarray or None
        Optional output array for derivatives [npts, nparams]

    Returns:
    - Source: float or np.ndarray
        Brightness per particle (MJy/sr)
    - df_out (if provided): np.ndarray
        Partial derivatives w.r.t. parameters
    - dScatt (if provided): float or np.ndarray
        Derivative of source function w.r.t. scattered light
    - dTherm (if provided): float or np.ndarray
        Derivative w.r.t. thermal emission
    """

    # want_partials = df_out is not None
    a = np.asarray(a) # 3 albedos, rest are emissivities 
    npar = len(a)

    if phase_type == 'kelsall':
        # Albedo assignment (0–3 per detector, rest share a[3])
        AlbedoDet = np.concatenate([a[0:4], np.full(6, a[3])])

        # Emissivity assignment (detectors 0–1: emiss = 1.0, rest from a[4:11])
        EmissDet = np.concatenate([np.ones(2), a[4:12]])

        Albedo = AlbedoDet[det]
        Emiss = EmissDet[det]

    if phase_type == 'skysurf':

        Albedo = a[0]
        Emiss = np.concatenate([np.ones(2), a[4:12]])[2] # Use the albedo for COBE 3rd detector (3.5 micron) for now

    if phase_type == 'jwst_fitting':

        Albedo = a[0]
        Emiss = a[4]


    # Total brightness
    Source = Albedo * Scatt + Emiss * (1.0 - Albedo) * Therm

    return Source
