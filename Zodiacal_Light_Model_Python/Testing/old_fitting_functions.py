def zkernel(data, a, phase_type='kelsall', L2 = False, indxpar=None, losinfo=False, no_colcorr=True, solar_irr=None, 
            new_iso_comp=False, iso_comp_only=False, 
            fitting_hong = None, fitting_albedo = None, fitting_em_band = None, fitting_em_smooth = None, fitting_em_ring = None):


    want_los_info = losinfo

    RMAX = 5.2
    NUMBER_OF_STEPS = 50

    npts = len(data)

    d2r = np.pi / 180.0
    eps = 1e-20
    if phase_type == 'kelsall':
        dbwave = [1.25, 2.2, 3.5, 4.9, 12., 25., 60., 100., 140., 240.]
    nmsg = 50000


    # DELETED by Rosalia
    # windx = zpar_wiring()
    # Apply wiring to zpars (critical for Kelsall mode)
    # This reindexes the zpars array to wire parameters together
    # For example: C_E3>C_A3 makes emissivity equal to albedo at 3.5 microns
    # Note: When no wiring is specified, this is just identity mapping
    # a = a[windx]


    FuncIndx = a[:,0]


    # Check if Hong parameters are stored at index 183+ (SKYSURF mode)
    # Only use Hong params if we're in skysurf mode AND they're non-zero
    if (phase_type == 'skysurf') | (phase_type == 'jwst_fitting'):
        hong_params = a[:,183:189]  # Get exactly 6 Hong parameters

        if np.any(hong_params != 0):  # Only use if non-zero
            aScatt = hong_params
    if phase_type == 'kelsall':
        nScatt = 9
        iScatt = np.arange(1, nScatt + 1)
        aScatt = a[:, iScatt]

    nTherm = 4
    iTherm = np.arange(10, 10 + nTherm)
    aTherm = a[:, iTherm]

    nDens_C = 19
    iDens_C = np.arange(14, 14 + nDens_C)
    aDens_C = a[:, iDens_C]
    nSrc_C = 12
    iSrc_C = np.arange(33, 33 + nSrc_C)
    aSrc_C = a[:, iSrc_C]


    nDens_B1 = 14
    iDens_B1 = np.arange(45, 45 + nDens_B1)
    aDens_B1 = a[:, iDens_B1]
    nSrc_B1 = 12
    iSrc_B1 = np.arange(59, 59 + nSrc_B1)
    aSrc_B1 = a[:, iSrc_B1]

    nDens_B2 = 14
    iDens_B2 = np.arange(71, 71 + nDens_B2)
    aDens_B2 = a[:, iDens_B2]
    nSrc_B2 = 12
    iSrc_B2 = np.arange(85, 85 + nSrc_B2)
    aSrc_B2 = a[:, iSrc_B2]

    nDens_B3 = 14
    iDens_B3 = np.arange(97, 97 + nDens_B3)
    aDens_B3 = a[:, iDens_B3]
    nSrc_B3 = 12
    iSrc_B3 = np.arange(111, 111 + nSrc_B3)
    aSrc_B3 = a[:, iSrc_B3]

    nDens_B4 = 14
    iDens_B4 = np.arange(123, 123 + nDens_B4)
    aDens_B4 = a[:, iDens_B4]
    nSrc_B4 = 12
    iSrc_B4 = np.arange(137, 137 + nSrc_B4)
    aSrc_B4 = a[:, iSrc_B4]

    nDens_RB = 21
    iDens_RB = np.arange(149, 149 + nDens_RB)
    iDens_RB = np.append(iDens_RB, 182)
    aDens_RB = a[:, iDens_RB]
    nSrc_RB = 12
    iSrc_RB = np.arange(170, 170 + nSrc_RB)
    aSrc_RB = a[:, iSrc_RB]

    SolElong, Earth_Dis, Earth_Lon, Earth_Mean_Lon, L2_Dis_from_sun = earthsun(data['day1990'], data['longitude'], data['latitude'])

    f = np.zeros(npts)

    if indxpar is not None:
        want_partials = True
        npar = len(a)
        if indxpar is None:
            nterms = npar
            indxpar = np.arange(nterms)
        else:
            nterms = len(indxpar)

        df = np.zeros((npts, nterms))
        parnum = np.arange(npar) - 1
        parnum[indxpar] = np.arange(nterms)
    else:
        want_partials = False

    if want_los_info:
        losdata = []

    # Detector assignment logic - only use for COBE data!
    # Still set for HST and JWST data to help with debugging
    if phase_type == 'skysurf':
        Det = -99
    if phase_type == 'jwst_fitting':
        Det = -99
    if phase_type == 'kelsall':
        detnum = np.zeros(npts, dtype=int) - 1
        for i in range(10):
            detnum += (np.array([d['wave_len'] for d in data]) == dbwave[i]) * (i + 1)

    # Loop through every position
    for ilos in range(npts):

        lat = data[ilos]['latitude'] * d2r
        lon = data[ilos]['longitude'] * d2r

        if L2 == True:
            Re = L2_Dis_from_sun[ilos]
        else:
            Re = Earth_Dis[ilos]
        Theta = Earth_Lon[ilos]

        COSTHETA = np.cos(Theta)
        SINTHETA = np.sin(Theta)
        COSLON = np.cos(lon)
        SINLON = np.sin(lon)
        COSLAT = np.cos(lat)
        SINLAT = np.sin(lat)

        X0 = Re * COSTHETA
        Y0 = Re * SINTHETA
        B = 2.0 * (X0 * COSLAT * COSLON + Y0 * COSLAT * SINLON)
        C = Re**2 - RMAX**2
        Q = -0.5 * B * (1.0 + np.sqrt(B**2 - 4.0 * C) / abs(B))
        RANGE = max([Q, C / Q])

        # CRITICAL FIX: IDL code uses simpint for ALL cases (gaussint is commented out)
        # The original Python code had conditional logic, but IDL always uses Simpson integration
        los, gqwts = simpint(RANGE)

        Sxy = los * COSLAT
        X = Re * COSTHETA + Sxy * COSLON
        Y = Re * SINTHETA + Sxy * SINLON
        Z = los * SINLAT
        R = np.sqrt(X**2 + Y**2 + Z**2)

        Lambda = data[ilos]['wave_len']

        if phase_type == 'kelsall':
            Det = detnum[ilos]

        # You have a unique solar irradiance value for each wavelength! 
        # So if its not None, you need to grab the solar_irr corresopnding to ilos
        # Same for the fitting functions
        if solar_irr is not None:
            solar_irr_ilos = solar_irr[ilos]
        else:
            solar_irr_ilos = None

        if fitting_hong is not None:
            fitting_hong_ilos = fitting_hong[ilos]
        else:
            fitting_hong_ilos = None

        if fitting_albedo is not None:
            fitting_albedo_ilos = fitting_albedo[ilos]
        else:
            fitting_albedo_ilos = None

        if fitting_em_smooth is not None:
            fitting_em_smooth_ilos = fitting_em_smooth[ilos]
        else:
            fitting_em_smooth_ilos = None

        if fitting_em_band is not None:
            fitting_em_band_ilos = fitting_em_band[ilos]
        else:
            fitting_em_band_ilos = None

        if fitting_em_ring is not None:
            fitting_em_ring_ilos = fitting_em_ring[ilos]
        else:
            fitting_em_ring_ilos = None

        # print(fitting_hong_ilos)
        # print(fitting_albedo_ilos)
        # print(fitting_em_ilos)

        Scatt = scattfunc(phase_type, Det, Lambda, los, R, Re, SolElong[ilos], aScatt[ilos], solar_irr=solar_irr_ilos, fitting_hong = fitting_hong_ilos)

        Therm = thermfunc(phase_type, Det, Lambda, R, aTherm[ilos], no_colcorr=no_colcorr)

        Dens_C = zcloud(X, Y, Z, R, aDens_C[ilos], FuncIndx=FuncIndx[ilos])

        Src_C = zsrcfunc(Det, Scatt, Therm, aSrc_C[ilos], phase_type=phase_type, fitting_albedo = fitting_albedo_ilos, fitting_em = fitting_em_smooth_ilos)

        # Isotropic component
        if new_iso_comp:
            Dens_new = new_isocloud(X, Y, Z, R, aDens_C[ilos])
        else:
            Dens_new = np.zeros_like(R)

        Dens_B1, dDens_B1 = migband(X, Y, Z, R, aDens_B1[ilos])
        Src_B1 = zsrcfunc(Det, Scatt, Therm, aSrc_B1[ilos], phase_type=phase_type, fitting_albedo = fitting_albedo_ilos, fitting_em = fitting_em_band_ilos)

        Dens_B2, dDens_B2 = migband(X, Y, Z, R, aDens_B2[ilos])
        Src_B2 = zsrcfunc(Det, Scatt, Therm, aSrc_B2[ilos], phase_type=phase_type, fitting_albedo = fitting_albedo_ilos, fitting_em = fitting_em_band_ilos)

        Dens_B3, dDens_B3 = migband(X, Y, Z, R, aDens_B3[ilos])
        Src_B3 = zsrcfunc(Det, Scatt, Therm, aSrc_B3[ilos], phase_type=phase_type, fitting_albedo = fitting_albedo_ilos, fitting_em = fitting_em_band_ilos)

        Dens_B4, dDens_B4 = migband(X, Y, Z, R, aDens_B4[ilos])
        Src_B4 = zsrcfunc(Det, Scatt, Therm, aSrc_B4[ilos], phase_type=phase_type, fitting_albedo = fitting_albedo_ilos, fitting_em = fitting_em_band_ilos)

        Dens_RB, dDens_RB = solring(X, Y, Z, R, Earth_Mean_Lon[ilos], aDens_RB[ilos])
        Src_RB = zsrcfunc(Det, Scatt, Therm, aSrc_RB[ilos], phase_type=phase_type, fitting_albedo = fitting_albedo_ilos, fitting_em = fitting_em_ring_ilos)

        # Isotropic component only mode
        if iso_comp_only:
            Dens_C = np.zeros_like(R)
            Dens_B1 = np.zeros_like(R)
            Dens_B2 = np.zeros_like(R)
            Dens_B3 = np.zeros_like(R)
            Dens_B4 = np.zeros_like(R)
            Dens_RB = np.zeros_like(R)

        Flux = (Src_C * Dens_C + Src_B1 * Dens_B1 + Src_B2 * Dens_B2 + Src_B3 * Dens_B3 + Src_B4 * Dens_B4 + Src_RB * Dens_RB + Src_C * Dens_new)

        # Sum the flux along the line of sight
        f[ilos] = np.sum(gqwts * Flux)

        # print(np.sum(Src_B1 * Dens_B1 + Src_B2 * Dens_B2 + Src_B3 * Dens_B3 + Src_B4 * Dens_B4 + Src_RB * Dens_RB + Src_C * Dens_new))

        if want_los_info:
            Dens = Dens_C + Dens_B1 + Dens_B2 + Dens_B3 + Dens_RB + Dens_new
            losdata.append({
                'wave_len': Lambda,
                'pixel_no': data[ilos]['pixel_no'],
                'lon': lon / d2r,
                'lat': lat / d2r,
                'elong': SolElong[ilos] / d2r,
                'day': data[ilos]['day1990'],
                'earth_lon': Theta / d2r,
                'source_rad': Re,
                'element': {
                    'earth_dist': los,
                    'solar_dist': R,
                    'ecliptic_z': Z,
                    'density': Dens,
                    'flux': Flux,
                    'int_wts': gqwts
                },
                'zodi': f[ilos]
            })

    if want_los_info:
        return losdata

    return f

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
            raise ValueError(f"{name} must be a scalar or match the length of the longest input ({npts}).")

    # Broadcast inputs
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

def get_zmod(lambda_, phase_type, day, lon, lat, L2 = False, zpar=None, solar_irr=None, new_iso_comp=False, iso_comp_only=False, 
             fitting_hong = None, fitting_albedo = None, fitting_em_band = None, fitting_em_smooth = None, fitting_em_ring = None):
    """
    Compute the Kelsall et al. 1998 (ApJ,508,44) ZODI model intensity 
    for a given set of LOS, Wavelength, and Time.
    
    Parameters:
    -----------
    lambda_ : float or array
        Wavelength in microns
    phase_type : str
        "kelsall" or "skysurf" - Whether to use the standard Kelsall+1998 
        phase function or the O'Brien+2025 phase function
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
    solar_irr : float, optional
        Solar irradiance (in MJy/sr) corresponding to bandpass
    no_colcorr : bool, optional
        If True, returns actual intensity instead of quoted intensity
    new_iso_comp : bool, optional
        If True, include isotropic component (default: False)
    iso_comp_only : bool, optional
        If True, model ONLY isotropic component (default: False)

    Returns:
    --------
    array_like
        Scalar or array of zodi model intensity in MJy sr^-1
    """
    
    # Input validation
    lambda_ = np.atleast_1d(lambda_)
    day = np.atleast_1d(day)
    lon = np.atleast_1d(lon)
    lat = np.atleast_1d(lat)

    # If lambda_ is a single value, broadcast it to the shape of the coordinates
    if lambda_.size == 1 and lambda_.shape != day.shape:
        lambda_ = np.full_like(day, lambda_[0], dtype=float)

    if solar_irr is not None:
        solar_irr = np.atleast_1d(solar_irr)

    # Check if lambda_ matches the shape of day, lon, and lat
    if not (lambda_.shape == day.shape == lon.shape == lat.shape):
        raise ValueError(
            f"Shape mismatch: lambda_ {lambda_.shape} must match "
            f"day {day.shape}, lon {lon.shape}, and lat {lat.shape}."
        )

    # Check wavelength restriction for current implementation
    if phase_type == 'skysurf':
        if any(element > 10 for element in lambda_):
            raise ValueError("This code does not work at lambda > 3.5 micron.")
    
    # if phase_type not in ['kelsall', 'skysurf']:
    #     raise ValueError("phase_type must be either 'kelsall' or 'skysurf'")

    if np.any(lon < 0) or np.any(lon > 360):
        raise ValueError("Longitude (lon) must be in the range 0 to 360 degrees.")
    if np.any(lat < -90) or np.any(lat > 90):
        raise ValueError("Latitude (lat) must be in the range -90 to 90 degrees.")

    # Load default parameters if noprovided
    if zpar is None:
        zpar = read_zpars()
    else:
        zpar = np.array(zpar, dtype=float)

    if phase_type == 'kelsall':
        # Array of DIRBE wavelengths
        dbwave = np.array([1.25, 2.2, 3.5, 4.9, 12, 25, 60, 100, 140, 240])
        zpar = np.array([zpar]*len(lambda_), dtype=float)
        # zpar = np.tile(zpar, (len(lambda_), 1))
        # print(zpar)
    
    # Use the SKYSURF (O'Brien+2025) albedo and phase function
    if (phase_type == 'skysurf') | (phase_type == 'jwst_fitting'):
        
        # Get SKYSURF parameters
        albedo_arr = get_albedo(lambda_)
        hong_params_arr = get_hong_params(lambda_)
        
        # if lambda_[0] > 1.6:
        # Apply multiplicative factor for longer wavelengths
        mult_arr = get_mult(lambda_)
        
        hong_params_arr[:,-3:] = hong_params_arr[:,3:] * mult_arr

        # Get emissivity
        emiss_arr = get_emiss(lambda_)
        
        # Update parameters
        zpar = put_zpar(zpar, 0, 0, 0, albedo_arr, det1=0, hg3=hong_params_arr, E1=emiss_arr)


    # Create data structure
    data = mk_zdata(lambda_, day, lon, lat)
    
    # Call kernel with phase_type parameter
    zodi = zkernel(data, zpar, phase_type=phase_type, L2 = L2, no_colcorr=True, solar_irr=solar_irr,
                   new_iso_comp=new_iso_comp, iso_comp_only=iso_comp_only, 
                   fitting_hong = fitting_hong, fitting_albedo = fitting_albedo, 
                   fitting_em_band = fitting_em_band, fitting_em_smooth = fitting_em_smooth, fitting_em_ring = fitting_em_ring)

    rounded_zodi = np.round(zodi, 5)

    return rounded_zodi