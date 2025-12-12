import numpy as np
import json
from scipy.interpolate import interp1d

from zodi_model.zkernelpy import zkernel, get_parnames
from zodi_model.skysurf_params import get_albedo, get_hong_params, get_mult, get_emiss, put_zpar
from zodi_model.solar_sp import solar_sp

def read_zpars(filename=None):
    """
    Reads the ZODI model parameters from a JSON file and returns them in the exact order expected by the model.
    Uses a mapping from model parameter names to JSON keys. Any missing parameters are filled with 0.0 for albedo/emissivity, 1.0 for emissivity if needed, otherwise 0.0.
    """
    if filename is None:
        import os
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
        'B3_No': 'ignored_BZ0', 'B3_Dz': 'ignored_BZ1', 'B3_Dr': 'ignored_BZ2', 'B3_Ro': 'ignored_BZ3',
        'B3_Vi': 'ignored_BZ4', 'B3_Vr': 'ignored_BZ5', 'B3_Pi': 'ignored_BZ6', 'B3_Pr': 'ignored_BZ7', 'B3_P2r': 'ignored_BZ8',
        'B3_Omega': 'ignored_BZ9', 'B3_Incl': 'ignored_BZ10',
        'B3_Xo': 'ignored_BZ11', 'B3_Yo': 'ignored_BZ12', 'B3_Zo': 'ignored_BZ13',
        'B3_A1': 'ignored_BZ14', 'B3_A2': 'ignored_BZ15', 'B3_A3': 'ignored_BZ16', 'B3_A4': 'ignored_BZ17',
        'B3_E3': 'ignored_BZ18', 'B3_E4': 'ignored_BZ19', 'B3_E5': 'ignored_BZ20', 'B3_E6': 'ignored_BZ21',
        'B3_E7': 'ignored_BZ22', 'B3_E8': 'ignored_BZ23', 'B3_E9': 'ignored_BZ24', 'B3_E10': 'ignored_BZ25',
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
        'LB_Theta': 'leading_blob_theta', 'LB_dTheta': 'leading_blob_dtheta', 'LB_dZ': 'leading_blob_dZ',
        # Trailing Blob
        'TB_No': 'trailing_blob_n', 'TB_R': 'R_TB', 'TB_dR': 'sigma_rTB',
        'TB_Theta': 'theta_TB', 'TB_dTheta': 'sigma_thetaTB', 'TB_dZ': 'sigma_zTB',
        # RB Geometry
        'RB_Omega': 'OmegaRB', 'RB_Incl': 'iRB', 'RB_Xo': 'blob_X0', 'RB_Yo': 'blob_Y0', 'RB_Zo': 'blob_Z0',
        # RB Albedo and Emissivity
        'RB_A1': 'A1_SR', 'RB_A2': 'A2_SR', 'RB_A3': 'A3_SR', 'RB_A4': 'A4_SR',
        'RB_E3': 'E3_SR', 'RB_E4': 'E4_SR', 'RB_E5': 'E5_SR', 'RB_E6': 'E6_SR', 'RB_E7': 'E7_SR',
        'RB_E8': 'E8_SR', 'RB_E9': 'E9_SR', 'RB_E10': 'E10_SR'
    }
    try:
        with open(filename, 'r') as f:
            data = json.load(f)
        parnames = get_parnames()
        zpar = []
        for name in parnames:
            json_key = param_map.get(name, name)
            # Albedo and emissivity logic
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

def get_zmod(lambda_, phase_type, day, lon, lat, zpar=None, solar_irr=None, no_colcorr=False, new_iso_comp=False, iso_comp_only=False):
    """
    Compute the Kelsall et al. 1998 (ApJ,508,44) ZODI model intensity 
    for a given set of LOS, Wavelength, and Time.
    
    Parameters:
    -----------
    lambda_ : float
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
    
    # Check wavelength restriction for current implementation
    if phase_type == 'skysurf':
        if lambda_ > 3.5:
            raise ValueError("This code does not work at lambda > 3.5 micron.")
    
    # Input validation
    lambda_ = np.atleast_1d(lambda_)
    day = np.atleast_1d(day)
    lon = np.atleast_1d(lon)
    lat = np.atleast_1d(lat)
    
    if phase_type not in ['kelsall', 'skysurf']:
        raise ValueError("phase_type must be either 'kelsall' or 'skysurf'")

    if np.any(lon < 0) or np.any(lon > 360):
        raise ValueError("Longitude (lon) must be in the range 0 to 360 degrees.")
    if np.any(lat < -90) or np.any(lat > 90):
        raise ValueError("Latitude (lat) must be in the range -90 to 90 degrees.")

    # Load default parameters if not provided
    if zpar is None:
        zpar = read_zpars()
    else:
        zpar = np.array(zpar, dtype=float)

    # Array of DIRBE wavelengths (will be modified for skysurf)
    dbwave = np.array([1.25, 2.2, 3.5, 4.9, 12, 25, 60, 100, 140, 240])
    
    # Use the SKYSURF (O'Brien+2025) albedo and phase function
    if phase_type == 'skysurf':
        # Modify dbwave to accommodate custom wavelength
        dbwave[0] = lambda_[0]  # Use first wavelength value
        
        # Get SKYSURF parameters
        albedo = get_albedo(lambda_[0])
        hong_params = get_hong_params(lambda_[0])
        
        if lambda_[0] > 1.6:
            # Apply multiplicative factor for longer wavelengths
            mult = get_mult(lambda_[0])
            hong_params[-3:] = hong_params[-3:] * mult
            
            # Get emissivity
            emiss = get_emiss(lambda_[0])
            
            # Update parameters
            zpar = put_zpar(zpar, 0, 0, 0, albedo, det1=0, hg3=hong_params, E1=emiss)
        else:
            # Update parameters without emissivity
            zpar = put_zpar(zpar, 0, 0, 0, albedo, det1=0, hg3=hong_params)

    # Create data structure
    data = mk_zdata(lambda_, day, lon, lat)
    
    # Call kernel with phase_type parameter
    zodi = zkernel(data, zpar, phase_type=phase_type, no_colcorr=True,
                   dbwave=dbwave, solar_irr=solar_irr,
                   new_iso_comp=new_iso_comp, iso_comp_only=iso_comp_only)
    # NOTE: Correction factors were removed after discovering that Python
    # matches the raw IDL zkernel output exactly (0.11165 MJy/sr).
    # The discrepancy was with IDL's get_zmod wrapper (0.11494 MJy/sr),
    # which appears to have a calibration factor built in.
    # See WAVELENGTH_DEPENDENCY_ISSUE.md for full investigation.

    rounded_zodi = np.round(zodi, 5)

    return rounded_zodi

# Example usage
if __name__ == "__main__":
    lambda_val = 1.25
    phase_type = 'kelsall'  # or 'skysurf'
    day = np.array([1.0, 2.0, 3.0])
    lon = np.array([1.0, 10.0, 20.0])
    lat = np.array([1.0, 5.0, 10.0])

    zodi = get_zmod(lambda_val, phase_type, day, lon, lat)
    print("Zodiacal light intensity:", zodi)