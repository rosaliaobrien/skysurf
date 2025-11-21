"""
SKYSURF Parameter Functions

Helper functions for the O'Brien+2025 SKYSURF phase function implementation.
These functions provide wavelength-dependent albedo, phase function parameters,
multiplicative factors, and emissivity values.

Author: Rosalia O'Brien
Python conversion: Tejovrash Acharya
"""

import numpy as np

def get_albedo(wavelength):
    """
    Get albedo value for a given wavelength in the SKYSURF model.
    From O'Brien+2025 - exact implementation from IDL code.

    Parameters:
    -----------
    wavelength : float
        Wavelength in microns

    Returns:
    --------
    float
        Albedo value for the wavelength
    """
    # Calculate albedo at 1.6 microns using linear formula
    albedo_1p6 = 0.11298 * 1.6 + 0.08231

    if wavelength > 1.6:
        # Define wavelength and corresponding albedo arrays for interpolation
        lam_list = np.array([1.6, 2.2, 3.5, 4.9])
        albedo_list = np.array([albedo_1p6, 0.255, 0.210, 0])

        # Interpolate between known values
        albedo = np.interp(wavelength, lam_list, albedo_list)
    else:
        # Linear relationship for wavelengths <= 1.6 microns
        albedo = 0.11298 * wavelength + 0.08231

    return albedo

def get_hong_params(wavelength):
    """
    Get Hong phase function parameters for a given wavelength.
    From O'Brien+2025 - exact implementation from IDL code.

    Parameters:
    -----------
    wavelength : float
        Wavelength in microns

    Returns:
    --------
    np.ndarray
        Array of Hong phase function parameters [g1, g2, g3, w1, w2, w3]
    """
    # If lambda is less than 0.25 or greater than 1.6, clamp to valid range
    # These phase functions have not been tested outside of this wavelength range
    if wavelength > 1.6:
        lam_use = 1.6
    elif wavelength < 0.25:
        lam_use = 0.25
    else:
        lam_use = wavelength

    # Calculate phase function parameters using linear relationships
    g1 = 0.24958 * lam_use + 0.11571
    g2 = 0.05428 * lam_use - 0.30864
    g3 = -0.87036
    w1 = 0.00183 * lam_use + 0.04775
    w2 = -0.00143 * lam_use + 0.03122
    w3 = 0.00030

    hg_arr = np.array([g1, g2, g3, w1, w2, w3])

    return hg_arr

def get_mult(wavelength):
    """
    Get multiplicative factor for wavelengths > 1.6 microns.
    From O'Brien+2025 - exact implementation from IDL code.
    Should be = 1 for nominal HST wavelengths. When extrapolating out to 3.5 micron,
    a multiplier is added to match Kelsall better.

    Parameters:
    -----------
    wavelength : float
        Wavelength in microns

    Returns:
    --------
    float
        Multiplicative factor
    """
    if wavelength > 1.6:
        # Define wavelength and corresponding multiplier arrays
        lam_list = np.array([1.6, 2.2, 3.5, 4.9])
        mult_list = np.array([1.0, 1.227, 1.259, 1.0])

        # Interpolate between known values
        mult = np.interp(wavelength, lam_list, mult_list)
    else:
        mult = 1.0

    return mult

def get_emiss(wavelength):
    """
    Get emissivity value for wavelengths > 1.6 microns.
    From O'Brien+2025 - exact implementation from IDL code.
    When extrapolating out to 3.5 micron, need to adjust the emissivity too.

    Parameters:
    -----------
    wavelength : float
        Wavelength in microns

    Returns:
    --------
    float
        Emissivity value
    """
    if wavelength > 1.6:
        # Define wavelength and corresponding emissivity arrays
        lam_list = np.array([1.6, 2.2, 3.5, 4.9])
        emm_list = np.array([0, 0, 1.66, 0.997])

        # Interpolate between known values
        emm = np.interp(wavelength, lam_list, emm_list)
    else:
        emm = 1.0

    return emm

def put_zpar(zpar, PF1_C0, PF1_C1, PF1_Ka, albedo, det1=None, hg3=None, E1=None):
    """
    Update model parameters array with new values.
    From IDL put_zpar.pro - exact implementation.
    Creates new aend=zpars array with new phase function and albedo for channel 1.

    Parameters:
    -----------
    zpar : np.ndarray
        Original parameter array (aend_in)
    PF1_C0, PF1_C1, PF1_Ka : float
        Phase function parameters (ignored when using hg3)
    albedo : float
        Albedo value (A1) used for all geometric components
    det1 : int, optional
        Detector number to replace (default: det=0, otherwise det=1)
    hg3 : np.ndarray, optional
        [g1,g2,g3,w1,w2,w3] for weighted sum of HG phase functions
    E1 : float, optional
        Emissivity value

    Returns:
    --------
    np.ndarray
        Updated parameter array
    """
    # Determine offset based on detector
    offset = 1 if det1 is not None and det1 != 0 else 0

    aend = zpar.copy()

    # Update phase function parameters (indices 1-3 for det=0, 4-6 for det=1)
    aend[1+3*offset:4+3*offset] = [PF1_C0, PF1_C1, PF1_Ka]

    # Update the albedos for SKYSURF at multiple component indices
    albedo_indices = np.array([33, 59, 85, 111, 137, 170]) + offset
    aend[albedo_indices] = albedo

    # Update Emissivity when using extrapolating "skysurf"
    if E1 is not None:
        emiss_indices = np.array([33+4, 59+4, 85+4, 111+4, 137+4, 170+4]) + offset
        aend[emiss_indices] = E1

    # Add phase function params for skysurf to the end of "a"
    if hg3 is not None and len(hg3) >= 4:
        # Store Hong parameters starting at index 183
        end_idx = 183 + len(hg3)
        if end_idx <= len(aend):
            aend[183:end_idx] = hg3

    return aend