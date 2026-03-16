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
    Get albedo value for a given wavelength (or array of wavelengths) 
    in the SKYSURF model.
    """
    # Ensure input is a numpy array for vectorization
    wavelength = np.atleast_1d(wavelength)
    
    # 1. Calculate the linear relationship for ALL wavelengths first
    # (This handles the wavelength <= 1.6 case)
    albedo = 0.11298 * wavelength + 0.08231

    # 2. Identify where wavelength > 1.6
    mask = wavelength > 1.6
    
    if np.any(mask):
        # Calculate albedo at 1.6 microns for the interpolation start point
        albedo_1p6 = 0.11298 * 1.6 + 0.08231
        
        lam_list = np.array([1.6, 2.2, 3.5, 4.9])
        albedo_list = np.array([albedo_1p6, 0.255, 0.210, 0])

        # Apply interpolation only to the elements that need it
        albedo[mask] = np.interp(wavelength[mask], lam_list, albedo_list)

    # Return a scalar if the input was a scalar, otherwise return the array
    return albedo[0] if albedo.size == 1 and not isinstance(wavelength, np.ndarray) else albedo

def get_hong_params(wave_arr):
    """
    Get Hong phase function parameters for a given wavelength array.
    From O'Brien+2025 - exact implementation from IDL code.

    Parameters:
    -----------
    wave_arr : np.ndarray or float
        Wavelength in microns.

    Returns:
    --------
    np.ndarray
        Array of Hong phase function parameters with shape (N, 6) 
        containing [g1, g2, g3, w1, w2, w3].
    """
    # Ensure input is a numpy array for vectorization
    wave_arr = np.atleast_1d(wave_arr)
    
    # Clamp to valid range [0.25, 1.6] without modifying the original input
    lam_use = np.clip(wave_arr, 0.25, 1.6)

    # Calculate phase function parameters using linear relationships
    g1 = 0.24958 * lam_use + 0.11571
    g2 = 0.05428 * lam_use - 0.30864
    g3 = np.full_like(lam_use, -0.87036)  # Ensure constants are same shape as input
    w1 = 0.00183 * lam_use + 0.04775
    w2 = -0.00143 * lam_use + 0.03122
    w3 = np.full_like(lam_use, 0.00030)

    # Stack into a (N, 6) array
    # axis=-1 ensures that for a 1D input of length N, we get (N, 6)
    hg_arr = np.stack([g1, g2, g3, w1, w2, w3], axis=-1)

    return hg_arr

def get_mult(wavelength):
    """
    Get multiplicative factor for wavelengths. Vectorized version.
    From O'Brien+2025 - exact implementation from IDL code.
    
    Parameters:
    -----------
    wavelength : float or ndarray
        Wavelength in microns
    
    Returns:
    --------
    float or ndarray
        Multiplicative factor(s)
    """
    # Ensure input is a numpy array for consistency
    wavelength = np.atleast_1d(wavelength)
    
    # Define interpolation points
    lam_list = np.array([1.6, 2.2, 3.5, 4.9])
    mult_list = np.array([1.0, 1.227, 1.259, 1.0])
    
    # np.interp handles the array and linear interpolation.
    # left=1.0 ensures any wavelength < 1.6 returns 1.0.
    # right=1.0 ensures any wavelength > 4.9 returns 1.0.
    mult = np.interp(wavelength, lam_list, mult_list, left=1.0, right=1.0)
    
    return mult.reshape(-1, 1) if mult.size > 1 else [mult.item()]

def get_emiss(wavelength):
    """
    Get emissivity value for wavelengths > 1.6 microns.
    From O'Brien+2025 - exact implementation from IDL code.
    When extrapolating out to 3.5 micron, need to adjust the emissivity too.

    Parameters:
    -----------
    wavelength : arr
        Wavelength in microns

    Returns:
    --------
    float
        Emissivity value
    """
    wavelength = np.atleast_1d(wavelength)
    
    # Define wavelength and corresponding emissivity arrays
    lam_list = np.array([1.6, 2.2, 3.5, 4.9])
    emm_list = np.array([0, 0, 1.66, 0.997])

    emm = np.interp(wavelength, lam_list, emm_list)

    # This ensures anything <= 1.6 is 1.0, overriding the 0.0 interp result.
    emm = np.where(wavelength > 1.6, emm, np.nan)

    # Return as a scalar if input was a scalar, otherwise as an array
    return emm if emm.size > 1 else emm.item()

def put_zpar(zpar, PF1_C0, PF1_C1, PF1_Ka, albedo, det1=None, hg3=None, E1=None):
    """
    Update model parameters array with new values, supporting multiple albedos.
    
    Returns:
    --------
    np.ndarray
        Updated parameter matrix of shape (N, len(zpar))
    """
    # Ensure albedo is at least 1D to determine the number of objects (N)
    albedo = np.atleast_1d(albedo)
    n_objects = len(albedo)
    
    # Create a 2D array by tiling the template zpar: shape (N, len(zpar))
    aend = np.tile(zpar, (n_objects, 1))

    # Determine offset based on detector
    offset = 1 if det1 is not None and det1 != 0 else 0

    # update Phase Function (broadcasts the same 3 values to all N rows)
    aend[:, 1+3*offset : 4+3*offset] = [PF1_C0, PF1_C1, PF1_Ka]

    # We use np.newaxis to align the (N,) albedo array with the (N, 6) slice
    albedo_indices = np.array([33, 59, 85, 111, 137, 170]) + offset
    aend[:, albedo_indices] = albedo[:, np.newaxis]

    # 3. Update Emissivity
    if E1 is not None:
        E1 = np.atleast_1d(E1)
        emiss_indices = albedo_indices + 4

        # E1 is set to NaN during get_emiss for positions that do not need to be updated
        mask = ~np.isnan(E1)

        # Update emissivities only where E1 != NaN
        aend[np.ix_(mask, emiss_indices)] = E1[mask, np.newaxis]

    # 4. Add Hong Parameters (hg3)
    if hg3 is not None:
        hg3 = np.atleast_1d(hg3)
        # If hg3 is 1D [g1, g2...], same params for all N objects.
        # If hg3 is 2D (N, 6), unique params for each object.
        start_idx = 183
        end_idx = start_idx + hg3.shape[-1]
        
        if hg3.ndim == 1:
            aend[:, start_idx:end_idx] = hg3
        else:
            aend[:, start_idx:end_idx] = hg3

    return aend