import numpy as np
from zodi_model.get_zmod import get_zmod, read_zpars, mk_zdata

def check_wavelength(wavelength=1.6, day=1.0, lon=0.0, lat=45.0):
    """
    Check ZODI model values for a specific wavelength.
    
    Args:
        wavelength (float): Wavelength in microns (default: 1.25)
        day (float): Day since 1990 (default: 1.0)
        lon (float): Longitude in degrees (default: 0.0)
        lat (float): Latitude in degrees (default: 0.0)
    """
    lambda_val = np.array([wavelength])
    day = np.array([day])
    lon = np.array([lon])
    lat = np.array([lat])
    
    print(f"\nChecking ZODI model for:")
    print(f"Wavelength: {wavelength} microns")
    print(f"Day: {day[0]}")
    print(f"Longitude: {lon[0]}°")
    print(f"Latitude: {lat[0]}°")
    
    # Test both phase types
    print("\n" + "="*50)
    print("KELSALL MODEL:")
    zodi_kelsall = get_zmod(lambda_val, "kelsall", day, lon, lat)
    print(f"Result: {zodi_kelsall}")

    print("\n" + "="*50)
    print("SKYSURF MODEL:")
    zodi_skysurf = get_zmod(lambda_val, "skysurf", day, lon, lat)
    print(f"Result: {zodi_skysurf}")

    print("\n" + "="*50)
    print("COMPARISON:")
    print(f"Kelsall:  {zodi_kelsall[0]:.8f}")
    print(f"SKYSURF:  {zodi_skysurf[0]:.8f}")
    if zodi_kelsall[0] != 0:
        diff_pct = abs(zodi_skysurf[0] - zodi_kelsall[0]) / zodi_kelsall[0] * 100
        print(f"Difference: {diff_pct:.2f}%")

if __name__ == "__main__":
    # Check default values (1.25 microns)
    check_wavelength()
    
    # You can also check other wavelengths
    # check_wavelength(wavelength=2.0)
    # check_wavelength(wavelength=3.5) 