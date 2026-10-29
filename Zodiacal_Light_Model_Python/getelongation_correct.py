import astropy.units as u
from astropy.coordinates import SkyCoord, get_body_barycentric, CartesianRepresentation
from astropy.time import Time
import numpy as np

def getelongation_correct(obstime,elon,elat,obsxyz,timeformat='isot'):
    target = SkyCoord(lat=elat*u.deg,lon=elon*u.deg, frame='barycentricmeanecliptic')

    # Convert barycentric coordinates from km to AU directly during initialization
    cart_au = CartesianRepresentation(
        (obsxyz[0] * u.km).to(u.AU),
        (obsxyz[1] * u.km).to(u.AU),
        (obsxyz[2] * u.km).to(u.AU)
    )

    # Observation times
    obs_time = Time(obstime, format=timeformat, scale='utc')
	
    # Get the Sun's barycentric position at the observation time (returns in AU)
    sun_bary = get_body_barycentric('sun', obs_time)

    #transform sun position from barycentric frame to telescope frame
    sun_jwst_cart = sun_bary - cart_au
    sun_jwst_coord = SkyCoord(sun_jwst_cart, frame='icrs')

    #elongation is angle separating sun from target
    solar_elongation = target.separation(sun_jwst_coord)

    return solar_elongation
