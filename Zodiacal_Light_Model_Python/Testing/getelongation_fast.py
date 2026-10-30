import numpy as np
KMINAU = 149597870.7

def getelongation_fast(elon,elat,obsxyz):
    d2r=np.pi/180

    #earth obliquity (for tranforming between ircs and ecliptic coordiantes)
    obliquity = 23.4392794 * np.pi/180
    
    #obsxyz in heliocentric equitorial, want heliocentric equitorial
    pos_ec = np.array([obsxyz[0]/KMINAU,obsxyz[1]/KMINAU*np.cos(obliquity) + obsxyz[2]/KMINAU*np.sin(obliquity),-obsxyz[1]/KMINAU*np.sin(obliquity) + obsxyz[2]/KMINAU*np.cos(obliquity)])

    obs_Dis = np.sqrt((obsxyz[0]/KMINAU)**2 + (obsxyz[1]/KMINAU)**2 + (obsxyz[2]/KMINAU)**2)

    #dot product to get angle between observatory-sun and observatory-target vectors
    dot_product = (
        -pos_ec[0] * np.cos(elat*d2r) * np.cos(elon*d2r) +
        -pos_ec[1] * np.cos(elat*d2r) * np.sin(elon*d2r) +
        -pos_ec[2] * np.sin(elat*d2r)
        )

    elongation = np.arccos(dot_product/obs_Dis)
    
    return elongation

    
