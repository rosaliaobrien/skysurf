from astropy.io import fits
from astropy.stats import sigma_clipped_stats, sigma_clip, biweight
#import modest
from scipy.stats import median_abs_deviation

import numpy as np
import copy
import matplotlib.pyplot as plt

#Timothy Carleton
#2020

'''
This code is measures the typcial pixel values within a subregion
of a fits image

The default method for that is to return the biweight median
of the sigma-clipped values

To use:
>import measureskyregion
>measureskyregion.measureskyregion(data,axis=0)

data is a numpy array
axis=0 refers to the data axis that the statistics will be calculated

'''

def measureskybin(data,axis=0):
    newdat=copy.copy(data)
    w=np.where((data<-126.5) & (data>895.5))
    newdat[w]=np.nan
    loci=np.nanmean(newdat)
    sigi=np.nanstd(newdat)
    nout=1
    while nout>0:
        wi=np.where((newdat<loci-5*sigi) | (newdat>loci+3*sigi))
        nout=len(wi[0])
        newdat[wi]=np.nan
        loci=np.nanmedian(newdat)
        sigi=.5*(np.nanpercentile(newdat,84.1)-np.nanpercentile(newdat,15.9))

    return np.nanmedian(newdat),sigi
    
def measureskyregion(data,axis=0):
    vals=[]
    sigs=[]
    
    for i in range(len(data)):
        vali,sigi=measureskybin(data[i])
        vals.append(vali)
        sigs.append(sigi)
    return np.array(vals),np.array(sigs),np.array(sigs),np.array(sigs)
