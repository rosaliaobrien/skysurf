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
    
# def measureskyregion(data,axis=0):
#     #first, perform 1 iteration of 3-sigma clipping of the data (with respect to the median)
#     sigclipped1=sigma_clip(data,masked=False,axis=axis,cenfunc='mean',sigma_lower=5,sigma_upper=3,maxiters=1)
#     sigclipped2=sigma_clip(data,masked=False,axis=axis,cenfunc='median',sigma_lower=5,sigma_upper=3,maxiters=10,stdfunc=lambda dt,axis=0: .5*(np.nanpercentile(dt,84.1,axis=axis)-np.nanpercentile(dt,15.9,axis=axis)))

#     #return the median and standard deviation of the sigma-clipped data
#     #return np.nanmedian(sigclipped1,axis=axis), np.nanstd(sigclipped1,axis=axis)
#     print(np.nanmedian(sigclipped2,axis=axis),.5*(np.nanpercentile(sigclipped2,84.1)-np.nanpercentile(sigclipped2,15.9)))
#     return np.nanmedian(sigclipped2,axis=axis),.5*(np.nanpercentile(sigclipped2,84.1)-np.nanpercentile(sigclipped2,15.9)) ,np.nanpercentile(sigclipped1,30.8,axis=axis),np.nanpercentile(sigclipped1,69.1,axis=axis)

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

    wf=np.where(np.isfinite(newdat))
    #plt.hist(newdat[wf],bins=1000)
    #plt.show()
    #print(len(np.where(np.isfinite(newdat))[0]))
    #return modest.mode(newdat[wf]),sigi
    return np.nanmedian(newdat),sigi
    
def measureskyregion(data,axis=0):
    vals=[]
    sigs=[]
    
    for i in range(len(data)):
        vali,sigi=measureskybin(data[i])
        vals.append(vali)
        sigs.append(sigi)
    return np.array(vals),np.array(sigs),np.array(sigs),np.array(sigs)
