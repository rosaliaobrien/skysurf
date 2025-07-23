#Author:Tim Carleton, based on Zak Goisman's Javascript HST completeness calculator
#If you use this, please cite Goisman et al. 2025

import numpy as np
from scipy.special import erfc, erfcinv
p9arg=erfcinv(2*0.9)
def completeness(camera,filt,exp,size,background=-1):

    if type(size)==np.ndarray:
        wpoint=np.where(size<0.01)
        if len(wpoint[0])>0:
            print('Fitting function results unreliable for size<0.01", defaulting to point-source size of 0.01"')
            usesize=np.array(size)
            usesize[w]=0.01
        else:
            usesize=np.array(size)
    else:
        if size<.01:
            print('Fitting function results unreliable for size<0.01", defaulting to point-source size of 0.01"')
            usesize=0.01
        else:
            usesize=size
    dat=np.genfromtxt('completeness_data/'+camera.upper()+'table.csv',delimiter=',',dtype='str')
    date=np.genfromtxt('completeness_data/'+camera.upper()+'_elephant_data.csv',delimiter=',',dtype='str')
    rowi=np.where(dat[1:,0]=='ACS '+filt.upper())[0]
    if len(rowi)==0:
        print('Filter not included')
        return np.nan

    else:
        rowi+=1 # account for header row

    
    a_e=float(date[rowi,1])
    b_e=float(date[rowi,2])
    x0_e=float(date[rowi,3])
    c_e=float(date[rowi,4])

    expa=float(dat[rowi,1])
    backa=float(dat[rowi,2])
    avge=float(dat[rowi,22])
    avgb=float(dat[rowi,23])
        
    logsize=np.log10(usesize)
    
    magconstant=a_e+(b_e*(logsize-x0_e-1))/(1+np.exp(-c_e*(logsize-x0_e))) #magnitude limit for avgerage exposure at given size

    if background<0:
        finalmag=expa*(np.log10(exp)-avge)+magconstant #correct for exposure time and background
    else:
        finalmag=expa*(np.log10(exp)-avge)+backa*(np.log10(background)-avgb)+magconstant #correct for exposure time and background

    return finalmag

def logistic(a,b,x0,c,val):
    lgval=np.log10(val)
    return a+(b*(lgval-x0-1))/(1+np.exp(-c*(lgval-x0)))

def probdetect(camera,filt,exp,size,magnitude,background=-1):
    if type(size)==np.ndarray:
        wpoint=np.where(size<0.01)
        if len(wpoint[0])>0:
            print('Fitting function results unreliable for size<0.01", defaulting to point-source size of 0.01"')
            usesize=np.array(size)
            usesize[w]=0.01
        else:
            usesize=np.array(size)
    else:
        if size<.01:
            print('Fitting function results unreliable for size<0.01", defaulting to point-source size of 0.01"')
            usesize=0.01
        else:
            usesize=size

    dat=np.genfromtxt('completeness_data/'+camera.upper()+'table.csv',delimiter=',',dtype='str')
    date=np.genfromtxt('completeness_data/'+camera.upper()+'_elephant_data.csv',delimiter=',',dtype='str')
    rowi=np.where(dat[1:,0]=='ACS '+filt.upper())[0]
    if len(rowi)==0:
        print('Filter not included')
        return np.nan

    else:
        rowi+=1 # account for header row

    
    a_50=float(date[rowi,1])
    b_50=float(date[rowi,2])
    x0_50=float(date[rowi,3])
    c_50=float(date[rowi,4])

    a_90=float(date[rowi,5])
    b_90=float(date[rowi,6])
    x0_90=float(date[rowi,7])
    c_90=float(date[rowi,8])

    a_95=float(date[rowi,9])
    b_95=float(date[rowi,10])
    x0_95=float(date[rowi,11])
    c_95=float(date[rowi,12])

    mag50_e=logistic(a_50,b_50,x0_50,c_50,usesize)
    mag90_e=logistic(a_90,b_90,x0_90,c_90,usesize)
    mag95_e=logistic(a_95,b_95,x0_95,c_95,usesize)
    
    expa=float(dat[rowi,1])
    backa=float(dat[rowi,2])
    avge=float(dat[rowi,22])
    avgb=float(dat[rowi,23])
        
    logsize=np.log10(usesize)

    if background<0:
        mag50=expa*(np.log10(exp)-avge)+mag50_e #correct for exposure time and background
    else:
        mag50=expa*(np.log10(exp)-avge)+backa*(np.log10(background)-avgb)+mag50_e #correct for exposure time and background

    delta=mag50 - mag50_e

    mag50_s=mag50_e+delta
    mag90_s=mag90_e+delta
    mag95_s=mag95_e+delta

    m0=mag50_s
    sigma=(mag90_s-m0)/(np.sqrt(2)*p9arg)

    comp=0.5*erfc((magnitude-m0)/(np.sqrt(2)*sigma))*100
    

    return comp
