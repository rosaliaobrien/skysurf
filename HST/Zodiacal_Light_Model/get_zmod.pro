@zmodel_init
;---------------------------------------------------------------------
; Function GET_ZMOD
;
; Function to compute the Kelsall et al. 1998 (ApJ,508,44) ZODI model 
; intensity for a given set of LOS, Wavelength, and Time.
;
; Written By:   BA Franz
; Date:   23 January 1995
; 
; 
; Modified By: Rosalia O'Brien
; Date: 05 May 2025
;
; Usage:
;     zodi = GET_ZMOD(lambda,day,lon,lat)
;
; Inputs:
;     lambda     - [float] wavelength in microns
;     phase_type - [str] "kelsall" or "skysurf"
;                  Whether to use the standard Kelsall+1998 phase function 
;                  or the O'Brien+2025 phase function. To use the Kelsall 
;                  phase function, the wavelength must be one of DIRBE's 
;                  nominal wavelengths: 1.25, 2.2, 3.5, 4.9, 12, 25, 60, 
;                  100, 140, 240
;     day        - [arr] 1990 day number(s), where 1.0 = 1 Jan 1990
;     lon        - [arr] ecliptic longitude(s)
;     lat        - [arr] ecliptic latitude(s)
;
;     note: lambda, day, lon, lat can be any mixture of scalars or arrays, 
;           but all arrays must be of same length.
;
; Outputs:
;
;     zodi - scalar or array of zodi model intensity in MJy sr-1
;            For comparison with DIRBE data, this is 'quoted' intensity
;            at the nominal wavelength(s) for an assumed source
;            spectrum nu*F_nu = constant
;
; Optional Keyword Inputs:
;
;     zpar=zpar - array of zodi model params.  If not specified, the
;                 default will be restored from zpars.xdr
;              
;---------------------------------------------------------------------
function get_zmod,lambda,phase_type,day,lon,lat,zpar=zpar,solar_irr=solar_irr

if (n_elements(zpar) eq 0) then begin
    aend    = 0
    astart  = 0
    chisqr  = 0
    covar   = 0
    sigma_a = 0
    freevars= 0
    wiring  = 0
    restore,'zpars.xdr'
    zpar = aend
endif

; Array of DIRBE wavelengths
dbwave = [1.25,2.2,3.5,4.9,12,25,60,100,140,240]

; Use the SKYSURF (O'Brien+2025) albedo and phase function
if phase_type eq 'skysurf' then begin
  
  
  dbwave[0] = lambda ; Modify dbwave this way to avoid breaking the code elsewhere
  
  albedo = get_albedo(lambda)
  hong_params = get_hong_params(lambda)
  
  if lambda gt 1.6 then begin
    
    mult = get_mult(lambda)
    hong_params[-3:*] = hong_params[-3:*]*mult
    
    emm = get_emiss(lambda)
    
    zpar = put_zpar(zpar,0,0,0,albedo,det1=0,hg3=hong_params,E1=emm)
  
  endif else begin
  
    zpar = put_zpar(zpar,0,0,0,albedo,det1=0,hg3=hong_params)
  
  endelse
  
endif

data = mk_zdata(lambda,day,lon,lat)

;     no_colcorr=no_colcorr - if set, the result will be actual
;                  intensity at the specified wavelength(s) instead of
;                  quoted intensity for a spectrum nu*F_nu = constant
if phase_type eq 'skysurf' then begin
  zkernel,data,zpar,phase_type,z,no_colcorr=1,dbwave=dbwave,losinfo=0,solar_irr=solar_irr
endif else if phase_type eq 'kelsall' then begin
  zkernel,data,zpar,phase_type,z,no_colcorr=1,dbwave=dbwave,losinfo=0,solar_irr=solar_irr
endif

; Only return up to 6 decimal places
rounded_z_str = STRING(z, FORMAT='(F0.7)')
rounded_z_float = FLOAT(rounded_z_str)

return,rounded_z_float
end
