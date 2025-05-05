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
;     no_colcorr=no_colcorr - if set, the result will be actual
;                  intensity at the specified wavelength(s) instead of
;                  quoted intensity for a spectrum nu*F_nu = constant
;              
;---------------------------------------------------------------------
function get_zmod,lambda,phase_type,day,lon,lat,zpar=zpar,no_colcorr=no_colcorr,solar_irr=solar_irr

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

;if not isa(lambda_str,/string) then begin
;  lambda = lambda_str
;  ;print,lambda
;endif

;if phase_type eq 'kelsall' then begin
dbwave = [1.25,2.2,3.5,4.9,12,25,60,100,140,240]
;endif

if phase_type eq 'skysurf' then begin
  albedo = get_albedo(lambda)
  hong_params = get_hong_params(lambda)
  print,lambda,albedo,hong_params
  zpar = put_zpar(zpar,0,0,0,albedo,det1=0,hg3=hong_params)
endif

; HG3 parameters from mpfit_zodiscat.pro
;paramshg3 = [$
;  [  0.625472,  0.211267, -0.399582],$
;  [  0.643291,  0.215381, -0.336206],$
;  [  0.635869,  0.190329, -0.316795],$
;  [0.00494688, 0.0565555, 0.0181312],$
;  [0.00524653, 0.0570165, 0.0173618],$
;  [0.00504638, 0.0600067, 0.0145624]$
;  
;  ]

;if phase_type eq 'skysurf' then begin
;  
;; If "lambda_test" is a string, the continue
;if isa(lambda_str,/string) then begin
;  
;  ;print,lambda_str
;
;  ; Depending on the value of lambdra, do something specific
;  ; Make lambda is all uppercase
;  case (lambda_str) of
;
;    'BAND1'  : begin
;      lambda = 1.25
;      zpar = put_zpar(zpar,-0.9421,0.1214,-0.1648,0.20411940,det1=0,hg3=[paramshg3[*,0],paramshg3[*,3]])
;      
;    end
;    
;;    'acswfc/f814w'  : begin
;;      lambda = 0.8045
;;      hg3_params = 
;;      zpar = put_zpar(zpar,0,0,0,albedo,det1=0,hg3=hg3_params)
;;    end
;;    
;;    'acswfc/f850lp'  : begin
;;      lambda = 0.90332
;;      hg3_params = 
;;      zpar = put_zpar(zpar,0,0,0,albedo,det1=0,hg3=hg3_params)
;;    end
;    
;    'wfc3ir/f098m'  : begin
;      lambda = 0.98647
;      albedo = 0.1482474433986153
;      hg3_params = [0.009999999776482582,-0.4955938520979887,-0.8955820710274301,0.07143416511169320,0.007934861581750215,0.0002084439024563433]
;      zpar = put_zpar(zpar,0,0,0,albedo,det1=0,hg3=hg3_params)
;    end
;    
;    'wfc3ir/f105w'  : begin
;      lambda = 1.0551
;      albedo = 0.1698644022678380
;      hg3_params = [0.47,-0.23,-0.86,0.05,0.03,0.00035]
;      zpar = put_zpar(zpar,0,0,0,albedo,det1=0,hg3=hg3_params)
;    end
;    
;    ; Bad fit
;;    'wfc3ir/f110w'  : begin
;;      lambda = 1.15345
;;      albedo = 0.1953394674340974
;;      hg3_params = [0.6138175510672645,0.2054723503224905,-0.3634739407969614,-0.8641126615473673,0.004930593358436438,0.05724403628584696,0.01714930010504230,0.0002532280641650158]
;;      zpar = put_zpar(zpar,0,0,0,albedo,det1=0,hg3=hg3_params)
;;    end
;    
;    'wfc3ir/f125w'  : begin
;      lambda = 1.24861
;      albedo = 0.1891574685552602
;      hg3_params = [0.2487920277988392,-0.2913825123625323,-0.8671190992280190,0.05359972225980716,0.02566414150864571,0.0003136064545263614]
;      zpar = put_zpar(zpar,0,0,0,albedo,det1=0,hg3=hg3_params)
;    end
;      
;    'wfc3ir/f140w'  : begin
;      lambda = 1.39229
;      albedo = 0.1961424921714678
;      hg3_params = [0.2741750549855834,-0.2792601659354249,-0.8992027520615717,0.05158500586528676,0.02776766154508247,0.0002248007702315055]
;      zpar = put_zpar(zpar,0,0,0,albedo,det1=0,hg3=hg3_params)
;    end
;    
;    'wfc3ir/f160w'  : begin
;      lambda = 1.53692
;      albedo = 0.2183097777740185
;      hg3_params = [0.3131078001014243,-0.2504466316217165,-0.8586292987743477,0.04957965051501052,0.02963613933584503,0.0003616822039085169]
;      zpar = put_zpar(zpar,0,0,0,albedo,det1=0,hg3=hg3_params)
;    end
;
;  endcase
;  
;;  dbwave[0] = lambda
;  
;endif

data = mk_zdata(lambda,day,lon,lat)

dbwave[0] = lambda
zkernel,data,zpar,z,no_colcorr=no_colcorr,dbwave=dbwave,losinfo=0,solar_irr=solar_irr

; Only return up to 6 decimal places
rounded_z_str = STRING(z, FORMAT='(F0.7)')
rounded_z_float = FLOAT(rounded_z_str)

return,rounded_z_float
end
