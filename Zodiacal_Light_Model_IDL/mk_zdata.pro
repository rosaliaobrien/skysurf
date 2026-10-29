function ZODI_str, Nstruct
 
  common RZODI_str, defined
 
  if N_elements( defined ) NE 1 then begin
 
       ZODI = { RZODI               , $  ; Structure Name
                PIXEL_NO:    0L     , $  ; DIRBE Pixel Number
                DAY1990:     0.0    , $  ; 1990 Day Number
                WAVE_LEN:    0.0    , $  ; Wave Length (microns)
                LONGITUDE:   0.0    , $  ; Geo Eclip Lon of LOS (deg)
                LATITUDE:    0.0    , $  ; Geo Eclip Lot of LOS (deg)
                PHOTOMETRY:  0.0      }  ; Observed MJy/sr

       defined = 1 

    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {RZODI}, Nstruct )
end

;---------------------------------------------------------------------
; Function MK_ZDATA
;
; Function to generat a zmodel input dataset for a set of LOS, Wavelength, 
; and Time.
;
; Written By: 	BA Franz
; Date:		23 January 1995
;
; Usage:
;     zdata = MK_ZDATA(lambda,day,lon,lat)
;
; Inputs:
;     lambda - wavelength(s) in um 
;     day    - 1990 day number(s), where 1.0 = 1 Jan 1990
;     lon    - ecliptic longitude(s)
;     lat    - ecliptic latitude(s)
;
;     note: lambda, day, lon, lat can be any mixture of scalars or arrays,
;           but all arrays must be of same length.
;
; Outputs:
;
;     zdata - data structure for input into ZMODEL or one of its kernel
;             functions (ZFUNC, ZTFUNC, ZPCFUNC, DZDTFUNC)
;
;     where:    
;         zdata(i).wave_len   = input wavelength (um)
;         zdata(i).longitude  = input ecliptic longitude (deg)
;         zdata(i).latitude   = input ecliptic latitude (deg)
;         zdata(i).day1990    = input day number (1.0 = 1 Jan 1990 00:00)
;
;---------------------------------------------------------------------
function mk_zdata,lambda,day,lon,lat

nwave = n_elements(lambda)
nday  = n_elements(day)
nlon  = n_elements(lon)
nlat  = n_elements(lat)

npts  = max([nwave,nday,nlon,nlat])

if ( (nwave ne 1 and nwave ne npts) or $
     (nday  ne 1 and nday  ne npts) or $
     (nlon  ne 1 and nlon  ne npts) or $
     (nlat  ne 1 and nlat  ne npts) ) then begin
    message,'Inputs must be scalars or arrays of equal dimension.',/info
endif

data = zodi_str(npts)
data.wave_len  = lambda
data.day1990   = day
data.longitude = lon
data.latitude  = lat

return,data
end
