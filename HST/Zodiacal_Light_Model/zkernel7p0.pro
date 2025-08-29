;============================================================================
;============================================================================
;
;                               ZKERNEL
;
;    DIRBE Kernel for use by ZMODEL to Fit the Zodiacal Dust Cloud
;
;                       Last Revision: 17 June 1996  (mark 6p8)
;
;============================================================================
;============================================================================

;
;----------------------------------------------------------------------------
; Function GAUSSINT
;
; Returns Gauss-Laguerre abscissae and weights for 24 or 50 point gaussian 
; quadrature.
;
; Written By: BA Franz, ARC, 9/94
; Quadrature abscissae and weights computed by W. T. Reach
;
;----------------------------------------------------------------------------
pro gaussint,a,b,grid,wts,numpts=numpts

if (n_elements(numpts) eq 0) then numpts = 24

;
; Set-up static storage space
common gaussint_cmn,init,grid24,wts24,grid25,wts25,grid50,wts50

;
; Initialize abscissae and weights if first call
;
if (n_elements(init) eq 0) then begin

    init = 1

    grid24 = double( $ 
         [0.0567047755,0.2990108986,0.7359095554,1.3691831160, $
          2.2013260537,3.2356758036,4.4764966151,5.9290837627, $
          7.5998993100,9.4967492209,11.629014912,14.007957977, $
          16.647125597,19.56289801, 22.775241987,26.308772391, $
          30.194291163,34.471097572,39.190608804,44.422349336, $
          50.264574993,56.864967174,64.466670616,73.534234792] $
                   )/73.534234792

    wts24  = double( $ 
         [0.1455497377,0.3393497722,0.5347365921,0.7322248725, $
          0.9326159014,1.1367925904,1.3457293379,1.560519046,  $
          1.7824092263,2.0128491498,2.2535525026,2.5065825127, $
          2.7744704429,3.0603848697,3.3683805666,3.7037765832, $
          4.0737527888,4.4883345170,4.9621093140,5.5174318658, $
          6.19095,7.0486426246,8.2265276593,10.0758794961]     $
                   )/73.534234792

    grid50 = double( $
        [0.0286305183,0.1508829357,0.3709487815,0.6890906999, $
         1.1056250235,1.6209617511,2.2356103759,2.9501833666, $
         3.7653997744,4.6820893876,5.7011975748,6.8237909098, $
         8.0510636694,9.3843453083,10.825109031,12.374981608, $
         14.035754599,15.809397197,17.698070933,19.704146535, $
         21.830223306,24.079151444,26.454057841,28.958376011, $
         31.595880956,34.370729963,37.287510610,40.351297573, $
         43.567720270,46.943043991,50.484267963,54.199244880, $
         58.096828017,62.187054175,66.481373878,70.992944826, $
         75.737011547,80.731404802,85.997211136,91.559690412, $
         97.449565614,103.70489123,110.37385880,117.51919820, $
         125.22547013,133.61202792,142.85832548,153.26037197, $
         165.3856433,180.69834370] $
                   )/180.69834370

    wts50  = double( $
        [0.0734786263,0.1711133062,0.2690583985,0.3672778840, $
         0.4658590232,0.5648992928,0.6644999757,0.7647657788, $
         0.8658052545,0.9677314398,1.0706625889,1.1747230029, $
         1.2800439495,1.3867646993,1.4950336890,1.6050098453, $
         1.7168640913,1.8307810796,1.9469611911,2.0656228595, $
         2.1870052872,2.3113716402,2.4390128272,2.5702520038, $
         2.7054499706,2.8450116969,2.9893942558,3.1391165624, $
         3.2947714220,3.4570405806,3.6267137126,3.8047126447, $
         3.9921226303,4.1902332693,4.4005928365,4.6250816069, $
         4.8660126478,5.1262732767,5.4095283377,5.7205203736, $
         6.0655271405,6.4530854743,6.8951889633,7.4093807809, $
         8.0226692310,8.7795282602,9.7603031266,11.131390723, $
         13.324267692,18.114146002] $
                   )/180.69834370

endif

;
; Select grid length based on user input
;
case numpts of
    50   : begin
             grid = grid50
             wts  = wts50
           end
    else : begin
             grid = grid24
             wts  = wts24
           endelse
endcase

;
; Scale grid to user limits
;
grid = temporary(grid) * (b-a) + a
wts  = temporary(wts)  * (b-a)

return
end

;
;-------------------------------------------------------------------------------
;
; Procedure SimpInt
;
; Simpson Integration abscissae and weights as per HTF NEWESTZKERNEL0
; Necessary for toroidal bands which are very narrow 
; Added: JP 07 may 1996
;
pro simpint,a,b,grid,wts,stepsize=stepsize

; a is assumed to be zero!

if (n_elements(stepsize) eq 0) then stepsize = 0.025

h = stepsize
tsteps = round(b/h) + 1

grid = findgen(tsteps) * h

wts  = fltarr(tsteps)
wts(0) = .333333 * h
for i = 1,tsteps-1,2 do wts(i) = 1.333333 * h
for i = 2,tsteps-2,2 do wts(i) = 0.666667 * h
wts(tsteps-1) = .333333 * h

return
end

;
;-------------------------------------------------------------------------------
;
; Procedure EARTHSUN
;
; Returns Solar Elongation and Earth Position from Input Day and Lon, Lat
;
; Written By:  BA Franz, WT Reach
;
; Inputs:
;     Day   - 1990 Day Number (1 = 01 Jan 1990)
;     Lon   - Ecliptic Lon of Line-of-Sight (Deg)
;     Lat   - Ecliptic Lat of Line-of-Sight (Deg)
;-------------------------------------------------------------------------------
pro earthsun,day,lon,lat,SolElong,Earth_Dis,Earth_Lon,Earth_Mean_Lon

pi2   = 2*!pi
d2r   = !pi/180.
eccen = 0.01671254

lambda_solar = (-80.598349 + 0.98564736 * day + $
	       1.912 * cos(pi2/365.25 * (day-94.8))) * d2r

mean_anomaly = (356.637087 + 0.98560028 * day)*d2r mod pi2

Earth_Dis = (1.-eccen^2)/(1.+eccen*cos(mean_anomaly))  ; (AU)

Earth_Lon = -(!pi-lambda_solar) mod (2.*!pi)           ; (rad)

SolElong  = acos( cos(lat*d2r) * cos(lon*d2r - lambda_solar) )

Earth_Mean_Lon = ((99.403445 + 0.98564736*day)*d2r) mod pi2

return
end


;
;-------------------------------------------------------------------------------
;
; Function COLCORR
;
; Returns DIRBE color correction coefficients for a modified BB source.
;
; Performs a table look-up of temperature versus correction.  The color
; corrections are generated from the DIRBE Explanatory Supplement tables 
; from Appendix B, for a black body with an emissivity index of 0.  On first 
; call, the Exp. Supp. data is interpolated to an even temperature grid
; and stored in an array as described below:
;
;	Table(itemp,ival,idet)	
;
; 	itemp = temperature index
;	ival  = value index (0=Color Correction, 1=dCdT)
;	idet  = detector index (0=1.25um, 9=240um)
;
; Written By:  BA Franz, 10/95
;
; Warning:  Only applies to standard DIRBE wavelengths.
;
;-------------------------------------------------------------------------------

function colcorr,det,t,dcdt

;
; Create static space for colcor correction tables and supporting data
;
common colcor_cmn,init,temp,table,tmax,tmin,dt,nt

;
; If first call, initialize color correction tables
;
if (n_elements(init) eq 0) then begin

    init = 1

    ;
    ; Read Exp. Supp. Table for BB at P=0
    ;

;    infile = 'colcorr.tab'
;    print,'Reading color correction file: ',infile

    fmt    = '(f,f,f,f,f,f,f,f,f,f,f)'
    readcol,infile,form=fmt,intemp,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
    incc   = [[c1],[c2],[c3],[c4],[c5],[c6],[c7],[c8],[c9],[c10]]

    ;
    ; Define Interpolation Range
    ;
    tmax =  1000.
    tmin =  10.
    dt   =  5.
    nt   =  (tmax - tmin)/dt + 1.
    temp =  findgen(nt)*dt + tmin

    ;
    ; Build Table of CC and dCdT per Det
    ;
    table = fltarr(nt,2,10)
    for i=0,9 do begin
        cc  = spline(intemp,incc(*,i),temp)
        cc2 = spline(intemp,incc(*,i),temp+dt)
        table(*,0,i) = cc
        table(*,1,i) = (cc2-cc)/dt
    endfor

endif

;
; Compute CC and dCdT at Each Temperature in t
;
it      = fix((t - tmin)/dt) >  0  <  (nt-1)
delta_t = (t-temp(it))
dcdt    = table(it,1,det)
cc      = table(it,0,det) + dcdt * delta_t 

return,cc
end


;
;-------------------------------------------------------------------------------
;+
; NAME:
;    ZPLANCK
;
; PURPOSE:
;    IDL function to return the spectral radiance of a blackbody,
;    i.e. the Planck curve, in units of W/cm^2/sr.  This is 
;    nu*I_nu (=lambda*I_lambda).  The 
;    blackbody temperature and either frequency (in icm or GHz)
;    or wavelength (in microns) are inputs to the function.  The
;    routine also optionally returns the derivative with respect to 
;    temperature, in units of W/cm^2/sr/K.
;
; CALLING SEQUENCE:
;    result = planck (temperature, nu_or_lambda, [dBbT], [units=<units>])
;
; INPUTS: 
;    T              float scalar    Temperature of blackbody, in K.
;    nu_or_lambda   float [array]   Frequency or wavelength at which to
;                                   calculate spectrum.  Units are as
;                                   specified with "units" keyword.
;
; OPTIONAL INPUT:
;    units          keyword string  'Microns', 'icm', or 'GHz'  to identify 
;                                   units of nu_or_lambda.  Only first 
;                                   character is required.  If left out,
;                                   default is 'microns'.
;
; OUTPUTS:
;    planck         float [array]   Spectral radiance at each wavelength.
;                                   Units are W/cm^2/sr.
;
; OPTIONAL OUTPUT:
;    dB/dT          float [array]   (Optional) derivative of Planck 
;                                   with respect to temperature.  Units 
;                                   are W/cm^2/sr/K
;
; SIDE EFFECTS:
;    None known.
;
; PROCEDURE:
;    Identifies units using the "units" keyword, then converts the 
;    supplied independent variable into microns to evaluate the Planck 
;    function.  Uses Rayleigh-Jeans and Wien approximations for the low-
;    and high-frequency end, respectively.
;
; EXAMPLE:
;    To produce a 35 K spectrum at 2,4,6,8,...,100 microns:
;
;       wavelength = 2. + 2.*findgen(50) 
;       temp  = 35.
;       blackbody = planck (wavelength, temp, units='micron')
;
;    One could also get back the derivative by including it in the call:
;       blackbody = planck (wavelength, temp, deriv, units='m')
;
; RESTRICTIONS:
;    Temperature must be given as a real*4 input, NOT an integer!  
;    Routine also gives incorrect results for T < 1 microkelvin (so sue me).
;
; REVISION HISTORY:
;    Written by Rich Isaacman, General Sciences Corp.  17 April 1991
;    Allow array of Temperatures and single Lambda
;    Rename to avoid confusion  6/95
;
; REFERENCE:
;    See Allen, Astrophysical Quantities for the Planck formula.
;-
function zplanck, T, lambda, dBdT
;
; Define some necessary constants
;
c = 2.9979e14                            ;Speed of light (micron/sec)
h = 6.6262e-34                           ;Planck constant
k = 1.3807e-23                           ;Boltzmann constant
stephan_boltzmann = 5.6697e-12           ;Stephan-Boltzmann constant
coeff = 15. / !PI^5 * stephan_boltzmann  ;Appropriate combination
;
;  Introduce dimensionless variable chi, used to check whether we are on 
;  Wien or Rayleigh Jeans tails
;

chi = h * c / (lambda * k * T)
val = CHI*0
dBdT = CHI*0
;  Start on Rayleigh Jeans side
;
rj = where (chi lt 0.001)
if rj(0) ne -1 then begin
    val(rj) = coeff * T(rj)^4 * chi(rj)^3
    dBdT(rj) = val(rj)/T(rj) 
endif
;
;  Now do nonapproximate part
;
exact = where (chi ge 0.001 and chi le 50)
if exact(0) ne -1 then begin
    chi_ex = chi(exact)
    val(exact) = coeff * T(exact)^4 * chi_ex^4 / (exp(chi_ex) - 1.)
    dBdT(exact) = val(exact) * chi_ex / T(exact) / (1. - exp(-chi_ex))
endif
;
;  ...and finally the Wien tail
;
wien = where (chi gt 50.)
if wien(0) ne -1 then begin
    chi_wn = chi(wien)
    val(wien) = coeff * T(wien)^4 * chi_wn^4 * exp(-chi_wn)
    dBdT(wien) = val(wien) * chi_wn / T(wien) / (1. - exp(-chi_wn))
endif
;
if (n_elements(chi) eq 1) then val=val(0)
if (n_elements(chi) eq 1) then dBdT=dBdT(0)

return, val
;
end


;
;--------------------------------------------------------------
; Function THERMFUNC
;
; Returns the thermal component of the zodiacal dust source 
; function.
;
; Written By: BA Franz, ARC, 9/95
;
; Inputs: 
;   det      - DIRBE detector number (0-9) or -1 if non-DIRBE
;   Lambda   - Wavelength in um
;   R        - sqrt(x^2+y^2+z^2) (AU)
;   a        - Model parameters
;
; Outputs:
;   THERMFUNC - Brightness per particle in MJy/sr
;   df        - Partial derivative of source function wrt each 
;               parameter
;
;--------------------------------------------------------------
function thermfunc,det,Lambda,R,a,df,no_colcorr=no_colcorr

    if (n_params(0) gt 4) then want_partials=1 else want_partials=0
 
    d2r  = !pi/180.
    eps  = double(1e-20)

    ;
    ; Thermal Parameters
    ;
    To1      = a(0)
    To2      = a(1)
    Kappa    = a(2)
    Delta    = a(3) ; Temperature power law component 
    ; Remember: The temperature falls off as R^Delta, so the temperature for each iteration
    ; along the line of site is different!

    ;
    ; Conversion from W/cm^2/Sr (Nu Bnu) to MJy/Sr (Bnu)
    ; Lambda/3e14 Hz^-1 * 1e4 cm^2/m^2 * 1e20 (MJy)/(W/m^2)
    ;
    cfact = Lambda / 3.0e-10  

    ;
    ; First Temperature Component
    ;
    Temp1    = To1/R^Delta
    Bnu1     = zPlanck(Temp1,Lambda,dBdT1) * cfact
    dBdT1    = dBdT1 * cfact
    if (det ge 0 and not keyword_set(no_colcorr)) then begin
        print,'Doing color correction!'
        CCtherm1 = colcorr(det,temp1,dCdT1)
    endif else begin
        CCtherm1 = 1.0        
        dCdT1    = 0.0
    endelse

    ;
    ; Second Temperature Component
    ;
    if (Kappa ne 0.0) then begin
        Temp2    = To2/R^Delta
        Bnu2     = zPlanck(Temp2,Lambda,dBdT2) * cfact    
        dBdT2    = dBdT2 * cfact
        if (det ge 0 and not keyword_set(no_colcorr)) then begin
            CCtherm2 = colcorr(det,temp2,dCdT2)
        endif else begin
            CCtherm2 = 1.0        
            dCdT2    = 0.0
        endelse
    endif else begin
        Temp2    = 0.0
        Bnu2     = eps
        dBdT2    = 0.0
        CCtherm2 = 1.0
        dCdT2    = 0.0
    endelse

    ;
    ; Total Thermal Contribution
    ;
    Therm    = Bnu1*CCtherm1  + Kappa * Bnu2*CCtherm2
    dBdT     = dBdT1 + Kappa * dBdT2

    ;
    ; Calculate Derivatives of Varying Parameters if Requested
    ; --------------------------------------------------------
    if (want_partials) then begin

        npts = n_elements(R)
        npar = n_elements(a)
        if (n_elements(df) ne npts*npar) then df = dblarr(npts,npar)

        dT1 =  dBdT1*CCTherm1 + Bnu1*dCdT1
        dT2 = (dBdT2*CCTherm2 + Bnu2*dCdT2)*Kappa
        
        ; To1
        df(*,0) = dT1/R^Delta
 
        ; To2
        df(*,1) = dT2/R^Delta
 
        ; Kappa
        df(*,2) = Bnu2*CCtherm2

        ; Delta
        df(*,3) = -1.*alog(R)*(Temp1*dT1 + Temp2*dT2)
 
    endif

return,Therm
end    


;
;-------------------------------------------------------------------------------
; IDL Function PHASEFUNC.PRO
;
; Trial Phase Function for ZKERNEL
;
; Primary Author:  JP, using BAF template  Dec 1995
;
; Contributors:   
;
; Inputs:
;        x - array of scattering angles (radians)
;        a - Model parameters
;
; Optional Inputs:
;
; Outputs:
;    phi   - array of phase function
;     pder - partial derivatives of phase function w.r.t each parameter
;
; Optional Outputs:
;
; Last modification:   18 March 1996 JP -- three parameter PF
; 01 September RA -- Sum of 3 HG functions if (a.length eq 6)
;
;-------------------------------------------------------------------------------
function phasefunc,x,a,pder

if (n_elements(a) eq 3) then begin 
  C0    = a(0)
  C1    = a(1)
  ka    = a(2)
  ;
  ; Evaluate the function
  ;
  pisq = !pi^2
  picu = !pi^3
  q0 = 2.
  q1 = !pi
  qk = (exp(ka*!pi) + 1.)/(ka*ka + 1)
  v = (2.*!pi)*( C0*q0 + C1*q1 + qk)
  u = C0 + C1*x + exp(ka*x)
  phi = u/v
  norm_phi = phi
  ;
  ; Evaluate the partial derivatives
  ;
      npts = n_elements(x)
      pder = dblarr(npts,3)
      pder(*,0) =  (1.d0 - phi*4.*!pi)/v
      pder(*,1) =  (x - 2.d0*!pi*phi*q1)/v
      pder(*,2) =  (x*exp(ka*x) - 2.d0*!pi*phi*( (ka^2+1)*!pi*exp(ka*!pi) $
	  - 2.d0*ka*(exp(ka*!pi)+1))/((ka*ka+1.)^2))/v
endif

; change for any number of HG functions...
if (n_elements(a) gt 3) then begin

  phi = hong_phase_func(x,a)

endif

return,phi
end

;
;--------------------------------------------------------------
; Function SCATTFUNC
;
; Returns the scattered light component of the zodiacal dust 
; source function.
;
; Written By: BA Franz, ARC, 11/95
; Modified By: JP Dec 1995
;
; Inputs: 
;   det      - DIRBE detector number (0-9) or -1 if non-DIRBE
;   Lambda   - Wavelength (µm)
;   LOS      - Distance from Earth (AU)
;   R        - sqrt(x^2+y^2+z^2) (AU)
;   Re       - Distance from Earth to Sun (AU)
;   SolElong - Solar elongation (radians)
;   a        - Model parameters
;
; Outputs:
;   SCATTFUNC - Brightness per particle in MJy/sr
;   df        - Partial derivative of source function wrt each 
;               parameter
;
;--------------------------------------------------------------
function scattfunc,phase_type,det,Lambda,LOS,R,Re,SolElong,a,solar_irr=solar_irr

    if (n_params(0) gt 7) then want_partials=1 else want_partials=0
 
    d2r  = !pi/180.
    eps  = double(1e-20)

    ;
    ; Solar Flux at 1 AU from Sun (courtesy W. T. Reach)
    ;
    SolFlux1AU = [ 2.3405606e+08, $
                   1.2309874e+08, $
                   64292872.,     $
                   35733824.,     $
                   5763843.0,     $
                   1327989.4,     $
                   230553.73,     $
                   82999.336,     $
                   42346.605,     $
                   14409.608 ]
    if (n_elements(solar_irr) gt 0) then begin
        SolFlux = solar_irr / R^2
    endif else if ((det ge 0) and (phase_type eq 'skysurf')) then begin
        SolFlux = solar_sp(Lambda) / R^2
    endif else if ((det ge 0) and (phase_type eq 'kelsall')) then begin
        SolFlux = SolFlux1AU(det) / R^2
;        SolFlux = solar_sp(Lambda) / R^2
    endif else begin
        SolFlux = 0.0
    endelse

    if (det le 2) then begin
      ;
      ; Evaluate Phase Function
      ;
      PhaseAng   = asin( Re/R*sin(SolElong) < 1.0 )
      itest      = 1.0*(los ge Re*cos(SolElong))
      ScatAng    = (1.-itest)*PhaseAng + (itest)*(!pi - PhaseAng)
      if (n_elements(a) eq 9) then aaa        = a(3*det:3*det+2) ; original
      if (n_elements(a) mod 2 eq 0) then aaa  = a ; when using HG functions
      Phi        = phasefunc(ScatAng,aaa,dphi)

      ;
      ; Total Scattered Light Contribution
      ;
      CCscat     = 1.D0
      Scatt      = SolFLux * CCscat * Phi
    endif else begin
      Scatt      = 0.0d0
    endelse

    ;
    ; Calculate Derivatives of Varying Parameters if Requested
    ; --------------------------------------------------------
;    if (want_partials) then begin
;
;        npts = n_elements(R)
;        npar = n_elements(a)
;        if (n_elements(df) ne npts*npar) then df = dblarr(npts,npar) $
;                                         else df = temporary(df)*0.0
;
;        if (det le 2) then begin
;        ; PF_C0
;        df(*,3*det) = SolFlux * CCscat * dphi(*,0)
;
;        ; PF_C1
;        df(*,3*det+1) = SolFlux * CCscat * dphi(*,1)
;
;        ; PF_Ka
;        df(*,3*det+2) = SolFlux * CCscat * dphi(*,2)
;
;	endif
;
;    endif

return,Scatt
end    
;
;--------------------------------------------------------------
; Function ZSRCFUNC
;
; Returns the zodiacal dust source function, the brightness per 
; dust grain.
;
; Written By: BA Franz, ARC, 9/95
;
; Inputs: 
;   Scatt    - Scattered light contribution
;   Therm    - Thermal emission contribution
;   a        - Model parameters (Emissivities/Albedos)
;
; Outputs:
;   ZSRCFUNC  - Brightness per particle in MJy/sr
;   df        - Partial derivative of source function wrt each 
;               parameter
;
;--------------------------------------------------------------
function zsrcfunc,det,Scatt,Therm,a,df,dScatt,dTherm,phase_type

    if (n_params(0) gt 3) then want_partials=1 else want_partials=0
 
    npar = n_elements(a)

    ; 
    ; Albedo per Det
    ;
    AlbedoDet = [a(0:3),dblarr(6)+a(3)]

    ;
    ; Thermal Emissivity per Det
    ;
    EmissDet = [dblarr(2)+1.,a(4:11)] 

    ;
    ; Set Albedo and Emissivity for This Detector
    ; -------------------------------------------
    if (det ge 0) then begin
        Albedo = AlbedoDet(det)
        
        ; In put_zpar.pro, the emissitivity for skysurf is put into the slow for detector = 2 (which corresponds to 3.5 micron)
        ; So if using the 'skysurf' mode, use that emissivity
        if phase_type eq 'skysurf' then begin
          Emiss  = EmissDet(2)
          
        endif else if phase_type eq 'kelsall' then begin
          Emiss  = EmissDet(det)
          
        endif
        
    endif else begin
        Albedo   = 0.
        Emiss    = 1.
    endelse
    
    ;
    ; Calculate Total Source per LOS Element
    ; ------------------------------------------
    Source = Albedo * Scatt + Emiss * (1.0-Albedo) * Therm

    ;
    ; Calculate Derivatives of Varying Parameters if Requested
    ; --------------------------------------------------------
    if (want_partials) then begin

        npts = n_elements(Therm)
        npar = n_elements(a)
        if (n_elements(df) ne npts*npar) then  $
            df = dblarr(npts,npar)             $
        else                                   $
            df = temporary(df)*0.0             

        dScatt = Albedo
        dTherm = Emiss * (1.0-Albedo)

        ; Albedo
        sign1 = 1.0*(Albedo ge 0) - (Albedo lt 0)
        sign2 = 1.0*((1.-Albedo) ge 0) - ((1.-Albedo) lt 0)
        dA    = sign1*Scatt - sign2*Emiss*Therm
        df(*,(det < 3)) = dA

        ; Emissivity
        if (det ge 2) then $
            df(*,4 + (det-2)) = Therm/Emiss

    endif

return,Source
end    


;
;--------------------------------------------------------------
; Function ZCloud
;
; Returns the number density of the zodiacal cloud.
;
; Written By: BA Franz, ARC, 9/95
;
; Inputs: 
;   x,y,z - heliocentric cartesian coordinates
;   R     - sqrt(x^2+y^2+z^2)
;   a     - Model parameters
;
; Keyword Inputs:
;   FuncIndex - selector for latitudinal density distribution
;		0 : Good
;		1 : Modified Fan
;		2 : Widened Modified Fan
;		3 : Ellipsoid
;		4 : Sombrero
;
; Outputs:
;   ZCloud  - number density
;   df      - partial of density wrt each parameter
;
;--------------------------------------------------------------
function zcloud,x,y,z,R,a,df,FuncIndx=FuncIndx

if (n_params(0) gt 5) then want_partials=1 else want_partials=0

No = a(0)

if (No ne 0.0) then begin

    d2r = !pi/180.

    Alpha    = a(1)   
    Beta     = a(2) 
    Gamma    = a(3) 
    Mu       = a(4)
    Mu2      = a(5)
    Mu3      = a(6)
    Mu4      = a(7)
    Mu5      = a(8)
    Mu6      = a(9)
    Mu7      = a(10)
    Mu8      = a(11)
    Mu9      = a(12)
    Mu10     = a(13)

    Omega    = a(14) * d2r
    Incl     = a(15) * d2r
    Xo       = a(16)
    Yo       = a(17)
    Zo       = a(18)

    ;
    ; Some frequently used geometry terms
    ;
    sino = sin(Omega) 
    coso = cos(Omega)
    sini = sin(Incl)
    cosi = cos(Incl)

    ;
    ; Translate origin from Sun to cloud center.
    ;
    Xp  = X - Xo
    Yp  = Y - Yo
    Zp  = Z - Zo
    Rc  = sqrt( Xp^2 + Yp^2 + Zp^2 )

    ;
    ; Rotate into coord system of cloud to determine Z-height in cloud
    ; (Zc) and radial distance in cloud symetry plane (Rc).
    ; Perform rotation of Omega around Z-Axis to align X with cloud 
    ; line-of-nodes. Then, perform Incl rotation around X-Axis. 
    ;
    Zc   = sino*sini*Xp - coso*sini*Yp + cosi*Zp
    Zeta = abs(Zc/Rc)

    ;
    ; Compute Radial power-law parameter
    ; Mark6p2  sets this to Rc, not R    20 May 1996
    ; Mark6p7 allows the polynomial mods only beyond a cutoff radius
    ; Mark6p8 puts is back the way it was, but keeps insistence on positive
    ;         coeffs -- uses R, not Rc
    ;
    AlphaR = Alpha + ( Mu6^2*(R-1.) + Mu7^2*(R-1.)^2 + Mu8^2*(R-1.)^3 )
    ;protection against overflow 12May1996
    stuff = where(AlphaR gt 50.)
    if (stuff(0) ne -1) then AlphaR(stuff) = 50.

    ;
    ; Radial and Latitudinal Density Distribution
    ;
    case FuncIndx of

        ;
        ; John Good Density Function
        ;
        0 : begin
            Dc     = sqrt(Rc^2 - Zc^2)  
            ZoD    = abs(Zc/Dc)
            ZoD2G  = ZoD^Gamma
            Dens_Cloud_Vert = exp( -Beta * ZoD2G )
            Dens_Cloud_Rad  = 1./Dc^AlphaR
        end

        ;
        ; Modified Fan Density Function
        ;
        1 : begin
            Zeta2G  = Zeta^Gamma
            Dens_Cloud_Vert = exp( -Beta * Zeta2G )
            Dens_Cloud_Rad  = 1./Rc^AlphaR
        end

        ;
        ; Widened Modified Fan Density Function
        ;
        2 : begin
            GZR    = (0.5*Zeta^2/Mu) * (Zeta lt Mu) + $
                     (Zeta-0.5*Mu)   * (Zeta ge Mu)
            GZR2G  = GZR^Gamma
            Dens_Cloud_Vert = exp( -Beta * GZR2G )
            Dens_Cloud_Rad  = 1./Rc^AlphaR
        end

        ;
        ; Ellipsoidal Density Function
        ;
        3 : begin
            Bterm = (1. + (Beta*Zeta)^2)
            Dens_Cloud_Vert = 1./Bterm^Gamma
            Dens_Cloud_Rad  = 1./Rc^AlphaR
        end

        ;
        ; Cosine (Sombrero)
        ;
        4 : begin
            cosB   = sqrt(1-Zeta^2)
            cosB2G = cosB^Gamma
            Dens_Cloud_Vert = (1. + Beta * cosB2G)/(1. + Beta)
            Dens_Cloud_Rad  = 1./Rc^AlphaR
        end

        ;
        ; Modified Fan Density Function with Polynomial Z-height
        ;
        5 : begin
            aZ = abs(Z)
            BetaZ = Beta + Mu3*aZ + Mu4*aZ^2 + Mu5*aZ^3
            Dens_Cloud_Vert = exp( -BetaZ * Zeta )
            Dens_Cloud_Rad  = 1./Rc^AlphaR
        end

    endcase
    ;
    ; Total Density of Cloud Component
    ;
    Dens = No * Dens_Cloud_Vert * Dens_Cloud_Rad

    ;
    ; Partials of Density wrt Parameters
    ;
    if (want_partials) then begin

        npts = n_elements(x)
        npar = n_elements(a)
        if (n_elements(df) ne npts*npar) then df = dblarr(npts,npar)

        sZc  = 1.0*(Zc ge 0.) - (Zc lt 0.)

        ; No
        df(*,0) = Dens/No

        case FuncIndx of
            0    : lnR = alog(Dc)
            else : lnR = alog(Rc)
        endcase

        ; Alpha
        df(*,1)  = -Dens*lnR

        ; Mu6
        df(*,9) = -Dens*lnR*(R-1.) *2.*Mu6

        ; Mu7
        df(*,10) = -Dens*lnR*(R-1.)^2 *2.*Mu7

        ; Mu8
        df(*,11) = -Dens*lnR*(R-1.)^3 *2.*Mu8


        case FuncIndx of

            ;
            ; John Good Density Function
            ;
            0 : begin

                Zsqr = (1.-Zeta^2)

                dlnf = AlphaR*Zeta/Zsqr - $
                       Beta*Gamma*Zeta^(Gamma-1)*Zsqr^(-Gamma/2-1)

                ; Beta
                df(*,2) = -Dens*ZoD2G

                ; Gamma
                df(*,3) = -Beta*Dens*ZoD2G*alog(ZoD)

            end

            ;
            ; Modified Fan Density Function
            ;
            1 : begin

                dlnf = -Beta*Gamma*Zeta^(Gamma-1)

                ; Beta
                df(*,2) = -Dens*Zeta2G

                ; Gamma
                df(*,3) = -Beta*Dens*Zeta2G*alog(Zeta)

            end

            ;
            ; Widened Modified Fan Density Function
            ;
            2 : begin

                dlnf = -Beta*Gamma*GZR^(Gamma-1) * $
                    ( (Zeta/Mu) * (Zeta lt Mu) + 1.* (Zeta ge Mu) )

                ; Beta
                df(*,2) = -Dens*GZR2G

                ; Gamma
                df(*,3) = -Beta*Dens*GZR2G*alog(GZR)

                ; Mu
                df(*,4) = Beta*Dens*0.5* Gamma * GZR^(Gamma-1)* $
                             ( (Zeta/Mu)^2 * (Zeta lt Mu) + (Zeta ge Mu) )

            end

            ;
            ; Ellipsoidal Model
            ;
            3 : begin

                dlnf = -2*Gamma*Beta^2*Zeta/(1+(Beta*Zeta)^2)

                ; Beta
                df(*,2) = -2.*Gamma/Beta*Dens*(1.-1./Bterm)

                ; Gamma
                df(*,3) = -Dens*alog(Bterm)

            end

            ;
            ; Cosine Model
            ;
            4 : begin

                dlnf = -Gamma*Zeta*Beta/(1+Beta)/Dens_Cloud_Vert*cosB2G/cosB^2

                ; Beta
                df(*,2) = Dens/Dens_Cloud_Vert* $
                              (-Dens_Cloud_Vert + cosB2G)/(1+Beta)

                ; Gamma
                df(*,3) = Dens/Dens_Cloud_Vert* $
                              (Beta*cosB2G/(1+Beta)*alog(cosB))
            end

            ;
            ; Modified Fan Density Function
            ;
            5 : begin
              
                dlnf = -BetaZ

                ; Beta
                df(*,2) = -Dens*Zeta

                ; Mu3
                df(*,6) = -Dens*Zeta*aZ

                ; Mu4
                df(*,7) = -Dens*Zeta*aZ^2

                ; Mu5
                df(*,8) = -Dens*Zeta*aZ^3

            end

        endcase

        ; Omega
        df(*,14) = Dens*dlnf*sZc/Rc*(coso*sini*Xp+sino*sini*Yp)*d2r
 
        ; Incl
        df(*,15) = Dens*dlnf*sZc/Rc*(sino*cosi*Xp-coso*cosi*Yp-sini*Zp)*d2r

        ; Xo
        df(*,16) = Dens*(-dlnf*sZc/Rc*sini*sino + (dlnf*Zeta+AlphaR)*Xp/Rc^2)
 
        ; Yo
        df(*,17) = Dens*(dlnf*sZc/Rc*coso*sini + (dlnf*Zeta+AlphaR)*Yp/Rc^2)

        ; Zo
        df(*,18) = Dens*(-dlnf*sZc/Rc*cosi + (dlnf*Zeta+AlphaR)*Zp/Rc^2)
 
    endif

endif else Dens = x*0

return,Dens
end


;
;--------------------------------------------------------------
; Function MigBand
;
; Returns the number density of a dust band.
;
; The model is a parameterization of the migrating band concept
; of W. T. Reach, as provided by Reach.
;
; Written By: BA Franz, ARC, 9/95
; Modified By: JP 26 May 1996 based on HTF R cutoff idea
;              changes meaning of dR, and makes p2r and Vr useless
; Modified By: JP 11 June 1996 to protect against floating underflow
;
; Inputs: 
;   x,y,z - heliocentric cartesian coordinates
;   R     - sqrt(x^2+y^2+z^2)
;   a     - Model parameters
;
; Outputs:
;   migband - number density
;   df      - partial of density wrt each parameter
;
;--------------------------------------------------------------
function migband,x,y,z,R,a,df

if (n_params(0) gt 5) then want_partials=1 else want_partials=0

No = a(0)

if (No ne 0.0) then begin

    d2r = !pi/180.

    Dz       = a(1) * d2r
    Dr       = a(2)
    Ro       = a(3)
    Vi       = a(4)
    Vr       = a(5)
    Pi       = a(6)
    Pr       = a(7)
    P2r      = a(8)
    Omega    = a(9)  * d2r
    Incl     = a(10) * d2r
    Xo       = a(11)
    Yo       = a(12)
    Zo       = a(13)

    sino = sin(Omega)
    coso = cos(Omega)
    sini = sin(Incl)
    cosi = cos(Incl)

    ;
    ; Translate origin from Sun to cloud center.
    ;
    Xp  = X - Xo
    Yp  = Y - Yo
    Zp  = Z - Zo

    ;
    ; Rotate into coord system of cloud to determine Z-height in cloud
    ; (Zc) and radial distance in cloud symetry plane (Rc).
    ; Perform rotation of Omega around Z-Axis to align X with cloud 
    ; line-of-nodes. Then, perform Incl rotation around X-Axis. 
    ;
    Zc = sino*sini*Xp - coso*sini*Yp + cosi*Zp
    Rc = sqrt( Xp^2 + Yp^2 + Zp^2 )

    Zeta   = abs(Zc/Rc)
    ZDz    = Zeta / Dz
    RDr    = Rc / Dr
    ViTerm = 1. +  ZDz^Pi  / Vi
; protect against underflow in exponential terms
    WtTerm = 1. + 0.*RDr
    Dens   = 0.0 * Rc 
    arg1   = RDr^20
    arg2   = ZDz^6
    ican   = where ((arg1 le 86) and (arg2 le 86))
    if (ican(0) ne -1) then begin
         WtTerm(ican) = 1. -  exp(-arg1(ican))
         Dens(ican) = No * (Ro/Rc(ican))^Pr * exp(-arg2(ican))  $
                 * Viterm(ican) * WtTerm(ican)
    endif
    ican   = where((arg2 le 86) and (arg1 gt 86))
    if (ican(0) ne -1) then $
         Dens(ican) = No * (Ro/Rc(ican))^Pr * exp(-arg2(ican))  $
                 * Viterm(ican)      ;WtTerm should = 1.
;
    if (want_partials) then begin

        npts = n_elements(x)
        npar = n_elements(a)
        if (n_elements(df) ne npts*npar) then df = dblarr(npts,npar)

        sZ   = 1.0*(Zc ge 0.) - (Zc lt 0.)
        dlnf = (-6. * ZDz^5 + Pi/Vi * ZDz^(Pi-1) / ViTerm) / Dz

        ; No
        df(*,0) = Dens/No       ; No

        ; Dz
	df(*,1) = -Dens*ZDz*dlnf*d2r

        ; Dr   --- new
; now protect against underflow
	df(*,2) = 0. * Rc
	arg1 = ZDz^6
        arg2 = RDr^20/Dr
        ican = where((arg1 le 86) and (arg2 le 86) and (Dens ne 0.0))
	if (ican(0) ne -1) then $
        df(ican,2) =  No * (Ro/Rc(ican))^Pr * exp(-arg1(ican))*ViTerm(ican)  $
                 * (-20 * RDr(ican)^20 * exp(-arg2(ican)) )

        ; Ro

        ; Vi
        df(*,4) = -Dens*Vi^(-2) * ZDz^Pi / ViTerm

        ; Vr

        ; Pi
        df(*,6) = Dens * ZDz^Pi * alog(ZDz) / ViTerm / Vi

        ; Pr
        df(*,7) = Dens*alog(Ro/Rc)

        ; P2r

        ; Omega
        df(*,9) = Dens*dlnf*sZ/Rc*(coso*sini*Xp + sino*sini*Yp)*d2r

        ; Incl
        df(*,10) = Dens*dlnf*sZ/Rc*(sino*cosi*Xp - coso*cosi*Yp - sini*Zp)*d2r

        ; Xo
        df(*,11) = Dens*(-dlnf*sZ/Rc*sini*sino + (dlnf*Zeta+Pr)*Xp/Rc^2)
 
        ; Yo
        df(*,12) = Dens*(dlnf*sZ/Rc*coso*sini + (dlnf*Zeta+Pr)*Yp/Rc^2)

        ; Zo
        df(*,13) = Dens*(-dlnf*sZ/Rc*cosi + (dlnf*Zeta+Pr)*Zp/Rc^2)

    endif

endif else Dens = x*0

return,Dens
end


;
;--------------------------------------------------------------
; Function SolRing
;
; Returns the number density of the Earth's resonant dust ring.
; The model is a gaussian toroidal ring with Earth-leading and
; Earth-trailing gaussian blobs.  It is a parameterization of
; W. T. reach based on numerical images published by S. Dermott.
;
; Written By: BA Franz, ARC, 9/95
; Modified By: JP, 01-Apr-1996 to allow 2nd radius thickness to ring
;                  10-Apr-1996 to change R and Z laws in Dens_SR
;                  24-May-1996 to change R law back to r^-2 from r-4
;                  11-Jun-1996 to avoid unecessary component calcs
;                              and protect against underflow
;                  14-Jun-1996 speed up omega and incl derivs.
;                       
;                            
;
; Inputs: 
;   x,y,z - heliocentric cartesian coordinates
;   R     - sqrt(x^2+y^2+z^2)
;   Theta - Earth mean longitude
;   a     - Model parameters
;
; Outputs:
;   solring - number density
;   df      - partial of density wrt each parameter
;
;--------------------------------------------------------------
function solring,x,y,z,R,Theta,a,df

if (n_params(0) gt 6) then want_partials=1 else want_partials=0

SR_No = a(0)
LB_No = a(4)
TB_No = a(10)

if (SR_No ne 0.0 or LB_No ne 0.0 or TB_No ne 0) then begin

    d2r = !pi/180.

    ; 
    ; Set Parameters
    ;
    ; Ring
    SR_R        = a(1)
    SR_dR       = a(2)
;    SR_dR2      = a(21)                ;JP -- outer thickness
    SR_dZ       = a(3)
    ; Leading Blob
    LB_R        = a(5)
    LB_dR       = a(6)
    LB_Theta    = a(7) * d2r
    LB_dTheta   = a(8) * d2r
    LB_dZ       = a(9)
    ; Trailing Blob
    TB_R        = a(11)
    TB_dR       = a(12)
    TB_Theta    = a(13) * d2r
    TB_dTheta   = a(14) * d2r
    TB_dZ       = a(15)
    ; Symmetry plane
    Omega       = a(16) * d2r
    Incl        = a(17) * d2r
    Xo          = a(18)
    Yo          = a(19)
    Zo          = a(20)

    ;
    ; Some frequently used geometry terms
    ;
    sino = sin(Omega) 
    coso = cos(Omega)
    sini = sin(Incl)
    cosi = cos(Incl)

    ;
    ; Compute Density
    ;
    SR_Z     = sino*sini*X - coso*sini*Y + cosi*Z
    dTheta   = (atan(Y,X) - Theta)
;
;  JP changes --no more use of SR_dR2
;
    if (SR_No ne 0.0) then begin
	arg = (R-SR_R)^2/(SR_dR^2) + abs(SR_Z)/(SR_dZ)
	Dens_SR = 0. * R
        iloc = where (arg le 86)
	if (iloc(0) ne -1) then Dens_SR(iloc)  = SR_No * exp(-arg(iloc))
    endif else begin
	Dens_SR  = 0.0  * R
    endelse
            
;
; JP mod - make blobs also exp(-z)
;
    if (LB_No ne 0.0) then begin
         LB_Delta = (dTheta-LB_Theta) mod (2*!pi) 
         LB_Delta = LB_Delta + 2*!pi*((LB_Delta lt -!pi)-(LB_Delta gt !pi))
	 Dens_LB  = 0. * R
	 arg = (R-LB_R)^2/(LB_dR^2)         $
               +(LB_Delta)^2/(LB_dTheta^2)   $
               + abs(SR_Z)/(LB_dZ) 
         iloc = where(arg le 86)
         if (iloc(0) ne -1) then Dens_LB(iloc)  = LB_No * exp(-arg(iloc))
    endif else begin
	 Dens_LB  = 0.0 * R
    endelse

    if (TB_No ne 0.0) then begin
         TB_Delta = (dTheta-TB_Theta) mod (2*!pi)
         TB_Delta = TB_Delta + 2*!pi*((TB_Delta lt -!pi)-(TB_Delta gt !pi))
         Dens_TB  = 0. * R
	 arg =(R-TB_R)^2/(TB_dR^2)        $
              +(TB_Delta)^2/(TB_dTheta^2)  $
              +abs(SR_Z)/(TB_dZ) 
         iloc = where(arg le 86)
         if (iloc(0) ne -1) then Dens_TB(iloc)  = TB_No * exp(-arg(iloc))
    endif else begin
	 Dens_TB = 0.0 * R
    endelse
;
; end JP mod

    Dens = Dens_SR + Dens_LB + Dens_TB

    ;
    ; Partials of Density wrt Parameters
    ;
    if (want_partials) then begin

        npts = n_elements(x)
        npar = n_elements(a)
        if (n_elements(df) ne npts*npar) then df = dblarr(npts,npar)

        if (SR_No ne 0.0) then begin
          ; SR_No
          df(*,0) = Dens_SR/SR_No

          ; SR_R
          df(*,1) = Dens_SR * 2.*((R-SR_R))/(SR_dR^2)

          ;SR_dR
          df(*,2) = Dens_SR * 2.*((R-SR_R)^2)/(SR_dR^3)
       
          ; SR_dZ
          df(*,3) = Dens_SR * abs(SR_Z)/(SR_dZ^2)

        endif

        if (LB_No ne 0.0) then begin
          ; LB_No
          df(*,4) = Dens_LB/LB_No

          ; LB_R
          df(*,5) = Dens_LB * 2.*(R-LB_R)/LB_dR^2

          ; LB_dR
          df(*,6) = Dens_LB * 2.*(R-LB_R)^2/LB_dR^3

          ; LB_dTheta
          df(*,8) = Dens_LB * 2.*LB_Delta^2/LB_dTheta^3 * d2r

          ; LB_dZ
          df(*,9) = Dens_LB * abs(SR_Z)/LB_dZ^2

        endif

        if (TB_No ne 0.0) then begin
          ; TB_No
          df(*,10) = Dens_TB/TB_No

          ; TB_R
          df(*,11) = Dens_TB * 2.*(R-TB_R)/TB_dR^2

          ; TB_dR
          df(*,12) = Dens_TB * 2.*(R-TB_R)^2/TB_dR^3

          ; TB_dTheta
          df(*,14) = Dens_TB * 2.*TB_Delta^2/TB_dTheta^3 * d2r

          ; TB_dZ
          df(*,15) = Dens_TB * abs(SR_Z)/TB_dZ^2

        endif

; set up some common factors
;
        signZ = 1.0*(SR_Z ge 0.) - (SR_Z lt 0.)
        dfdz  = -signZ * (Dens_SR/SR_dZ + Dens_LB/LB_dZ + Dens_TB/TB_dZ)

        ; RB_Omega   
        df(*,16) = dfdz*(coso*sini*X + sino*sini*Y) * d2r

        ; RB_Incl    
        df(*,17) = dfdz*(sino*cosi*X - coso*cosi*Y - sini*Z) * d2r


    endif


endif else Dens = x*0

return,Dens
end

;--------------------------------------------------------------
; Function new_isocloud
;
; Returns the number density of an isotropic-like cloud.
; Copied the function for the smooth cloud, and modified.
;
; Added by: Rosalia O'Brien, Fall 2025
;
; Inputs:
;   x,y,z - heliocentric cartesian coordinates
;   R     - sqrt(x^2+y^2+z^2)
;   Theta - Earth mean longitude
;   a     - Model parameters
;
;--------------------------------------------------------------
function new_isocloud,x,y,z,R,a,No,Alpha

  ;  No = 1e-42
  
  FuncIndx = 2

  if (No ne 0.0) then begin

    d2r = !pi/180.

    ;    Alpha    = -50.0
    Beta     = 0.0
    Gamma    = 0.0
    Mu       = a(4)
    Mu2      = a(5)
    Mu3      = a(6)
    Mu4      = a(7)
    Mu5      = a(8)
    Mu6      = a(9)
    Mu7      = a(10)
    Mu8      = a(11)
    Mu9      = a(12)
    Mu10     = a(13)

    Omega    = 0.0 * d2r
    Incl     = 0.0 * d2r
    Xo       = 0.0
    Yo       = 0.0
    Zo       = 0.0

    ;
    ; Some frequently used geometry terms
    ;
    sino = sin(Omega)
    coso = cos(Omega)
    sini = sin(Incl)
    cosi = cos(Incl)

    ;
    ; Translate origin from Sun to cloud center.
    ;
    Xp  = X - Xo
    Yp  = Y - Yo
    Zp  = Z - Zo
    Rc  = sqrt( Xp^2 + Yp^2 + Zp^2 )

    ;
    ; Rotate into coord system of cloud to determine Z-height in cloud
    ; (Zc) and radial distance in cloud symetry plane (Rc).
    ; Perform rotation of Omega around Z-Axis to align X with cloud
    ; line-of-nodes. Then, perform Incl rotation around X-Axis.
    ;
    Zc   = sino*sini*Xp - coso*sini*Yp + cosi*Zp
    Zeta = abs(Zc/Rc)

    ;
    ; Compute Radial power-law parameter
    ; Mark6p2  sets this to Rc, not R    20 May 1996
    ; Mark6p7 allows the polynomial mods only beyond a cutoff radius
    ; Mark6p8 puts is back the way it was, but keeps insistence on positive
    ;         coeffs -- uses R, not Rc
    ;
    AlphaR = Alpha + ( Mu6^2*(R-1.) + Mu7^2*(R-1.)^2 + Mu8^2*(R-1.)^3 )
    ;protection against overflow 12May1996
    stuff = where(AlphaR gt 50.)
    if (stuff(0) ne -1) then AlphaR(stuff) = 50.

    ;
    ; Radial and Latitudinal Density Distribution
    ;
    case FuncIndx of

      ;
      ; John Good Density Function
      ;
      0 : begin
        Dc     = sqrt(Rc^2 - Zc^2)
        ZoD    = abs(Zc/Dc)
        ZoD2G  = ZoD^Gamma
        Dens_Cloud_Vert = exp( -Beta * ZoD2G )
        Dens_Cloud_Rad  = 1./Dc^AlphaR
      end

      ;
      ; Modified Fan Density Function
      ;
      1 : begin
        Zeta2G  = Zeta^Gamma
        Dens_Cloud_Vert = exp( -Beta * Zeta2G )
        Dens_Cloud_Rad  = 1./Rc^AlphaR
      end

      ;
      ; Widened Modified Fan Density Function
      ;
      2 : begin
        GZR    = (0.5*Zeta^2/Mu) * (Zeta lt Mu) + $
          (Zeta-0.5*Mu)   * (Zeta ge Mu)
        GZR2G  = GZR^Gamma
        Dens_Cloud_Vert = exp( -Beta * GZR2G )
        Dens_Cloud_Rad  = 1./Rc^AlphaR
      end

      ;
      ; Ellipsoidal Density Function
      ;
      3 : begin
        Bterm = (1. + (Beta*Zeta)^2)
        Dens_Cloud_Vert = 1./Bterm^Gamma
        Dens_Cloud_Rad  = 1./Rc^AlphaR
      end

      ;
      ; Cosine (Sombrero)
      ;
      4 : begin
        cosB   = sqrt(1-Zeta^2)
        cosB2G = cosB^Gamma
        Dens_Cloud_Vert = (1. + Beta * cosB2G)/(1. + Beta)
        Dens_Cloud_Rad  = 1./Rc^AlphaR
      end

      ;
      ; Modified Fan Density Function with Polynomial Z-height
      ;
      5 : begin
        aZ = abs(Z)
        BetaZ = Beta + Mu3*aZ + Mu4*aZ^2 + Mu5*aZ^3
        Dens_Cloud_Vert = exp( -BetaZ * Zeta )
        Dens_Cloud_Rad  = 1./Rc^AlphaR
      end

    endcase
    ;
    ; Total Density of Cloud Component
    ;
    Dens = No * Dens_Cloud_Vert * Dens_Cloud_Rad

  endif else Dens = x*0

  return,Dens
end


;
;----------------------------------------------------------------------------
; Function GET_PARNAMES
;
; Assign names to ZKERNEL parameters
;
; Written By: BA Franz, ARC, 2/94
; Modified By :  JP March 1996 to add parnames band 1,2,3 phase func
;                JP April 1996 to add second ring radial law
;
;----------------------------------------------------------------------------
function get_parnames

npars        = 256
parname      = strarr(npars)
parname(*)   = 'Unused'

parname(  0) = 'FuncIndx'

parname(  1) = 'PF1_C0'
parname(  2) = 'PF1_C1'
parname(  3) = 'PF1_Ka'
parname(  4) = 'PF2_C0'
parname(  5) = 'PF2_C1'
parname(  6) = 'PF2_Ka'
parname(  7) = 'PF3_C0'
parname(  8) = 'PF3_C1'
parname(  9) = 'PF3_Ka'

parname( 10) = 'To1'
parname( 11) = 'To2'
parname( 12) = 'K'
parname( 13) = 'Delta'

parname( 14) = 'C_No'
parname( 15) = 'C_Alpha'
parname( 16) = 'C_Beta'
parname( 17) = 'C_Gamma'
parname( 18) = 'C_Mu'
parname( 19) = 'C_Mu2'
parname( 20) = 'C_Mu3'
parname( 21) = 'C_Mu4'
parname( 22) = 'C_Mu5'
parname( 23) = 'C_Mu6'
parname( 24) = 'C_Mu7'
parname( 25) = 'C_Mu8'
parname( 26) = 'C_Mu9'
parname( 27) = 'C_Mu10'
parname( 28) = 'C_Omega'
parname( 29) = 'C_Incl'
parname( 30) = 'C_Xo'
parname( 31) = 'C_Yo'
parname( 32) = 'C_Zo'

parname( 33) = 'C_A1'
parname( 34) = 'C_A2'
parname( 35) = 'C_A3'
parname( 36) = 'C_A4'
parname( 37) = 'C_E3'
parname( 38) = 'C_E4'
parname( 39) = 'C_E5'
parname( 40) = 'C_E6'
parname( 41) = 'C_E7'
parname( 42) = 'C_E8'
parname( 43) = 'C_E9'
parname( 44) = 'C_E10'

parname( 45) = 'B1_No'
parname( 46) = 'B1_Dz'
parname( 47) = 'B1_Dr'
parname( 48) = 'B1_Ro'
parname( 49) = 'B1_Vi'
parname( 50) = 'B1_Vr'
parname( 51) = 'B1_Pi'
parname( 52) = 'B1_Pr'
parname( 53) = 'B1_P2r'
parname( 54) = 'B1_Omega'
parname( 55) = 'B1_Incl'
parname( 56) = 'B1_Xo'
parname( 57) = 'B1_Yo'
parname( 58) = 'B1_Zo'

parname( 59) = 'B1_A1'
parname( 60) = 'B1_A2'
parname( 61) = 'B1_A3'
parname( 62) = 'B1_A4'
parname( 63) = 'B1_E3'
parname( 64) = 'B1_E4'
parname( 65) = 'B1_E5'
parname( 66) = 'B1_E6'
parname( 67) = 'B1_E7'
parname( 68) = 'B1_E8'
parname( 69) = 'B1_E9'
parname( 70) = 'B1_E10'

parname( 71) = 'B2_No'
parname( 72) = 'B2_Dz'
parname( 73) = 'B2_Dr'
parname( 74) = 'B2_Ro'
parname( 75) = 'B2_Vi'
parname( 76) = 'B2_Vr'
parname( 77) = 'B2_Pi'
parname( 78) = 'B2_Pr'
parname( 79) = 'B2_P2r'
parname( 80) = 'B2_Omega'
parname( 81) = 'B2_Incl'
parname( 82) = 'B2_Xo'
parname( 83) = 'B2_Yo'
parname( 84) = 'B2_Zo'

parname( 85) = 'B2_A1'
parname( 86) = 'B2_A2'
parname( 87) = 'B2_A3'
parname( 88) = 'B2_A4'
parname( 89) = 'B2_E3'
parname( 90) = 'B2_E4'
parname( 91) = 'B2_E5'
parname( 92) = 'B2_E6'
parname( 93) = 'B2_E7'
parname( 94) = 'B2_E8'
parname( 95) = 'B2_E9'
parname( 96) = 'B2_E10'

parname( 97) = 'B3_No'
parname( 98) = 'B3_Dz'
parname( 99) = 'B3_Dr'
parname(100) = 'B3_Ro'
parname(101) = 'B3_Vi'
parname(102) = 'B3_Vr'
parname(103) = 'B3_Pi'
parname(104) = 'B3_Pr'
parname(105) = 'B3_P2r'
parname(106) = 'B3_Omega'
parname(107) = 'B3_Incl'
parname(108) = 'B3_Xo'
parname(109) = 'B3_Yo'
parname(110) = 'B3_Zo'

parname(111) = 'B3_A1'
parname(112) = 'B3_A2'
parname(113) = 'B3_A3'
parname(114) = 'B3_A4'
parname(115) = 'B3_E3'
parname(116) = 'B3_E4'
parname(117) = 'B3_E5'
parname(118) = 'B3_E6'
parname(119) = 'B3_E7'
parname(120) = 'B3_E8'
parname(121) = 'B3_E9'
parname(122) = 'B3_E10'

parname(123) = 'B4_No'
parname(124) = 'B4_Dz'
parname(125) = 'B4_Dr'
parname(126) = 'B4_Ro'
parname(127) = 'B4_Vi'
parname(128) = 'B4_Vr'
parname(129) = 'B4_Pi'
parname(130) = 'B4_Pr'
parname(131) = 'B4_P2r'
parname(132) = 'B4_Omega'
parname(133) = 'B4_Incl'
parname(134) = 'B4_Xo'
parname(135) = 'B4_Yo'
parname(136) = 'B4_Zo'

parname(137) = 'B4_A1'
parname(138) = 'B4_A2'
parname(139) = 'B4_A3'
parname(140) = 'B4_A4'
parname(141) = 'B4_E3'
parname(142) = 'B4_E4'
parname(143) = 'B4_E5'
parname(144) = 'B4_E6'
parname(145) = 'B4_E7'
parname(146) = 'B4_E8'
parname(147) = 'B4_E9'
parname(148) = 'B4_E10'

parname(149) = 'SR_No'
parname(150) = 'SR_R'
parname(151) = 'SR_dR'
parname(152) = 'SR_dZ'
parname(153) = 'LB_No'
parname(154) = 'LB_R'
parname(155) = 'LB_dR'
parname(156) = 'LB_Theta'
parname(157) = 'LB_dTheta'
parname(158) = 'LB_dZ'
parname(159) = 'TB_No'
parname(160) = 'TB_R'
parname(161) = 'TB_dR'
parname(162) = 'TB_Theta'
parname(163) = 'TB_dTheta'
parname(164) = 'TB_dZ'
parname(165) = 'RB_Omega'
parname(166) = 'RB_Incl'
parname(167) = 'RB_Xo'
parname(168) = 'RB_Yo'
parname(169) = 'RB_Zo'

parname(170) = 'RB_A1'
parname(171) = 'RB_A2'
parname(172) = 'RB_A3'
parname(173) = 'RB_A4'
parname(174) = 'RB_E3'
parname(175) = 'RB_E4'
parname(176) = 'RB_E5'
parname(177) = 'RB_E6'
parname(178) = 'RB_E7'
parname(179) = 'RB_E8'
parname(180) = 'RB_E9'
parname(181) = 'RB_E10'

parname(182) = 'SR_dR2'       ;JP

return,parname
end

;
;-------------------------------------------------------------------------------
; Function GET_PARINDEX
;-------------------------------------------------------------------------------
function get_parindex,parname

    n        = n_elements(parname)
    indx     = intarr(n)
    namelist = strupcase(get_parnames())

    for i=0,n-1 do begin
        name    = strupcase(strtrim(parname(i),2))
        ptr     = where(namelist eq name)
        indx(i) = ptr(0)
    endfor

    if (n eq 1) then indx = indx(0)

    return,indx

end


;
;-------------------------------------------------------------------------------
; Function WIRING2TOK
;-------------------------------------------------------------------------------
function wiring2tok,wiring,input,output

    ptr1 = 0
    ptr2 = strpos(wiring,'>',ptr1)

    if (ptr2 ne -1) then begin

        input = strtrim(strmid(wiring,ptr1,ptr2),2)

        len  = strlen(wiring)
        nout = 0
        done = 0

        while (not done) do begin

            nout = nout+1
            ptr1 = ptr2 + 1
            ptr2 = strpos(wiring,'+',ptr1)
            if (ptr2 eq -1) then begin
                ptr2 = len
                done = 1
            endif

            token = strtrim(strmid(wiring,ptr1,ptr2-ptr1),2)

            case nout of
              1    : output = token
              else : output = [output,token]
            endcase

        endwhile

    endif

    return, nout

end


;
;-------------------------------------------------------------------------------
; Function ZPAR_WIRING
;
; Builds the parameter wiring index from a list of wiring strings.  If 
; list is provided, it just returns the last installed index array.  If
; no wiring was installed previously, and no wiring diagram was supplied,
; the function just returns a clean index array.
;
; BAF, 11/95
;-------------------------------------------------------------------------------
function zpar_wiring,wiring,reset=reset

    common zpar_wiring, wiring_indx

    npar = 256

    ;
    ; If provided, convert wiring diagram to wiring index array
    ;
    if ( ((nset = n_elements(wiring))) gt 0 ) then begin

	;
	; Start with a cleared index
        ;
        wiring_indx = indgen(npar)

        ;
        ; Loop through each string, parsing to extract the input parameter
        ; index and the list of output parameter indices.
        ;
        for iset = 0,nset-1 do begin
            nout = wiring2tok(wiring(iset),input,output)
            if (nout ne 0) then begin

                ; Get indices
                indx_in  = get_parindex(input)
                indx_out = get_parindex(output)

                ; Validate
                if (min([indx_in,indx_out]) lt 0) then begin
                    message,'Invalid parameter name',/info
                    stop
                endif

                ; Update the wiring index array.  We must account for
                ; dependencies from previous strings, so we have to
                ; search the wiring index for each output index.
                for i=0,nout-1 do begin
                    s = where( wiring_indx eq indx_out(i) )
                    wiring_indx(s) = indx_in
                 endfor

            endif
        endfor

    endif else if (n_elements(wiring_indx) eq 0) then begin
         ;
         ; If it doesn't exist, create a dummy
         ;
         wiring_indx = indgen(npar)
    endif else if (keyword_set(reset)) then begin
         ;
         ; Clear the wiring index
         ;
         wiring_indx = indgen(npar)
    endif

    return,wiring_indx

end



;;
;;-------------------------------------------------------------------------------
;; Function PARFIX
;;
;; Function to fix parameters of ZKERNEL based on wavelength.
;;
;; Written By: BA Franz, ARC, 4/94
;; Modified By: JP March 1996 for wavelength-dependent phase function
;; Modified By: JP May 26 1996 for migband changes
;;
;;-------------------------------------------------------------------------------
;function parfix,wl,a,indxpar
;
;npar  = n_elements(a)
;indxa = indgen(npar)
;if (n_elements(indxpar) eq 0) then indxpar = indxa
;
;windx = zpar_wiring()
;
;;
;; Names of parameters which are currently inoperable or not optimizable
;;
;pfix = 'FuncIndx'
;pfix = [pfix,'B1_Ro','B1_Vr','B1_P2r']
;pfix = [pfix,'B2_Ro','B2_Vr','B2_P2r']
;pfix = [pfix,'B3_Ro','B3_Vr','B3_P2r']
;pfix = [pfix,'B4_Ro','B4_Vr','B4_P2r']
;pfix = [pfix,'LB_Theta']
;pfix = [pfix,'TB_Theta']
;pfix = [pfix,'RB_Xo','RB_Yo','RB_Zo']
;
;;
;; Names of parameters which require a specific wavelength
;;
;s = where(wl eq 1.25)  &  if (s(0) eq -1) then begin
;    pfix = [pfix,'C_A1','B1_A1','B2_A1','B3_A1','B4_A1','RB_A1']
;    pfix = [pfix,'PF1_C0','PF1_C1','PF1_Ka']
;endif
;s = where(wl eq 2.20)  &  if (s(0) eq -1) then begin
;    pfix = [pfix,'C_A2','B1_A2','B2_A2','B3_A2','B4_A2','RB_A2']
;    pfix = [pfix,'PF2_C0','PF2_C1','PF2_Ka']
;endif
;s = where(wl eq 3.50)  &  if (s(0) eq -1) then begin
;    pfix = [pfix,'C_A3','B1_A3','B2_A3','B3_A3','B4_A3','RB_A3']
;    pfix = [pfix,'C_E3','B1_E3','B2_E3','B3_E3','B4_E3','RB_E3']
;    pfix = [pfix,'PF3_C0','PF3_C1','PF3_Ka']
;endif
;s = where(wl eq 4.90)  &  if (s(0) eq -1) then begin
;    pfix = [pfix,'C_A4','B1_A4','B2_A4','B3_A4','B4_A4','RB_A4']
;    pfix = [pfix,'C_E4','B1_E4','B2_E4','B3_E4','B4_E4','RB_E4']
;endif
;s = where(wl eq 12.)  &  if (s(0) eq -1) then begin
;    pfix = [pfix,'C_E5','B1_E5','B2_E5','B3_E5','B4_E5','RB_E5']
;endif
;s = where(wl eq 25.)  &  if (s(0) eq -1) then begin
;    pfix = [pfix,'C_E6','B1_E6','B2_E6','B3_E6','B4_E6','RB_E6']
;endif
;s = where(wl eq 60.)  &  if (s(0) eq -1) then begin
;    pfix = [pfix,'C_E7','B1_E7','B2_E7','B3_E7','B4_E7','RB_E7']
;endif
;s = where(wl eq 100.)  &  if (s(0) eq -1) then begin
;    pfix = [pfix,'C_E8','B1_E8','B2_E8','B3_E8','B4_E8','RB_E8']
;endif
;s = where(wl eq 140.)  &  if (s(0) eq -1) then begin
;    pfix = [pfix,'C_E9','B1_E9','B2_E9','B3_E9','B4_E9','RB_E9']
;endif
;s = where(wl eq 240.)  &  if (s(0) eq -1) then begin
;    pfix = [pfix,'C_E10','B1_E10','B2_E10','B3_E10','B4_E10','RB_E10']
;endif
;
;;
;; Convert names to indices of parameters we want fixed
;;
;pfixIndx = get_parindex(pfix)
;
;;
;; Fix params which are unused
;;
;s = where(a eq -1)
;if (s(0) ne -1) then pfixIndx = [pfixIndx,s]
;
;;
;; Fix params for components which are turned-off
;;
;s = get_parindex(['K','To2'])
;if (a(s(0)) eq 0.0) then $
;    pfixIndx = [pfixIndx,s]
;
;s = get_parindex('B1_No')
;if (a(s(0)) eq 0.0) then $
;    pfixIndx = [pfixIndx,s(0)+indgen(26)]
;s = get_parindex('B2_No')
;if (a(s(0)) eq 0.0) then $
;    pfixIndx = [pfixIndx,s(0)+indgen(26)]
;s = get_parindex('B3_No')
;if (a(s(0)) eq 0.0) then $
;    pfixIndx = [pfixIndx,s(0)+indgen(26)]
;s = get_parindex('B4_No')
;if (a(s(0)) eq 0.0) then $
;    pfixIndx = [pfixIndx,s(0)+indgen(26)]
;
;s = get_parindex('SR_No')
;if (a(s(0)) eq 0.0) then $
;    pfixIndx = [pfixIndx,s(0)+indgen(4)]
;s = get_parindex('LB_No')
;if (a(s(0)) eq 0.0) then $
;    pfixIndx = [pfixIndx,s(0)+indgen(6)]
;s = get_parindex('TB_No')
;if (a(s(0)) eq 0.0) then $
;    pfixIndx = [pfixIndx,s(0)+indgen(6)]
;s = get_parindex(['SR_No','LB_No','TB_No'])
;if (max(a(s)) eq 0.0) then $
;    pfixIndx = [pfixIndx,s(0)+indgen(33)]
;
;;
;; Enforce unique, sorted order
;;
;if (n_elements(pfixIndx) gt 1) then begin
;    pfixIndx = pfixIndx(uniq(pfixIndx,sort(pfixIndx)))
;endif
;if (n_elements(indxpar) gt 1) then begin
;    indxpar = indxpar(uniq(indxpar,sort(indxpar)))
;endif
;
;;
;; Apply pfixIndx to indxpar
;;
;s = nowhere(indxpar,pfixIndx)
;if (s(0) ne -1) then indxpar = indxpar(s)
;
;;
;; Fix parameters which are wired to other parameters
;;
;s = anywhere(indxpar,windx)
;if (s(0) ne -1) then		$
;    indxpar = indxpar(s)	$
;else				$
;    print,'Warning: parfix found no free parameters'
;
;return,indxpar
;end
;
;
;
;-------------------------------------------------------------------------------
; IDL Function ZKERNEL.PRO
;
; Zodiacal Dust Model
;
; Primary Authors: B.A. Franz, ARC, (301) 513-7776
;                  W.T. Reach, USRA (301) 286-4255
; Modified By:     JP March 1996 for phase function/scattering function
;                  JP April 1996 for second radial law
;                  JP April 1996 as per HTF 4/16/96 code to include max. dist
;                                from sun for integration limits
;                                HTF CODE IN CAPS.
;
; Inputs:
;    x - n-element structured array of independent variables which at least
;        contains the following fields.  (See ZODI_STR.PRO)        
;        x.longitude  : Heliocentric Ecliptic longitude of LOS in degrees
;        x.latitude   : Heliocentric Ecliptic latitude of LOS in degrees
;        x.day1990    : 1990 Day Number of Observation
;        x.wave_len   : Wavelength in um
;
;    a - m-element floating array of model parameters
;
; Optional Inputs:
;    indxpar=indxpar - array of indices of variable parameters (Def=0..m-1)
;    dbwave=alternate array of wavelengths
;
; Outputs:
;    f  - n-element floating array of Model Flux in MJy/Sr
;
; Optional Outputs:
;    df - nxm array of partial derivatives for model parameters at each x
;
;-------------------------------------------------------------------------------
;
pro zkernel,data,a,phase_type,f,df,indxpar=indxpar,losinfo=losinfo,no_colcorr=no_colcorr,dbwave=dbwave,solar_irr=solar_irr,new_iso_comp=new_iso_comp,iso_comp_only=iso_comp_only

if (keyword_set(losinfo)) then want_los_info=1 else want_los_info=0

;
; Define some useful params
;
RMAX = 5.2                      ; JUPITER'S ORBIT RADIUS
NUMBER_OF_STEPS = 50            ; THIS MANY FUNCTION EVALUATIONS PER L.O.S.

npts = n_elements(data)  				; 
d2r  = !pi/180.D         				; degrees to radians
eps  = double(1e-20)     				; > zero
if (n_elements(dbwave) ne 10) then dbwave = [1.25,2.2,3.5,4.9,12.,25.,60.,100.,140.,240.]
nmsg = 50000

;
; Set wiring
;
windx = zpar_wiring()
a = a(windx)

;
; Set density function selection flag
;
FuncIndx = fix(a(0))

;
; Get Scatt Function Parameters (JP)
;
nScatt = 9
iScatt = 1 + indgen(nScatt)
aScatt = a(iScatt)
if (a[183] ne -1) then begin 
  aScatt = a[183:*] 
  aScatt = aScatt[where(aScatt ne -1)]; location of [g1,g2,g3,w1,w2,w3] for 3 HG phase function
endif

;
; Get Therm Function Parameters
;
nTherm = 4
iTherm = 10 + indgen(nTherm)
aTherm = a(iTherm)

;
; Get Cloud Parameters
;
nDens_C = 19				; Number of density params
iDens_C = 14 + indgen(nDens_C) 		; Index of density params in a
aDens_C = a(iDens_C)			; Vector of density params
nSrc_C  = 12				; Number of emissivity params
iSrc_C  = 33 + indgen(nSrc_C) 		; Index of emissivity params in a
aSrc_C  = a(iSrc_C)			; Vector of emissivity params

;
; Get Dust Band 1 Parameters
;
nDens_B1 = 14				; Number of density params
iDens_B1 = 45 + indgen(nDens_B1)	; Index of density params in a
aDens_B1 = a(iDens_B1)			; Vector of density params
nSrc_B1  = 12				; Number of emissivity params
iSrc_B1  = 59 + indgen(nSrc_B1) 	; Index of emissivity params in a
aSrc_B1  = a(iSrc_B1)			; Vector of emissivity params

;
; Get Dust Band 2 Parameters
;
nDens_B2 = 14				; Number of density params
iDens_B2 = 71 + indgen(nDens_B2)	; Index of density params in a
aDens_B2 = a(iDens_B2)			; Vector of density params
nSrc_B2  = 12				; Number of emissivity params
iSrc_B2  = 85 + indgen(nSrc_B2) 	; Index of emissivity params in a
aSrc_B2  = a(iSrc_B2)			; Vector of emissivity params

;
; Get Dust Band 3 Parameters
;
nDens_B3 = 14				; Number of density params
iDens_B3 = 97 + indgen(nDens_B3)	; Index of density params in a
aDens_B3 = a(iDens_B3)			; Vector of density params
nSrc_B3  = 12				; Number of emissivity params
iSrc_B3  = 111 + indgen(nSrc_B3) 	; Index of emissivity params in a
aSrc_B3  = a(iSrc_B3)			; Vector of emissivity params

;
; Get Dust Band 4 Parameters
;
nDens_B4 = 14				; Number of density params
iDens_B4 = 123 + indgen(nDens_B4)	; Index of density params in a
aDens_B4 = a(iDens_B4)			; Vector of density params
nSrc_B4  = 12				; Number of emissivity params
iSrc_B4  = 137 + indgen(nSrc_B4) 	; Index of emissivity params in a
aSrc_B4  = a(iSrc_B4)			; Vector of emissivity params

;
; Get Solar Ring Parameters
;
nDens_RB = 21				; Number of density params
iDens_RB = 149 + indgen(nDens_RB)	; Index of density params in a
nDens_RB = nDens_RB + 1                 ; JP add in SR_dR2
iDens_RB = [iDens_RB,182]               ; JP add in SR_dR2
aDens_RB = a(iDens_RB)			; Vector of density params
nSrc_RB  = 12				; Number of emissivity params
iSrc_RB  = 170 + indgen(nSrc_RB) 	; Index of emissivity params in a
aSrc_RB  = a(iSrc_RB)			; Vector of emissivity params


;
; Associate wavelengths with fullband detector numbers
; 1a=1b=1c=0 .... 9=b10
;
if phase_type eq 'kelsall' then begin
detnum = bytarr(npts) - 1
for i=0,9 do $
    detnum = temporary(detnum) + ( data.wave_len eq dbwave(i) ) * (i+1)
endif

;
; Calculate position of Earth in heliocentric ecliptic coords and
; Solar elongation of LOS
;
earthsun,data.day1990,data.longitude,data.latitude, $
         SolElong,Earth_Dis,Earth_Lon,Earth_Mean_Lon

; THIS PART GETS COMMENTED OUT WHEN CHANGE TO HELIOCENTRIC INTEGRATION RANGES
;------------------------------------------------
; Set grid, weights for Gauss-Leguerre quadrature
; Define LOS Integration Points (from Earth in AU)
;
;gaussint,0.,5.,los,gqwts,numpts=50
;------------------------------------------------
;
; Set-up for calculation of integrated flux
;
if (n_elements(f) ne npts) then f = dblarr(npts)

;
; Set-up for calculation of partial derivatives
;
if (n_params(0) gt 3) then begin

    want_partials = 1

    ;
    ; Create index of variable parameters, if not supplied.
    ;
    npar = n_elements(a)
    if (n_elements(indxpar) eq 0) then begin
        nterms  = npar
        indxpar = indgen(nterms)
    endif else $
        nterms = n_elements(indxpar)

    ;
    ; Create storage for partials, and create mask to expedite selection 
    ; of parameters which require partials.
    ;
    if (n_elements(df) ne npts*nterms) then $
        df = dblarr(npts,nterms)            $
    else                                    $
        df = df * 0.0

    parnum  = intarr(npar)-1  &  parnum(indxpar) = indgen(nterms)

endif else want_partials = 0

;
; Set-up to store line-of-sight info
;
if (want_los_info) then begin
    element_rec = { elemnt_str, $
                    earth_dist : 0.0, $
                    solar_dist : 0.0, $
                    ecliptic_z : 0.0, $
                    density    : 0.0, $
                    flux       : 0.0, $
                    int_wts    : 0.0  $
                  }

    losdata = replicate( $
              { wave_len   : 0.0, $
                pixel_no   : 0L,  $
                lon        : 0.0, $
                lat        : 0.0, $
                elong      : 0.0, $
                day        : 0.0, $
                earth_lon  : 0.0, $
                earth_rad  : 0.0, $
                ;element    : replicate( {elemnt_str},n_elements(los) ), $
                zodi       : 0.0 }, npts)
endif

;
; Loop over each LOS, Time, Wavelength
;
for ilos = 0L,npts-1 do begin

    ; 
    ; Compute Basic Geometry
    ; ----------------------
    lat  = double(data(ilos).latitude   * d2r )
    lon  = double(data(ilos).longitude  * d2r )
    ;
    ; Position of Earth
    ;
    Re    = Earth_Dis(ilos)
    Theta = Earth_Lon(ilos)

    COSTHETA = COS(THETA)
    SINTHETA = SIN(THETA)
    COSLON   = COS(LON)
    SINLON   = SIN(LON)
    COSLAT   = COS(LAT)
    SINLAT   = SIN(LAT)

    X0 = RE*COSTHETA  & Y0 = RE*SINTHETA        ; &Z0 = 0.0
    B  = 2.*( X0*COSLAT*COSLON + Y0*COSLAT*SINLON )
    C  = RE*RE - RMAX*RMAX
    Q  = -0.5*B*(1.0 + SQRT( B*B - 4.*C)/ABS(B))
    RANGE = MAX( [Q, C/Q] )

;  Set up Integration abscissae and wts, depending upon ecliptic latitude
;  Change by JP 08 May 1996  set to 20 rather than 15
;
betathresh = 20. * d2r
;
; THIS PART GETS INSIDE THE LOOP WHEN CHANGE TO HELIOCENTRIC INTEGRATION RANGES
;------------------------------------------------
; Set grid, weights for Gauss-Leguerre quadrature
; Define LOS Integration Points (from Earth in AU)
;
;if (abs(lat) gt betathresh) then begin
;    gaussint,0.,RANGE,los,gqwts,numpts=NUMBER_OF_STEPS
;endif else begin
simpint,0., RANGE,los,gqwts,stepsize=0.025
;endelse
;------------------------------------------------
    ;
    ; Position of LOS elements in Heliocentric Ecliptic coord system.
    ;
    Sxy   = los*COSLAT
    X     = Re*COSTHETA + Sxy*COSLON
    Y     = Re*SINTHETA + Sxy*SINLON
    Z     = los*SINLAT
    R     = sqrt( X^2 + Y^2 + Z^2 )

    ;
    ; Source Function per LOS element
    ; ------------------------------------------------------------------
    Lambda = double(data(ilos).wave_len)
    if phase_type eq 'kelsall' then begin
      Det    = detnum(ilos)
    endif
    if phase_type eq 'skysurf' then begin
      Det    = 0
    endif
    Scatt  = scattfunc(phase_type,Det,Lambda,los,R,Re,SolElong(ilos),aScatt, solar_irr=solar_irr);,dScatt)
    Therm  = thermfunc(Det,Lambda,R,aTherm,dTherm,no_colcorr=no_colcorr)

    ; 
    ; Calculate Density and Source Associated With Smooth Cloud Component
    ; -------------------------------------------------------------------
    Dens_C = zcloud(x,y,z,R,aDens_C,dDens_C,func=FuncIndx)
    Src_C  = zsrcfunc(Det,Scatt,Therm,aSrc_C,dSrc_C,dSrc_dScatt_C,dSrc_dTherm_C,phase_type)

    ;
    ; Calculate Density and Source Associated with Dust Bands
    ; -------------------------------------------------------------------
    Dens_B1 = migband(x,y,z,R,aDens_B1,dDens_B1)
    Src_B1  = zsrcfunc(Det,Scatt,Therm,aSrc_B1,dSrc_B1,dSrc_dScatt_B1,dSrc_dTherm_B1,phase_type)
    Dens_B2 = migband(x,y,z,R,aDens_B2,dDens_B2)
    Src_B2  = zsrcfunc(Det,Scatt,Therm,aSrc_B2,dSrc_B2,dSrc_dScatt_B2,dSrc_dTherm_B2,phase_type)
    Dens_B3 = migband(x,y,z,R,aDens_B3,dDens_B3)
    Src_B3  = zsrcfunc(Det,Scatt,Therm,aSrc_B3,dSrc_B3,dSrc_dScatt_B3,dSrc_dTherm_B3,phase_type)
    Dens_B4 = migband(x,y,z,R,aDens_B4,dDens_B4)
    Src_B4  = zsrcfunc(Det,Scatt,Therm,aSrc_B4,dSrc_B4,dSrc_dScatt_B4,dSrc_dTherm_B4,phase_type)

    ;
    ; Calculate Density and Source Associated with Solar Ring/Blobs
    ; -------------------------------------------------------------------
    Dens_RB = solring(x,y,z,R,Earth_Mean_Lon(ilos),aDens_RB,dDens_RB)
    Src_RB  = zsrcfunc(Det,Scatt,Therm,aSrc_RB,dSrc_RB,dSrc_dScatt_RB,dSrc_dTherm_RB,phase_type)
    
    ; Isotropic-like cloud component (optional)
    Dens_new = Dens_C*0
    if (keyword_set(new_iso_comp)) then begin
      No = 5.62e-10
      Alpha = -2.02
      Dens_new = new_isocloud(x,y,z,R,aDens_C,No,Alpha)
    endif
    
    if (keyword_set(iso_comp_only)) then begin
      Dens_C = Dens_C*0
      Dens_B1 = Dens_B1*0
      Dens_B2 = Dens_B2*0
      Dens_B3 = Dens_B3*0
      Dens_B4 = Dens_B4*0
      Dens_RB = Dens_RB*0
;      print,Dens_C
    endif

    ;
    ; Calculate Total Zodi Brightness per LOS Element
    ; --------------------------------------------------------------------
    Flux = Src_C  * Dens_C  + $
           Src_B1 * Dens_B1 + $
           Src_B2 * Dens_B2 + $
           Src_B3 * Dens_B3 + $
           Src_B4 * Dens_B4 + $
           Src_RB * Dens_RB + $
           Src_C  * Dens_new
   
    ;
    ; Integrate Along LOS to Get Total Intensity
    ; ------------------------------------------
    f(ilos) = total(gqwts*Flux)
    
    ;
    ; Store line-of-sight info (Needed by ZODI Explorer)
    ; --------------------------------------------------
    if (want_los_info) then begin
        Dens = Dens_C + Dens_B1 + Dens_B2 + Dens_B3 + Dens_RB
        losdata(ilos).wave_len             = Lambda
        losdata(ilos).pixel_no             = data(ilos).pixel_no
        losdata(ilos).lon                  = Lon / d2r
        losdata(ilos).lat                  = Lat / d2r
        losdata(ilos).elong                = SolElong(ilos) / d2r
        losdata(ilos).day                  = data(ilos).day1990
        losdata(ilos).earth_lon            = Theta / d2r
        losdata(ilos).earth_rad            = Re
;        losdata(ilos).element.earth_dist   = LOS
;        losdata(ilos).element.solar_dist   = R
;        losdata(ilos).element.ecliptic_z   = Z
;        losdata(ilos).element.density      = Dens
;        losdata(ilos).element.flux         = Flux
;        losdata(ilos).element.int_wts      = gqwts
        losdata(ilos).zodi = f(ilos)
    endif

    ;
    ; Calculate Derivatives of Varying Parameters if Requested
    ; --------------------------------------------------------
;    if (want_partials) then begin
;
;        ;
;        ; Scatt Function Partials
;        ;
;        for ii = 0,nScatt-1 do begin
;            ipar = parnum(windx(iScatt(ii)))
;            if (ipar ge 0) then begin
;                dSrcS = dScatt(*,ii) * ( $
;                        dSrc_dScatt_C  * Dens_C  + $
;                        dSrc_dScatt_B1 * Dens_B1 + $
;                        dSrc_dScatt_B2 * Dens_B2 + $
;                        dSrc_dScatt_B3 * Dens_B3 + $
;                        dSrc_dScatt_B4 * Dens_B4 + $
;                        dSrc_dScatt_RB * Dens_RB )
;                df(ilos,ipar) = df(ilos,ipar) + total(gqwts*dSrcS)
;            endif
;        endfor
;        ;
;        ; Therm Function Partials
;        ;
;        for ii = 0,nTherm-1 do begin
;            ipar = parnum(windx(iTherm(ii)))
;            if (ipar ge 0) then begin
;                dSrcT = dTherm(*,ii) * ( $
;                        dSrc_dTherm_C  * Dens_C  + $
;                        dSrc_dTherm_B1 * Dens_B1 + $
;                        dSrc_dTherm_B2 * Dens_B2 + $
;                        dSrc_dTherm_B3 * Dens_B3 + $
;                        dSrc_dTherm_B4 * Dens_B4 + $
;                        dSrc_dTherm_RB * Dens_RB )
;                df(ilos,ipar) = df(ilos,ipar) + total(gqwts*dSrcT)
;            endif
;        endfor
;
;        ;
;        ; Cloud Partials
;        ;
;        if (aDens_C(0) ne 0.0) then begin
;            for ii = 0,nDens_C-1 do begin
;                ipar = parnum(windx(iDens_C(ii)))
;                if (ipar ge 0) then $
;                    df(ilos,ipar) = df(ilos,ipar) $
;                                  + total(gqwts*Src_C*dDens_C(*,ii))
;            endfor
;            for ii = 0,nSrc_C-1 do begin
;                ipar = parnum(windx(iSrc_C(ii)))
;                if (ipar ge 0) then $
;                    df(ilos,ipar) = df(ilos,ipar) $
;                                  + total(gqwts*dSrc_C(*,ii)*Dens_C)
;            endfor
;        endif
;
;        ;
;        ; Band 1 Partials
;        ;
;        if (aDens_B1(0) ne 0.0) then begin
;            for ii = 0,nDens_B1-1 do begin
;                ipar = parnum(windx(iDens_B1(ii)))
;                if (ipar ge 0) then $
;                    df(ilos,ipar) = df(ilos,ipar) $
;                                  + total(gqwts*Src_B1*dDens_B1(*,ii))
;            endfor
;            for ii = 0,nSrc_B1-1 do begin
;                ipar = parnum(windx(iSrc_B1(ii)))
;                if (ipar ge 0) then $
;                    df(ilos,ipar) = df(ilos,ipar) $
;                                  + total(gqwts*dSrc_B1(*,ii)*Dens_B1)
;            endfor
;        endif
;
;        ;
;        ; Band 2 Partials
;        ;
;        if (aDens_B2(0) ne 0.0) then begin
;            for ii = 0,nDens_B2-1 do begin
;                ipar = parnum(windx(iDens_B2(ii)))
;                if (ipar ge 0) then $
;                    df(ilos,ipar) = df(ilos,ipar) $
;                                  + total(gqwts*Src_B2*dDens_B2(*,ii))
;            endfor
;            for ii = 0,nSrc_B2-1 do begin
;                ipar = parnum(windx(iSrc_B2(ii)))
;                if (ipar ge 0) then $
;                    df(ilos,ipar) = df(ilos,ipar) $
;                                  + total(gqwts*dSrc_B2(*,ii)*Dens_B2)
;            endfor
;        endif
;
;        ;
;        ; Band 3 Partials
;        ;
;        if (aDens_B3(0) ne 0.0) then begin
;            for ii = 0,nDens_B3-1 do begin
;                ipar = parnum(windx(iDens_B3(ii)))
;                if (ipar ge 0) then $
;                    df(ilos,ipar) = df(ilos,ipar) $
;                                  + total(gqwts*Src_B3*dDens_B3(*,ii))
;            endfor
;            for ii = 0,nSrc_B3-1 do begin
;                ipar = parnum(windx(iSrc_B3(ii)))
;                if (ipar ge 0) then $
;                    df(ilos,ipar) = df(ilos,ipar) $
;                                  + total(gqwts*dSrc_B3(*,ii)*Dens_B3)
;            endfor
;        endif
;
;        ;
;        ; Band 4 Partials
;        ;
;        if (aDens_B4(0) ne 0.0) then begin
;            for ii = 0,nDens_B4-1 do begin
;                ipar = parnum(windx(iDens_B4(ii)))
;                if (ipar ge 0) then $
;                    df(ilos,ipar) = df(ilos,ipar) $
;                                  + total(gqwts*Src_B4*dDens_B4(*,ii))
;            endfor
;            for ii = 0,nSrc_B4-1 do begin
;                ipar = parnum(windx(iSrc_B4(ii)))
;                if (ipar ge 0) then $
;                    df(ilos,ipar) = df(ilos,ipar) $
;                                  + total(gqwts*dSrc_B4(*,ii)*Dens_B4)
;            endfor
;        endif
;
;        ;
;        ; Solar Ring/Blob Partials
;        ;
;        if (aDens_RB(0) ne 0.0) then begin
;            for ii = 0,nDens_RB-1 do begin
;                ipar = parnum(windx(iDens_RB(ii)))
;                if (ipar ge 0) then $
;                    df(ilos,ipar) = df(ilos,ipar) $
;                                  + total(gqwts*Src_RB*dDens_RB(*,ii))
;            endfor
;            for ii = 0,nSrc_RB-1 do begin
;                ipar = parnum(windx(iSrc_RB(ii)))
;                if (ipar ge 0) then $
;                    df(ilos,ipar) = df(ilos,ipar) $
;                                  + total(gqwts*dSrc_RB(*,ii)*Dens_RB)
;            endfor
;        endif
;
;    endif

    if ((ilos+1) mod nmsg eq 0) then print,'LOS # ',ilos

endfor ; Loop over LOS

;
; Replace returned flux with LOS structure if requested
;
if (want_los_info) then f = temporary(losdata)

return
end

