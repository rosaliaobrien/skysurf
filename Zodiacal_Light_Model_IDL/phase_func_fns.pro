; This procedure file defines some functions necessary for the model to work.
; Author: Rosalia O'Brien
; 

; Define the Hong phase function
function hong_phase_func,x,a
  nf = n_elements(a)/2
  phi = fltarr(n_elements(x))
  for i=0,nf-1 do begin 
    g = a[i]
    w = a[i+nf]
    phi += w*(1-g^2)/(1+g^2-2*g*cos(x))^1.5
  endfor
return,phi
end

; Define a function that numerically integrates over the phase function
; Note: Phase function must integrate to 1 over 4pi steradians. So the integral
; is done in spherical coordinates.
;
; Parameters:
; -----------
; hg_params - arr
;   The Henyey-Greenstein parameters used in hong_phase_func
; phase_norm_fn - float
;   Number that is multiplied by the phase function to ensure it normalizes to 1.
;   This parameter is only set when TESTING. It should be ignored during normal runs of the model.
;
function get_phase_integral,hg_params,phase_norm_fn=phase_norm_fn

  ; Define the limits of integration
  theta_lower = 0.0
  theta_upper = !PI
  phi_lower = 0.0
  phi_upper = 2.0 * !PI

  ; Define the number of points for tabulation
  n_theta = 1000
  n_phi = 1000

  ; Create a grid of points
  theta = findgen(n_theta) / (n_theta - 1) * (theta_upper - theta_lower) + theta_lower
  phi = findgen(n_phi) / (n_phi - 1) * (phi_upper - phi_lower) + phi_lower

  phase_theta = phase_norm_fn*hong_phase_func(theta,hg_params)*sin(theta)
  phase_phi = phase_norm_fn*hong_phase_func(phi,hg_params)

  ; Perform the integration using INT_TABULATED
  theta_integration = fltarr(n_phi)+int_tabulated(theta, phase_theta)

  total_integral = int_tabulated(phi, theta_integration)

  ;PRINT, 'The integral is: ', total_integral

  return,total_integral

end

; Function to get albedo as a function of wavelength (lam)
; See O'Brien+2025
function get_albedo,lam

  albedo_1p6 = 0.11298*1.6+0.08231

  ; Define wavelength and corresponding albedo arrays
  lam_list    = [1.6, 2.2, 3.5, 4.9]
  albedo_list = [albedo_1p6, 0.255, 0.210, 0]

  if lam gt 1.6 then begin
    
    albedo = interpol(albedo_list,lam_list,lam)
    
  endif else begin 
    
    albedo = 0.11298*lam+0.08231
    
  endelse

  return, albedo
end

; Function to get multiplier as a function of wavelength (lam)
; Should be = 1 for nominal HST wavelengths. When trying to extrapolate out to 3.5 micron,
; a multiplier is added to match Kelsall better.
; I can't adjust the albedo directly, because the albedo is tied to the emissivity.
; See O'Brien+2025
function get_mult,lam

  ; Define wavelength and corresponding multiplier arrays
  lam_list    = [1.6, 2.2, 3.5, 4.9]
  mult_list   = [1.0, 1.227, 1.259, 1.0]

  if lam gt 1.6 then begin

    mult = interpol(mult_list,lam_list,lam)

  endif else begin

    mult = 1.0

  endelse

  return, mult
end

; When trying to extrapolate out to 3.5 micron,
; you need to adjust the emissivity too.
; See O'Brien+2025
function get_emiss,lam

  ; Define wavelength and corresponding multiplier arrays
  lam_list    = [1.6, 2.2, 3.5, 4.9]
  emm_list   = [0, 0, 1.66, 0.997]

  if lam gt 1.6 then begin

    emm = interpol(emm_list,lam_list,lam)

  endif else begin

    emm = 1.0

  endelse

  return, emm
end


; Function to get phase function parameters as a function of wavelength (lam)
; See O'Brien+2025
function get_hong_params,lam

  ; If lambda is less than 0.25 or greater than 1.8, then utilize phase function parameters for 0.25/ 1.8
  ; These phase functions have not been tested outside of this wavelength range
  lam_use = (lam gt 1.6) ? 1.6 : ((lam lt 0.25) ? 0.25 : lam)
  
  g1 = 0.24958*lam_use+0.11571
  g2 = 0.05428*lam_use-0.30864
  g3 = -0.87036
  w1 = 0.00183*lam_use+0.04775
  w2 = -0.00143*lam_use+0.03122
  w3 = 0.00030
  
  hg_arr = [g1, g2, g3, w1, w2, w3]
  
;  inte = get_phase_integral(hg_arr, phase_norm_fn=1)
;  print,'Integral =', inte
  
  return, hg_arr
  
end