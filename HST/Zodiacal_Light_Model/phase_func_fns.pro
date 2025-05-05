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

function get_albedo,lam

  ; Define wavelength and corresponding albedo arrays
  lam_list    = [2.2, 3.5]
  albedo_list = [0.255, 0.21]

  if lam gt 1.8 then begin
    ; Find the index of the closest wavelength in lam_list
    idx = (min(abs(lam_list - lam), minval))
    albedo = albedo_list[idx]
  endif else begin
    albedo = 0.113*lam+0.082
  endelse

  return, albedo
end

function get_hong_params,lam

  ; If lambda is less than 0.25 or greater than 1.8, then utilize phase function parameters for 0.25/ 1.8
  ; These phase functions have not been tested outside of this wavelength range
  lam_use = (lam gt 1.8) ? 1.8 : ((lam lt 0.25) ? 0.25 : lam)
  
  g1 = 0.2496*lam_use+0.1157
  g2 = 0.0543*lam_use-0.3086
  g3 = -0.8704
  w1 = 0.0018*lam_use+0.0478
  w2 = -0.0014*lam_use+0.031
  w3 = 0.0003
  
  hg_arr = [g1, g2, g3, w1, w2, w3]
  
  inte = get_phase_integral(hg_arr, phase_norm_fn=1)
  print,'Integral =', inte
  
  return, hg_arr
  
end