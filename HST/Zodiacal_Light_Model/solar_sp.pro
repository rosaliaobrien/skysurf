; Function to get the solar spectrum estimate for any HST filter.
; The input parameter is "wl" or the filter pivot wavelength. This must match
; the wavelength defined in get_zmod.pro
; 
; The solar spectrum values were estimated using Meftah+2018 solar irradiance spectrum and finding the 
; irradiance that HST would measure for each filter, using synphot in Python.
; Meftah+2018: https://www.aanda.org/articles/aa/full_html/2018/03/aa31316-17/aa31316-17.html
; 
; Author: Rosalia O'Brien
; 
function solar_sp, wl

  hst_flux   = [111835192.98326491, 143236136.2232679, 181404234.08055428, 200760536.90616906, 217740136.09056026, 237769811.5880849, 238564767.8989514, 240585677.00454003, 1324505.632015205, 5749844.193788173, 10869822.237373352, 33408097.604545005, 70728631.80980875, 114534021.98386078, 145482363.320424, 148777018.56162617, 177230311.92482492, 199488651.33783838, 215400106.40923926, 237706012.85011584, 238548500.14506322, 241000606.87217727, 241905098.72034737, 238604899.58276495, 234736463.26037943, 230445168.40117618, 222781719.74012354, 212495709.71050605]
  hst_wave   = [0.43292, 0.47462, 0.53609, 0.5922, 0.6312, 0.76932, 0.8045, 0.90332, 0.23733, 0.27104, 0.28229, 0.33548, 0.39239, 0.43261, 0.47721, 0.4938, 0.53064, 0.58865, 0.62416, 0.76541, 0.80525, 0.91918, 0.98647, 1.0551, 1.15345, 1.24861, 1.39229, 1.53692]
  hst_filter = ['f435w', 'f475w', 'f555w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp', 'f225w', 'f275w', 'f300x', 'f336w', 'f390w', 'f438w', 'f475w', 'f475x', 'f555w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp', 'f098m', 'f105w', 'f110w', 'f125w', 'f140w', 'f160w']
  hst_camera = ['acswfc', 'acswfc', 'acswfc', 'acswfc', 'acswfc', 'acswfc', 'acswfc', 'acswfc', 'wfc3uvis', 'wfc3uvis', 'wfc3uvis', 'wfc3uvis', 'wfc3uvis', 'wfc3uvis', 'wfc3uvis', 'wfc3uvis', 'wfc3uvis', 'wfc3uvis', 'wfc3uvis', 'wfc3uvis', 'wfc3uvis', 'wfc3uvis', 'wfc3ir', 'wfc3ir', 'wfc3ir', 'wfc3ir', 'wfc3ir', 'wfc3ir']

  ; Original Kelsall98 fluxes
  dirbe_flux = [ 2.3405606e+08, 1.2309874e+08, 64292872., 35733824., 5763843.0, 1327989.4, 230553.73, 82999.336, 42346.605, 14409.608 ]
  dirbe_wave = [1.25,           2.2,           3.5,       4.9,       12,        25,        60,        100,       140,       240]
  
  wave_list = [hst_wave, dirbe_wave]
  flux_list = [hst_flux, dirbe_flux]
  
  ; Find closest wavelength
  diff_arr = abs(wave_list - wl)
  mn = (min(diff_arr, minval))
  idx = array_indices(diff_arr, minval)

;  print,'Sol flux info', wl, flux_list[idx]

  return, flux_list[idx]

end