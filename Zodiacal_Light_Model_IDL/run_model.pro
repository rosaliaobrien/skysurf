; Generate random arrays of ecliptic latitude, ecliptic longitude, and day number
random_arr = randomu(seed, 100)
lat = random_arr*180-90
lon = random_arr*360
day = random_arr*365

; Predict zodiacal light for a specific wavelength
wavelength = 1.25

; Run the zodiacal light model, using the SKYSURF prescription and the Kelsall prescription
zodi_result_skysurf = get_zmod(wavelength,"skysurf",day,lon,lat)
zodi_result_kelsall = get_zmod(wavelength,"kelsall",day,lon,lat)

; Plot and compare the two models
plt = SCATTERPLOT(lat,zodi_result_skysurf,XTITLE='Ecliptic Latitude [deg]',YTITLE='Zodiacal Light [MJy sr!u-1!n]', $
                  font_size = 15, SYM_COLOR = 'blue')
plt = SCATTERPLOT(lat,zodi_result_kelsall, overplot = 1, SYM_COLOR = 'red')





; Try running it on a range of wavelengths
; 
; Note that the SKYSURF version is only reliable for HSTs nominal wavelengths (0.3--1.7 microns) 
; and Sun angles greater than 80 deg

; Define array of wavelengths
wave_arr = [0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0]

; Initialize empty array to populate zodi values
zodi_skysurf = []

; Loop through all wavelengths
FOR i = 0, N_ELEMENTS(wave_arr)-1 DO BEGIN
  wave = wave_arr[i]

  ; Assume day = 100, lon = 100, lat = 90
  z_skysurf = get_zmod(wave, "skysurf", 100, 180, 90)

  zodi_skysurf = [zodi_skysurf, z_skysurf]
ENDFOR

; Plot and compare the two models
plt = plot(wave_arr,zodi_skysurf,XTITLE='Ecliptic Latitude [deg]',$
  YTITLE='Zodiacal Light [MJy sr!u-1!n]', font_size = 15, color = 'green')






; Plot the isotropic-like cloud!

isocloud = []
; Loop through all wavelengths
FOR i = 0, N_ELEMENTS(wave_arr)-1 DO BEGIN
  wave = wave_arr[i]

  ; Assume day = 100, lon = 100, lat = 90
  iso_skysurf = get_zmod(wave, "skysurf", 100, 180, 90,new_iso_comp=1,iso_comp_only=1)

  isocloud = [isocloud, iso_skysurf]
ENDFOR

; Plot and compare the two models
plt = plot(wave_arr,isocloud,overplot = 1, color = 'purple')


end