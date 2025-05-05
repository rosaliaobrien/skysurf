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

end