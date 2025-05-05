;Code to confirm that the new model works

; Define filter you want to fit to
camera = 'wfc3ir'
flter_str = 'wfc3ir/f160w'
flter_file = '/Users/rosaliaobrien/skysurf/zodi_modeling/scripts/for_hg_fitting_tables/wfc3ir_F160W.csv'
fit_type = 'hg4'

; Read skysurf measurements
skysurf_skys = READ_CSV(flter_file, HEADER=skysurf_hdr)

; Check that the columns are where you think they are
index_sky = WHERE(STRCMP(skysurf_hdr, 'sky-dgl'))+1
index_lat = WHERE(STRCMP(skysurf_hdr, 'Ecliptic_Latitude'))+1
index_err = WHERE(STRCMP(skysurf_hdr, 'Error_MJy/sr'))+1
if index_sky ne 48 then begin
  print,'!!! Columns are out of order !!!'
  return
endif
if index_lat ne 23 then begin
  print,'!!! Columns are out of order !!!'
  return
endif
if index_err ne 44 then begin
  print,'!!! Columns are out of order !!!'
  return
endif

; Assign certain columns to variables
lat = skysurf_skys.FIELD23
lon = skysurf_skys.FIELD24
day = skysurf_skys.FIELD46
skys = skysurf_skys.FIELD48 ; FIT to the skys MINUS DGL column
sky_err = skysurf_skys.FIELD44
sunang = skysurf_skys.FIELD30

zodi_result = get_zmod(flter_str,day,lon,lat)
zodi_final = zodi_result
diff = skys-zodi_final
;plt = SCATTERPLOT(sunang,diff,XTITLE='Sun Angle',YTITLE='Sky - Zodi [MJy sr!u-1!n]',yrange=[-0.03,0.1])

plt = SCATTERPLOT(sunang,diff,XTITLE='Sun Angle [deg]',YTITLE='HST Sky - Model Zodi [MJy sr!u-1!n]', $
                  xrange=[70,190],yrange=[-0.03,0.1], font_size = 15, margin = fig_margin)
plt = plot([70, 200], [0, 0], color='grey', thick=2, overplot = 1)
plt = plot([70, 200], [0.02, 0.02], color='grey', thick=2, overplot = 1)
std = stddev(diff)
t = text(150, 0.08, 'std = '+string(std, format = '(f0.4)'), /data, font_size = 15, color = 'red')

; Compare with other methods

day = findgen(365)
lon = 10
lat = 80

zodi_new = get_zmod('wfc3ir/f125w',day,lon,lat)
zodi_old = get_zmod(1.25,day,lon,lat)
zodi_rick = get_zmod('BAND1',day,lon,lat)

plt = plot(day, zodi_new, color='green', ytitle = 'Zodi Intensity [MJy/sr]', xtitle = 'Day', font_size = 15, xrange = [0,365])
plt = plot(day, zodi_old, color='red', overplot = 1)
plt = plot(day, zodi_rick, color='blue', overplot = 1)

t = text(160, 0.13, "Rick's BAND1 Test", /data, font_size = 15, color = 'blue')
t = text(160, 0.125, 'Original Kelsall+98 model', /data, font_size = 15, color = 'red')
t = text(160, 0.12, "Rosalia's New Model", /data, font_size = 15, color = 'green')

end