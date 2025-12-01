# SKYSURF Zodiacal Light Model

This is an IDL model of the SKYSURF Zodiacal Light model (see O'Brien et al. 2025). Please refer to run_model.pro for example usage.

To run:
```
zodi_intensity = get_zmod(wavelength,"skysurf",day,lon,lat)
```

### Inputs:
````
wavelength  - [float] wavelength in microns
phase_type  - [str] "kelsall" or "skysurf"
				Whether to use the standard Kelsall+1998 phase function 
				or the O'Brien+2025 phase function. To use the Kelsall 
				phase function, the wavelength must be one of DIRBE's 
				nominal wavelengths: 1.25, 2.2, 3.5, 4.9, 12, 25, 60, 
				100, 140, 240
day        - [arr] 1990 day number(s), where 1.0 = 1 Jan 1990
lon        - [arr] ecliptic longitude(s)
lat        - [arr] ecliptic latitude(s)

note: lambda, day, lon, lat can be any mixture of scalars or arrays, 
           but all arrays must be of same length.
````
### Outputs:
`````
zodi - scalar or array of zodi model intensity in MJy sr-1
		For comparison with DIRBE data, this is 'quoted' intensity
		at the nominal wavelength(s) for an assumed source
		spectrum nu\*F_nu = constant
`````

### Optional Keyword Inputs:

````
zpar=zpar - array of zodi model params.  If not specified, the
			default will be restored from zpars.xdr
solar_irr - [float] solar irradiance (in MJy/sr) corresponding to a 
			bandpass
````

### Caveats

This model was written so that it can work between wavelengths 0.2--4.9 microns. 

For ease of use, the solar spectrum incorporated into this model has a set bandwidth of 0.1 microns. For the best accuracy or finer spectral resolution, we highly recommend that the user inputs their own spectral irradiance (in units of MJy) for each wavelength bin. 

The SKYSURF model is only reliable at:
	- wavelengths between 0.2 and 2.7 microns.
	- Sun angles greater than 80 deg.
