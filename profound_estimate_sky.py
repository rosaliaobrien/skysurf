import numpy as np
import os
from astropy.io import fits
from measureskyregion_mean import measureskybin as measureskybin_mean
from measureskyregion import measureskybin
from split_image import bin_image

# FUNCTION TO WORK ON PROFOUND SKY MAPS where the rms data is in a separate extension
def calculate_sky_profound(sky_data, rms_data, bin_size = 64, dq_data = None, dq_good_list = [0], has_DQ = True, dq_fraction = 0.2, percentile = 50):

	'''
	Calculte the Percentile-clip sky surface brightness from a single science array from an Hubble Space Telescope flt/ flc image.
	Created for Project SKYSURF. Described in O'Brien et al. (2022): https://ui.adsabs.harvard.edu/abs/2022arXiv221008010O/abstract.

	Params
	------
	sky_data - arr
		Profound sky map
	rms_data - arr
		Profound sky rms map
	bin_size - int (optional)
		Size of sub-regions used to estimate the sky-SB. For SKYSURF, bin_size = 64 for ACS/WFC and WFC3/UVIS
		and bin_size = 39 for WFC3/IR
	dq_data - arr (optional)
		Data Quality data from 'DQ' extension
	dq_good_list - list (optional)
		List of DQ flags to be considered good when masking. MUST BE A BINARY NUMBER (0,1,2,4,8,16,etc)
	has_DQ - bool (optional)
		True if DQ will be used
	dq_fraction - float
		Fraction (between 0 and 1.0) of pixels flagged from DQ array within a sub-region for the subregion
		to be considered "bad" and subsequently masked
	percentile - int
		Percentile to be used when estimating the sky-SB. A value of 50 is recommended.

	Outputs
	-------
	'calculate_sky' returns a dictionary with the keys described below. Each key is within a list (so that creating pandas)
	dataframes is easier, so to extract a single value, a user must do dict[key][0].

	calc_sky - float
		The calculated sky-SB in the native units of the image
	calc_rms - float
		The calculated sky-SB rms in the native units of the image
	N_total_regions - int
		Number of sub-regions
	N_bad_regions - int
		Number of sub-regions where the algorithm detected an object
	N_bad_px_regions - int
		Number of sub-regions where the number of DQ-flagged pixels is greater than 20%
	N_good_regions - int
		Number of sub-regions used to calculated calc_sky and calc_rms
	calc_sky_mean - float 
		Same as calc_sky, but using a mean instead of a median during the iterative clipping
	calc_rms_mean - float 
		Same as calc_rms, but using a mean instead of a median during the iterative clipping
	sky_arr - arr
		Array containing all sky levels corresponding to individual sub-regions (used for creating plots)
	rms_arr - arr
		Array containing all sky rms levels corresponding to individual sub-regions (used for creating plots)
	cutouts - list
		List of Cutout2D objects corresponding to each sub-region (used for creating plots)
	lowest5perc_ind 
		Indices corresponding to the darkest 5% of sub-regions (used for creating plots)
	bad_ind - list
		Indices corresponding to sub-regions where the algorithm detected an object (used for creating plots)
	badpx_ind - list
		Indices corresponding to sub-regions where the number of DQ-flagged pixels is greater than 20% (used for creating plots)
	mean_x_pos and mean_y_pos - float
		Mean x/y position of all "good" sub-regions, in pixel units
	std_x_pos and std_y_pos - float
		Standard deviation of all "good" sub-regions, in pixel units
	'''

	sky_data = np.copy(sky_data)

	# Mask pixels flagged in DQ array 
	if has_DQ == True:
		# Define list of all possible bad flags
		dq_bad_list = 2**np.arange(0,20)

		# Loop trhough each bad flag and mask corresponding pixels in the SCIENCE array
		# Note: (dq_data | FLAG) returns dq_data if FLAG is within dq_bad_list
		# Example: (6|4) == 6
		# Example: (6|2) == 6
		for dq_bad in dq_bad_list:
			if dq_bad not in dq_good_list:
				sky_data[(dq_data | dq_bad) == dq_data] = np.nan

	### Make cutouts ###
	# cutout_shape, cutouts = bin_image(use_array=True, data_array=data, bin_size=(bin_size,bin_size),
	# 	bin_origin = 'lower left', show_image = False, ignore_borders = False, border = 0)

	#Define list of sky and RMS values so it can easily be appended
	all_skys = []
	all_rms = []

	# all_skys_mean = []
	# all_rms_mean = []
	
	#Define list of bad pixel regions (will be based on DQ array)
	badpx_ind = []

	#### SKY MAP #####

	# Make cutouts of SKY MAP
	cutout_shape, sky_cutouts = bin_image(use_array=True, data_array=sky_data, bin_size=(bin_size,bin_size),
		bin_origin = 'lower left', show_image = False, ignore_borders = False, border = 0)

	# Calculate subreigon information from SKY MAP
	for ci, c in enumerate(sky_cutouts):

		# Calculate sky and sky rms
		sky, sky_rms = measureskybin(c.data,axis=0)

		# If more than 20% of pixels are NaN (from DQ mask), then don't use this subregion
		if np.count_nonzero(np.isnan(c.data))/8 > (cutout_shape[0]*cutout_shape[1])*0.2:
			sky = float('nan')
			sky_rms = float('nan')
			badpxregs.append(ci)

		all_skys.append(sky)



	#### RMS MAP #####

	# Make cutouts of SKYRMS MAP
	cutout_shape, rms_cutouts = bin_image(use_array=True, data_array=rms_data, bin_size=(bin_size,bin_size),
		bin_origin = 'lower left', show_image = False, ignore_borders = False, border = 0)

	# Calculate subreigon information from SKY MAP
	for ci, c in enumerate(rms_cutouts):

		# Calculate the sky of the subregion (which is really the sky RMS of the subregion) and the RMS of the RMS
		rms, rms_rms = measureskybin(c.data,axis=0)


		all_rms.append(rms)

	# for ci, c in enumerate(cutouts):

	# 	sky, rms = measureskybin(c.data,axis=0)
	# 	sky_mean, rms_mean = measureskybin_mean(c.data,axis=0)

	# 	if np.count_nonzero(np.isnan(c.data)) > (cutout_shape[0]*cutout_shape[1])*dq_fraction:
	# 		sky = float('nan')
	# 		rms = float('nan')
	# 		badpx_ind.append(ci)

	# 		sky_mean = float('nan')
	# 		rms_mean = float('nan')

	# 	all_skys.append(sky)
	# 	all_rms.append(rms)

	# 	all_skys_mean.append(sky_mean)
	# 	all_rms_mean.append(rms_mean)

	#Define all_skyprms as an array of each regions sky value plus its RMS value
	all_skyprms = np.array(all_skys)+np.array(all_rms)
	bad_ind = []
	#Get list of bad indicies (that arent bad pixel values) and set these sky values to nan
	for sky_i, sky in enumerate(all_skys):
		if sky > np.nanmin(all_skyprms):
			bad_ind.append(sky_i)
			all_skys[sky_i] = float('nan')
			all_rms[sky_i] = float('nan')

	#Define sky and RMS arrays after finished editing them
	bkg_arr = np.array(all_skys)
	rms_arr = np.array(all_rms)

	#Good indices are where bkg array ISNT a NaN
	goodind = np.where(~np.isnan(bkg_arr))[0].tolist()

	#Define array without the nans
	bkg_arr_nonans = bkg_arr[goodind]
	calc_bkg = np.nanpercentile(bkg_arr_nonans,percentile) #Calc sky is Nth percentile of arr w/out nans
	calc_rms = np.nanmean(rms_arr) #Calc rms is just mean of all RMS values

	#Get INDICIES of lowest 5% of regions :)
	N_good = len(goodind) #Number of good regions
	N_5perc_of_good = int(N_good*0.05) #Number of lowest 5% of good regions
	lowest5perc_ind = np.argsort(bkg_arr)[:N_5perc_of_good] #Indices of lowest 5% of good regions

	#Find average x and y positions of useable subregions
	good_x_pos = [] #List of useable x positions
	good_y_pos = [] #List of useable y positions
	for c in np.array(sky_cutouts)[goodind]:
	    x,y = c.position_original #For each useable (not NaN) cutout, get position in original data array
	    good_x_pos.append(x)
	    good_y_pos.append(y)
	mean_x_pos = np.mean(good_x_pos)
	mean_y_pos = np.mean(good_y_pos)
	std_x_pos = np.std(good_x_pos)
	std_y_pos = np.std(good_y_pos)


	### Get MEAN sky value (for errors) ###
	# all_skyprms_mean = np.array(all_skys_mean)+np.array(all_rms_mean)
	# badind_mean = []
	# for sky_i, sky in enumerate(all_skys_mean):
	# 	if sky > np.nanmin(all_skyprms_mean):
	# 		badind_mean.append(sky_i)
	# 		all_skys_mean[sky_i] = float('nan')
	# 		all_rms_mean[sky_i] = float('nan')
	# bkg_arr_mean = np.array(all_skys_mean)
	# rms_arr_mean = np.array(all_rms_mean)
	# goodind_mean = np.where(~np.isnan(bkg_arr_mean))[0].tolist()
	# bkg_arr_nonans_mean = bkg_arr_mean[goodind_mean]
	# calc_bkg_mean = np.nanpercentile(bkg_arr_nonans_mean,percentile) #Calc sky is 5th percentile of arr w/out nans
	# calc_rms_mean = np.nanmean(rms_arr_mean) #Calc rms is just mean of all RMS values

	N_tot = len(sky_cutouts)
	N_bad = len(bad_ind)
	N_badpx = len(badpx_ind)
	N_good = N_tot-N_bad-N_badpx

	dic = {'calc_sky': [calc_bkg], 'calc_rms': [calc_rms], 
	       'N_total_regions': [N_tot], 'N_bad_regions': [N_bad], 'N_bad_px_regions': [N_badpx], 'N_good_regions': [N_good],  
	       'sky_arr': [bkg_arr], 'rms_arr': [rms_arr], 
	       'cutouts': [sky_cutouts], 'lowest5perc_ind': [lowest5perc_ind], 'bad_ind': [bad_ind], 'badpx_ind': [badpx_ind],
	       'mean_x_pos': [mean_x_pos], 'mean_y_pos': [mean_y_pos], 'std_x_pos': [std_x_pos], 'std_y_pos': [std_y_pos]}

	return dic
