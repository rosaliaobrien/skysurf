from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.nddata.utils import Cutout2D

def bin_image(file_path = None, ext = None, use_array = False, data_array = None, bin_size=None, 
    bin_number=None, bin_origin = 'lower left', fig_dimensions = (6,6), show_image = True, is_segm =  False, 
    segm_object = None, copy = True, ignore_borders = False, border = 0, saveplot = False):
    
    '''
    Bins fits images into bins of certain sizes of a specific number of bins then plots the resulting bins.
    UPDATED FOR UVIS but still works on IR! Aka new and improved version!!
    
    Parameters
    -----------
    file_path - str
        File path of fits file
    ext - int
        Extension of fits file that contains the data
    bin_size - tuple (int)
        (ysize, xsize) of bins
    bin_number - tuple (int)
        (y_number, x_number) of bins
    bin_origin - str 
        Either 'lower left', 'upper right', or 'center' to specify where the bins align
        
    Outputs
    -------
    data - arr
        Numpy array of original data
    cutouts - list
        List of Cutout2D objects that are copies of the original data
    '''
    
    if use_array == True:
        data = data_array
        
    else:
        #Read file and data from file
        file = fits.open(file_path)
        data = file[ext].data
    
    #Raise exception if both bin size and bin number are defined
    if (not bin_size == None) and (not bin_number == None):
        raise Exception('Choose only bin size of bin number, not both.')
        
    if ignore_borders == False:
        border = 0

    img_xsize = len(data[0])-border*2
    img_ysize = len(data)-border*2
    
    #Based on wanted bin size, round DOWN to the number of bins that is closest to that bin size but still fills the image
    #Use leftover pixels for final bin
    if not bin_size == None:

        # print(img_xsize, img_ysize)
  
        num_xbins = math.floor(img_xsize/bin_size[1])
        num_ybins = math.floor(img_ysize/bin_size[0])
        xsize = bin_size[1]
        ysize = bin_size[0]

        leftover_x_px = (img_xsize % xsize)
        leftover_y_px = (img_ysize % ysize)

        # print('Num x bins = {}, Num y bins = {}'.format(num_xbins, num_ybins))
        # print('Leftover x px: {}, Leftover y px: {}'.format(leftover_x_px, leftover_y_px))
        
    #Find size of bins based on bin number
    if not bin_number == None:
        num_xbins = bin_number[1]
        num_ybins = bin_number[0]
        xsize = img_xsize/bin_number[1]
        ysize = img_ysize/bin_number[0]
        bin_origin = 'lower left'
    
    #First cutout aligns with lower left corner
    if bin_origin == 'lower left':
        cutouts = [] #Make list of cutouts
        xstart = xsize/2 #Start cutout (center coord of cutout) so the left side of the cutout aligns with x=0
        for xbin_index in range(0, num_xbins):

            if xbin_index == (num_xbins - 1): #If on last x bin, take leftover x px into account!
                xsize += leftover_x_px #Size of last xbin includes leftover px
                xstart += leftover_x_px/2 #Since xstart denotes the x-center of this bin, add half of the leftover px to the starting x location

            ystart = ysize/2 #Start cutout so the bottom of the cutout aligns with y=0

            for ybin_index in range(0, num_ybins):

                if ybin_index != (num_ybins-1): #If NOT on last y bin, continue as normal
                    # print(ystart)
                    cutout = Cutout2D(data, (xstart, ystart), (ysize, xsize), copy=copy, mode='strict') #Make cutout
                    ystart += ysize #Redefine position of cutout for next cutout
                    cutouts.append(cutout)

                if ybin_index == (num_ybins-1): #If on last y bin, take leftover y px into account!
                    # print(ystart)
                    cutout = Cutout2D(data, (xstart, ystart+leftover_y_px/2), (ysize+leftover_y_px, xsize), copy=copy, mode='strict') #Make cutout
                    cutouts.append(cutout)

            xstart += xsize


       
    #First cutout aligns with upper right corner
    if bin_origin == 'upper right':
        xstart = img_xsize-xsize/2 
        cutouts = [] #Make list of cutouts
        for xbin_index in range(0, num_xbins):

            if xbin_index == (num_xbins - 1): #If on last x bin, take leftover x px into account!
                xsize += leftover_x_px
                xstart -= leftover_x_px/2

            ystart = img_ysize-ysize/2 #Start cutout so the top side of the cutout aligns with y=ymax

            for ybin_index in range(0, num_ybins):

                if ybin_index != (num_ybins-1): #If NOT on last y bin, continue as normal
                    cutout = Cutout2D(data, (xstart, ystart), (ysize, xsize), copy=copy, mode='strict') #Make cutout
                    ystart -= ysize #Redefine position of cutout for next cutout
                    cutouts.append(cutout)

                if ybin_index == (num_ybins-1): #If on last y bin, take leftover y px into account!
                    cutout = Cutout2D(data, (xstart, ystart-leftover_y_px/2), (ysize+leftover_y_px, xsize), copy=copy, mode='strict') #Make cutout
                    cutouts.append(cutout)

            xstart -= xsize
            
#     #Cutouts are centered on image
#     if bin_origin == 'center':
#         cutouts = []
        
#         if ignore_borders == False:
#             xstart = (img_xsize % xsize)/2 - xsize/2
#         xrange = num_xbins+1
#         yrange = num_ybins+1
#         if ignore_borders == True:
#             xstart = border
# #            xrange = num_xbins-2
# #            yrange = num_ybins-2


#         # print(xstart)
        
#         for xbin_index in range(0, xrange):
        
#             if ignore_borders == False:
#                 ystart = (img_ysize % ysize)/2 - ysize/2
#             if ignore_borders == True:
#                 ystart = border
            
#             for ybin_index in range(0, yrange):
# #                print(xstart, ystart)
                
#                 if ignore_borders == False:
#                     cutout = Cutout2D(data, (xstart, ystart), (ysize, xsize), copy=True, mode='partial')
#                 if ignore_borders == True:
#                     cutout = Cutout2D(data, (xstart, ystart), (ysize, xsize), copy=True, mode='partial')
#                 ystart += ysize
#                 cutouts.append(cutout)
#             xstart += xsize
          
    if show_image == True:
        if is_segm == False:
            plt.figure(figsize=fig_dimensions)
            maskimg2 = np.ravel(data)
            maskimg2 = maskimg2[maskimg2 < 0.5] 
            minv2 = np.percentile(data,5)
            maxv2 = np.percentile(data,95)
            plt.imshow(data, vmin=minv2, vmax=maxv2, cmap='Greys_r', origin='lower')  
            for c in cutouts:
                c.plot_on_original(color='black')
            if saveplot == True:
                plt.savefig('bin_uvis.pdf')
            plt.show()
        if is_segm == True:
            plt.figure(figsize=fig_dimensions)
            plt.imshow(segm_object, cmap=segm_object.cmap(random_state=12345), origin='lower')  
            for c in cutouts:
                c.plot_on_original(color='white')
            plt.show()
    
    return (ysize, xsize), cutouts
