import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
import numpy as np
from astropy.stats import sigma_clip

def make_plots(data, cutouts, goodind, badind, sky, rms, badpx = None, figsize = (10,10),
    title = None, save = False, savepath = None, show = False, text_pos = 0.12):
    '''
    Params
    ------
    data - arr
        Full array containing SCI image
    cutouts - list
        List of 2DCutout objects
    goodind - list
        List of indices for "good" cutout regions
    badind - list
        List of indices for "bad" cutout regions
    sky - float
        Sky-SB measurement
    rms - float
        Sky-SB rms measurement
    badpx - list
        List of indices where a specified fraction of pixels are flagged via the DQ array
    figsize - tuple
        Size of final figure
    title - str
        Title of figure
    save - bool
        True if figure will be saved
    savepath - str
        Filename of saved figure. Figure will only be saved is save == True
    show - bool
        Whether to display plot
    text_pos - float
        Defines the position under the main figure where text containing the sky, sky rms, and other 
        information will be plotted. This value represents the fraction of the total y-length of the figure.
    '''
    
    #Make sure path where image will save is specified if you are choosing to save the output
    if (save == True) & (savepath == None):
        raise Exception('Must specify savepath')

    fig, ax1 = plt.subplots(figsize = figsize)
    make_subplot(ax1, sky, rms, data, cutouts, goodind = goodind, badind = badind, badpx = badpx, text_pos = text_pos)

    #Add title to image
    if not title == None:
        fig.suptitle(title, y = 1.005, fontsize = 20)

    plt.tight_layout()
    
    if save == True:
        plt.savefig(savepath, bbox_inches = 'tight')
    if show == True:
        plt.show()
    if not show == True:
        plt.close()

def make_subplot(ax, skysurf_sky, skysurf_rms, data, cutouts, goodind = [], badind = [], badpx = [], text_pos = 0.12):

    # Save SKYSURF sky-SB and sky-SB RMS to show at bottom of plot
    calc_bkg = np.copy(skysurf_sky)
    calc_rms = np.copy(skysurf_rms)
    
    # If the median pixel value of the data is negative, make it positive
    # Run sigma_clip one time only to make the code more efficient
    clipped_data = sigma_clip(data, sigma = 3)
    clipped_median = np.ma.median(clipped_data)
    # if clipped_median < 0:
    #     data += np.abs(clipped_median)
    
    # Find the sky and rms values you will use for plotting
    plotting_data = np.abs(data-clipped_median)
    clipped_plotting_data = sigma_clip(plotting_data, sigma = 3)
    plot_sky = np.ma.median(clipped_plotting_data)
    plot_rms = np.ma.std(clipped_plotting_data)  
        
    # Set the minimum and maximum values for the colorbar
    minv = plot_sky-plot_rms
    maxv = plot_sky+10*plot_rms
        
    # Only plot if the data is not ALL NaNs
    if np.isnan(data).all() == False:
        ax.imshow(plotting_data, norm=LogNorm(vmin = minv, vmax = maxv), cmap='Greys', origin='lower')
    elif np.isnan(data).all() == True:
        ax.imshow(np.ones(np.shape(data)), cmap='Greys', origin='lower')
    
    #For each cutout object, plot it on the original image
    for ci, c in enumerate(cutouts):
        #Good regions are green
        if ci in goodind:
            c.plot_on_original(fill = False, color='green', label='Lowest 5% of good regions' if ci == goodind[0] else '', alpha = 0.5, ax = ax)
            c.plot_on_original(fill = True, facecolor='green', alpha = 0.09, ax = ax)
        #Bad regions are red
        if not np.shape(badind) == (0,): #Only run if there are bad regions
            if ci in badind:
                c.plot_on_original(fill = False, color='red', label='Bad regions' if ci == badind[0] else '', alpha = 0.5, ax = ax)
                c.plot_on_original(fill = True, facecolor='red', alpha = 0.09, ax = ax)
        if ci in badpx:
            c.plot_on_original(fill = False, color='purple', label='Regions with too many bad px' if ci == badpx[0] else '', alpha = 0.5, ax = ax)
            c.plot_on_original(fill = True, facecolor='purple', alpha = 0.09, ax = ax)
        #Neutral regions aren't colored
        if (not ci in goodind) and (not ci in badind) and (not ci in badpx):
            c.plot_on_original(color='white', alpha = 0.3, ax = ax) #White border

    reg_frac = len(badind)/len(cutouts)

    #Add target name, calculated sky, calculated RMS, and calculated gradient
    ax.text(0,-len(data)*text_pos,'sky: {:.3f} \nrms: {:.3f} \nF={:.3f}'.format(calc_bkg,calc_rms,reg_frac), fontsize = 15)

    #No tick marks
    ax.tick_params(labelbottom=False, labelleft = False)

    #Legend to indicate bad = red and good = green
    ax.legend(frameon = False, bbox_to_anchor=(1, -0.13), loc='lower right', fontsize = 13)

    if np.shape(data)[1] == 4096:
        ax.axvline(np.shape(data)[1]/2, color = 'black')
