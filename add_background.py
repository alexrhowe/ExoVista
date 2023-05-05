from astropy.io import fits
import numpy as np
from scipy import interpolate
import load_scene
import matplotlib.pyplot as plt
import os.path
from src.constants import *
from src import Settings
from src import generate_scene

def add_background(inputfile, time = 0):
    """
    This routine reads the output .fits file produced by exoVista

    --- INPUTS ---

    --- RETURNS ---
    """
    
    scene_data = load_scene.load_scene(inputfile, time)
    # scene_data = (lambdas, xystar, fstar, xyplanet, fplanet, diskimage, angdiam, pixscale)

    hdul = fits.open(inputfile)
    
    version = 0.
    if 'VERSION' in hdul[0].header:
        version = hdul[0].header['VERSION']

    if version <= 2.1:
        specstart = 15
        hstar = 3
    else:
        specstart = 16
        hstar = 4
    
    dist = hdul[hstar].header['DIST']
    parallax = 1./dist

    gallist = ['GALAXIES_10lat0.30-0.37um.fits', 'GALAXIES_10lat0.37-0.46um.fits',
               'GALAXIES_10lat0.46-0.57um.fits', 'GALAXIES_10lat0.57-0.70um.fits',
               'GALAXIES_10lat0.70-0.87um.fits', 'GALAXIES_10lat0.87-1.07um.fits']

    for g in gallist:
        if not os.path.isfile(g):
            print('Error Haystacks background files not found.')
            exit()
    
    galfits0 = fits.open("GALAXIES_10lat0.30-0.37um.fits")
    # It's not internally documented, but these cubes have a pixel scale of 9 mas.
    
    ngalwave = len(galfits0[-1].data)
    galpix = (galfits0[1].shape)
    galpixscale = 9.0
    
    evwave = hdul[1].data[1:]
    disk_data = scene_data[5]
    
    # Planned extension to allow upsampling the disk spectral resolution to the stellar resolution
    '''
    if res == 'high':
        evwave = scene_data[0]
        disk_data = upsample_disk(inputfile, scene_data[5])
    else:
        evwave = hdul[1].data[1:]
        disk_data = scene_data[5]
    '''
    evpix = scene_data[5][0].shape
    evpixscale = scene_data[-1]
    
    clippix = int(np.ceil(evpix[0]*evpixscale/galpixscale))+1
    
    print('Interpolating extragalactic background.')
    
    background_limit = galpix[0] - clippix
    start_x = np.random.randint(0,background_limit)
    start_y = np.random.randint(0,background_limit)
    
    galwave = np.zeros(ngalwave*6)
    background_clip = np.zeros((ngalwave*6,clippix,clippix))
    for i in range(0,len(gallist)):
        galfits = fits.open(gallist[i])
        galwave[i*ngalwave:(i+1)*ngalwave] = galfits[-1].data
        for j in range(0,ngalwave):
            background_clip[i*ngalwave+j,:,:] = galfits[i+1].data[start_x:start_x+clippix,start_y:start_y+clippix]
    
    # not yet accounting for parallax and proper motion
    
    # need to use 2D interpolation, which can be tricky with the way the method is written
    # zero point of the subgrid is the same
    # interpolate to the higher resolution
    
    background_temp = np.zeros((ngalwave*6,evpix[0],evpix[1]))
    
    gal_points = (galwave, np.linspace(0,clippix-1,clippix), np.linspace(0,clippix-1,clippix))
    
    ev_grid_convert = (evpix[0]-1)*evpixscale/galpixscale
    ev_points = (evwave, np.linspace(0,ev_grid_convert,evpix[0]), np.linspace(0,ev_grid_convert,evpix[0]))
    
    f = interpolate.RegularGridInterpolator(gal_points, background_clip)
    
    background_final = np.zeros((len(evwave),evpix[0],evpix[1]))
    for i in range(0,len(ev_points[1])):
        if i%16==0: print(i)
        for j in range(0,len(ev_points[2])):
            background_final[:,i,j] = f((ev_points[0],ev_points[1][i],ev_points[2][j])) * (evpixscale/galpixscale)**2
            
    # calculate expected number of background stars in frame; usually << 1; requires the "spectra.fits" Haystacks model
    '''
    starfits = fits.open("spectra.fits")
    totalstars = len(starfits[0].data[2])
    nstars = totalstars * (evpix[0]*evpixscale/galpixscale / galpix[0])**2
    '''
    
    # Test plotting routine
    '''
    fig = plt.figure(figsize=(10.8,10.8))
    ax = fig.add_subplot(111)
    ax.set_xlim(0,250)
    ax.set_ylim(0,250)
    plt.subplots_adjust(left=0.001, bottom=0.001, right=0.999, top=0.999, wspace=0, hspace=0)
    xc = np.arange(0,251)
    yc = np.arange(0,251)
    ax.pcolor(xc,yc,np.log10(background_final[-1]+disk_data[-1]),cmap='inferno')
    plt.show()
    '''

    return background_final
