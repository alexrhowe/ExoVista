import astropy.io.fits
import numpy as np
import scipy.interpolate

def load_scene(inputfile, time = 0):
    """
    This routine reads the output .fits file produced by exoVista
    and converts all quantities to flux at the same spectral resolution.
    Pixel coordinates of all objects are returned. Note that pixel coordinates
    are in 00LL format, where (0,0) is the lower-left corner of the lower-left
    pixel (and (0.5,0.5) is the center of the lower-left pixel). The star
    should be located at exactly (npix/2,npix/2)--the intersection of
    the central 4 pixels.

    --- INPUTS ---
    inputfile = filename and path of fits file containing scene
    time = desired time (default = 0)

    --- RETURNS --- 
    lambda: wavelength vector (microns)
    xystar: location of star (pixels)
    fstar: stellar flux vector (Jy)
    xyplanet: Nplanets x 2 array of planet locations (pixels)
    fplanet: Nplanets x nlambda array of planet fluxes (Jy)
    diskimage: npix x npix x nlambda disk image cube (Jy per pix)
    angdiam: angular diameter of star (mas)
    pixscale: pixel scale (mas)
    """
    
    # Constants defined by the FITS file structure
    hdul = astropy.io.fits.open(inputfile)
    if len(hdul) < 4:
        print('Error missing extensions in FITS file.')
        exit()
    notime = False
    if hdul[3].data.ndim == 1: notime==True
    nplanets = len(hdul)-4
    if(abs(hdul[-1].header['A']-1)<1.e-5): nplanets -= 1 # remove the extra Earth twin if it is present

    #Define extension numbers
    lam_ext = 0
    disklam_ext = 1
    disk_ext = 2
    star_ext = 3
    planet_ext = 4 #first planet extension
    h = astropy.io.fits.getheader(inputfile, ext=0) #read header of first extension
    n_ext = h['N_EXT'] #get the largest extension #

    #Get wavelength array
    lambdas, h = astropy.io.fits.getdata(inputfile, ext=lam_ext, header=True) #read wavelength extension
    nlambda = len(lambdas)
    
    #STEP 1: STAR
    #Need to determine x, y, and fstar
    xystar = np.zeros(2)
    fstar = np.zeros(nlambda)
    d, h = astropy.io.fits.getdata(inputfile, ext=star_ext, header=True)
    angdiam = h['ANGDIAM']
    pixscale = h['PXSCLMAS']
    
    if d.ndim == 1: d = np.expand_dims(d, 0)
    t = d[:,0] # time vector
    x = d[:,1] # heliocentric x location vector (pix)
    y = d[:,2] # heliocentric y location vector (pix)
    
    xystar[0] = x[0] # pick the first entry by default
    xystar[1] = y[0]
    fstar = d[0,15:15+nlambda] # grab the stellar flux of first time entry
    
    # If the fits file contains a vector of times, interpolate...
    if len(t) > 1: 
        x_interp = scipy.interpolate.interp1d(t, x, kind='quadratic')
        y_interp = scipy.interpolate.interp1d(t, y, kind='quadratic')
        xystar[0] = x_interp(time)
        xystar[1] = y_interp(time)
        
        for ii in range(nlambda):
            fstar_interp = scipy.interpolate.interp1d(t, d[:,15+ii], kind='quadratic')
            fstar[ii] = fstar_interp(time)
            
    #STEP 2: PLANETS
    #;Need to determine x, y, and fplanet
    xyplanet = np.zeros((nplanets,2))
    fplanet = np.zeros((nplanets,nlambda))
    for ip in range(nplanets): #loop over all planets
        d, h = astropy.io.fits.getdata(inputfile, ext=planet_ext+ip, header=True)
        if d.ndim == 1: d = np.expand_dims(d, 0)
        t = d[:,0] # time vector
        x = d[:,1] # heliocentric x position vector (pix)
        y = d[:,2] # heliocentric y position vector (pix)
        
        xyplanet[ip,0] = x[0] # pick the first entry by default 
        xyplanet[ip,1] = y[0]
        contrast = d[0,15:15+nlambda]
        fplanet[ip,:] = contrast * fstar #convert to flux
                
        if len(t) > 1:
            x_interp = scipy.interpolate.interp1d(t, x, kind='quadratic')
            y_interp = scipy.interpolate.interp1d(t, y, kind='quadratic')
            xyplanet[ip,0] = x_interp(time)
            xyplanet[ip,1] = y_interp(time)
            
            for ii in range(nlambda):
                contrast_interp = scipy.interpolate.interp1d(t, d[:,15+ii], kind='quadratic')
                contrast = contrast_interp(time)
                fplanet[ip, ii] = contrast * fstar[ii]
                
    #STEP 3: DISK
    lambdas_disk = astropy.io.fits.getdata(inputfile, ext=disklam_ext) # disk wavelengths
    nlambdas_disk = len(lambdas_disk)
    temp = astropy.io.fits.getdata(inputfile, ext=disk_ext)
    contrast = temp[0:nlambdas_disk,:,:] #3D contrast data cube
    cprecision = temp[nlambdas_disk,:,:] #2D contrast precision
    
    #Interpolate the disk image cube to the desired wavelength spacing
    lambda_indices = np.searchsorted(lambdas_disk, lambdas) - 1
    
    #index in log lambda space (have to add on fractional indices)
    frac_lambda_indices = (lambda_indices + 
            (np.log(lambdas) - np.log(lambdas_disk[lambda_indices])) /
            (np.log(lambdas_disk[lambda_indices + 1])
             - np.log(lambdas_disk[lambda_indices])))
    
    contrast_interp = scipy.interpolate.interp1d(np.arange(len(lambdas_disk)), contrast, axis=0, kind='cubic')
    diskimage = np.multiply(contrast_interp(frac_lambda_indices).T, fstar).T
    
    return (lambdas, xystar, fstar, xyplanet, fplanet, diskimage, angdiam, pixscale)
