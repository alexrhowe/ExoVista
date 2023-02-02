import numpy as np
import pandas as pd
import os
from astropy.io import fits
from datetime import datetime
import matplotlib.pyplot as plt
from src.constants import *
from src import Settings

settings = Settings.Settings()
rng = np.random.default_rng()

def generate_disks(stars, planets, settings, nexozodis=None, rand_components=False):

    settings = settings
    rng = np.random.default_rng(settings.seed)
    
    nstars = len(stars)
    if settings.ncomponents < mincomponents or settings.ncomponents > maxcomponents:
        print('ERROR: Maximum number of disk components is set to {0:d}-{1:d}'.format(mincomponents,maxcomponents))
        return

    if rand_components:
        ndisks = rng.integers(mincomponents,maxcomponents,nstars)
    else:
        ndisks = np.full(nstars,settings.ncomponents)
    disks = np.zeros((nstars,maxcomponents,len(dlabel)))
    disks[:,:,dlabel.index('n')] = -1
    compositions = np.full(nstars,'')
    
    '''
    i: inclination of disk relative to midplane
    longnode: longitude of ascending node of disk
    nzodis: density of each disk component (3 components max)
    r: radial location of each dust belt component
    dror: width of belt
    rinner: inner truncation radius of each component
    eta: PR drag time divided by collision time scale for blowout size
    hor: scale height of each component
    NOTE: each disk component is assigned a linear combination of
    3 Henyey-Greenstein scattering phase functions.in
    g: HG asymmetry parameters (3 values for each component)
    w: HG weighting factors for each function and component
    '''
  
    # Distribute exozodi levels
    if nexozodis is not None:
        nzl = nexozodis
    else:
        with fits.open('nominal_maxL_distribution-Dec2019.fits') as data: nexozodi_levels = pd.DataFrame(data[0].data)
        nzl = np.array(nexozodi_levels[nexozodi_levels.columns[0]].values)
        rng.shuffle(nzl,axis=0)  # Deterministic version for testing
    
    # all disks assigned the same orientation as system midplane
    disks[:,:,dlabel.index('longnode')] = 0.
    disks[:,:,dlabel.index('i')] = 0.

    if len(nzl) < nstars:
        print('Error: exozodi file too short for target list.')
        return
    disks[:,0,dlabel.index('nzodis')] = nzl[0:nstars] * 0.642
    # factor of 0.642 accounts for the fact that cold components will also contribute dust to reproduce the LBTI distribution
    
    maxfac = np.log10(settings.density_ratio)    # outer disks must be within some density factor of inner disk
    minfac = np.log10(1./settings.density_ratio)
    for i in range(0,nstars):
        if ndisks[i] > 1:
            for icomp in range(1,ndisks[i]):
                disks[:,icomp,dlabel.index('nzodis')] = 10.**(rng.random(nstars)*(maxfac-minfac)+minfac)*disks[:,0,dlabel.index('nzodis')]
    
    '''
    Assign the disks a composition.
    Currently we assign the same composition to all components, as 
    thread_the_scene currently calls load_lqsca once for all components.
    An update to randomize composition for each component would require
    restructuring thread_the_scene and disk_contrast_image, and would slow
    things down a bit.
    '''

    lqqdir = 'lqq_files/'
    complist = []
    for root, dirs, files in os.walk(lqqdir):
        for dir0 in dirs:
            complist.append(dir0)
    complist = np.array(complist)
    randcomplist = complist[rng.integers(0,len(complist)-1,nstars)]
    compositions = randcomplist
    
    # Define spatial distribution of dust
    # We do this on a star by star basis
    # This allows us to generate rules based on the underlying planetary system
    for i in range(0,nstars):
        
        s = stars.loc[i]
        p = planets[i]
        nplanets = len(np.where(p[:,pllabel.index('R')] > 0)[0])

        # Enforce Hill sphere stability rules
        # Calculate Hill radii of all planets
        # Note that all missing planets will have a Hill radius = 0
        hill_radius = p[:,pllabel.index('a')] * (1.-p[:,pllabel.index('e')]) * (p[:,pllabel.index('M')] / (3 * s['mass'] * 333000.))**(1./3.) # AU
        hinner = p[:,pllabel.index('a')] - hill_radius * settings.stability_factor
        houter = p[:,pllabel.index('a')] + hill_radius * settings.stability_factor
        j = [j for j in range(0,len(p)) if hinner[j]>0 and p[j,pllabel.index('M')]>1.e-3]
        hinner = hinner[j]
        houter = houter[j]
        
        # Convert from Hill spheres to gaps between planets, including inside and outside the planets
        minr = min(settings.r_min)
        maxr = max(settings.r_max)
        ginner = np.append(0,houter)
        gouter = np.append(hinner,1.e20)
        
        # Accounting for the overall min and max r for dust
        binner = []
        bouter = []
        for ip in range(0,len(ginner)):
            if gouter[ip] < minr: pass
            elif ginner[ip] > maxr: break
            elif ginner[ip] < minr and gouter[ip] > maxr:
                binner = [minr]
                bouter = [maxr]
                break
            elif ginner[ip] < minr and gouter[ip] > minr:
                binner.append(minr)
                bouter.append(gouter[ip])
            elif ginner[ip] < maxr and gouter[ip] > maxr:
                binner.append(ginner[ip])
                bouter.append(maxr)
            else:
                binner.append(ginner[ip])
                bouter.append(gouter[ip])
        binner = np.array(binner)
        bouter = np.array(bouter)

        # Now we know where we can put the dust
        # Let's distribute each component
        # First, the easy stuff that doesn't depend on the system architecture...
        for j in range(0,ndisks[i]):
            disks[i,j,dlabel.index('n')] = j
            disks[i,j,dlabel.index('hor')] = rng.random() * (settings.hor_max - settings.hor_min) + settings.hor_min
            disks[i,j,dlabel.index('g0')] = rng.random() * (0.995-0.8) + 0.8  # most forward scatt term
            disks[i,j,dlabel.index('g1')] = rng.random() * (0.8-0.35) + 0.35  # medium forward scatt term
            disks[i,j,dlabel.index('g2')] = rng.random() * (0.35+0.3) - 0.3   # least forward scatt term
            disks[i,j,dlabel.index('w0')] = rng.random() * (0.8-0.4) + 0.4    # most forward scatt term
            minw = 0.1
            maxw = 1.0 - disks[i,j,dlabel.index('w0')]
            disks[i,j,dlabel.index('w1')] = rng.random() * (maxw-minw) + minw # medium forward scatt term
            disks[i,j,dlabel.index('w2')] = 1. - disks[i,j,dlabel.index('w0')] - disks[i,j,dlabel.index('w1')] # least forward scatt term
                
            # Now the more complicated part...
            # Defining the belt based on the largest region of stability, max(dr/r)
            # Preferentially produces belts that are too wide.
            # Instead, select a logarithmically random point within the gaps
            # Subject to a minimum and maximum belt width
            dror = 0
            r0 = 0
            lbinner = np.log10(binner)
            lbouter = np.log10(bouter)
            measure = np.sum(lbouter-lbinner)
            n = 0
            
            while (dror < settings.dror_min) or (r0 < settings.r_min[j] or settings.r_max[j] < r0) and n<50:
                rand = rng.random()*measure
                ip = 0
                while ip < len(binner):
                    step = lbouter[ip]-lbinner[ip]
                    if rand <= step:
                        r0 = binner[ip] * 10**rand
                        break
                    else:
                        rand -= step
                        ip += 1

                dror = max(r0-binner[ip], bouter[ip]-r0)/r0
                n += 1

            dror = min(dror,settings.dror_max)

            # Note: these do not determine the edges of the disk, but breaking points in the geometry.
            disks[i,j,dlabel.index('r')] = r0
            disks[i,j,dlabel.index('dror')] = dror
            disks[i,j,dlabel.index('rinner')] = 0.0
                
            # Determine if a massive planet could truncate the inward migrating component
            k = np.array( [ jj for jj in range(0,len(p)) if p[jj,pllabel.index('a')] < disks[i,j,dlabel.index('r')] and p[jj,pllabel.index('M')] > settings.rinner_mass_threshold ] )
            if len(k)>0: disks[i,j,dlabel.index('rinner')] = min( (1.1 * np.max( p[k,pllabel.index('a')] * (1+p[k,pllabel.index('e')]) )), disks[i,j,dlabel.index('r')] )
            
            # If the while loop failed to find a valid r, set density to zero
            if dror < settings.dror_min: disks[i,j,dlabel.index('nzodis')] = 0.
            
            # Now that we have the locations we can determine collision rates...
            beta = 0.5
            alpha = beta * grav_const * s['mass'] / c                         # quantity related to PR drag to be used later
            tpr = disks[i,j,dlabel.index('r')]**2 / (4 * alpha)               # yr
            torbit = np.sqrt(disks[i,j,dlabel.index('r')]**3 / s['mass'])     # yr
            tcoll = torbit / (4*np.pi*disks[i,j,dlabel.index('nzodis')]*1e-7) # yr            
            disks[i,j,dlabel.index('eta')] = tpr / tcoll
    
    print('generate_disks() done')
    return disks, compositions
