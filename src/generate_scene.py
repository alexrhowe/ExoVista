import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from os import path
from os import mkdir
from os import system
from astropy.io import fits
from datetime import datetime
import cython
from src import coordinates as coord
from src import nbody
from src import wrapImage
from src.constants import *
from src import Settings

settings = Settings.Settings()

def load_lqsca(lqq_dir, composition, rdust, rdust_boundaries, lam):    
    '''
    lqq_dir: directory (including /) of lqq files
    composition: composition of dust (string)
    rdust: grain size in microns
    rdust_boundaries: upper and lower limits of rdust bin in microns
    lam: desired wavelengths in microns
    Qsca: returned nsizes x nlambda array
    '''
    
    nsizes = len(rdust)
    nlambda = len(lam)
    Qsca_array = np.zeros((nsizes,nlambda))

    for i in range(0,nsizes):
    
        lqq_file = lqq_dir + composition + '/' + composition + '-s_{0:.4f}-smin_{1:.4f}-smax_{2:.4f}-dnds~s^-3.5.lqq'.format(rdust[i],rdust_boundaries[i],rdust_boundaries[i+1])
        
        # read the file
        fin = open(lqq_file,'r')
        head = fin.readline()
        fin.readline()
        d = fin.readlines()
        dlen = len(d)

        l = np.zeros(dlen)
        qsca = np.zeros(dlen)
        for j in range(0,dlen):
            line = d[j].split()
            l[j] = float(line[0])
            qsca[j] = float(line[2])
            if np.isnan(qsca[j]): qsca[j] = qsca[j-1] # Handles errors that appear in the last line of some of the files.

        # interpolate to desired wavelengths
        f = interp1d(l,qsca,kind='cubic')
        Qsca_array[i] = f(lam)
        
        if np.isnan(Qsca_array[i,0]):
            print('Error in Qsca file:')
            print(Qsca_array[i])
            print(qsca)
            print(i, lqq_file)
            exit()
    return Qsca_array
        
def lambertian(beta):
    # Returns the value of the Lambert phase function
    # beta = phase angle in radians
    # See, e.g., R. Brown (2005) "Sing-Visit Photometric
    # and Obscurational Completeness" Eq 4
    
    phi = (np.sin(beta) + (np.pi-beta) * np.cos(beta)) / np.pi
    return phi

def generate_scene(stars, planets, disks, albedos, compositions, settings):
    '''
    Some important conventions on coordinates:
    i = +/- 90 is an edge-on orientation
    +x points to right (W)
    +y points up (N)
    +z coordinate points toward the observer (out of page)
    '''
    settings = settings
    
    pixscale_arcsec = settings.pixscale
    if settings.output_dir == 'output': print('Default output directory: ./output/')
    if settings.output_dir[-1] != '/': settings.output_dir = settings.output_dir + '/'
    if not path.isdir(settings.output_dir): mkdir(settings.output_dir)

    print('Generating scenes. (This may take a while.)')

    # Define time array
    timemin = 0.
    if settings.timemax<=0.: settings.timemax = 1.e-10
    ntimes = int(np.ceil((settings.timemax-timemin)/settings.dt)) + 1
    time_vector = np.linspace(timemin,settings.timemax,ntimes)

    # Define master wavelength array
    dlnlambda = 1./settings.specres
    lnlambdamax = np.log(settings.lambdamax)
    lnlambdamin = np.log(settings.lambdamin)
    nlambda = int(np.ceil((lnlambdamax-lnlambdamin)/dlnlambda))
    lam = np.linspace(lnlambdamin,lnlambdamax,nlambda)
    lam = np.exp(lam)
    transition_lambda0 = (lam[0]/lam[1])*lam[0]
    transition_lambda2 = (lam[-1]/lam[-2])*lam[-1]
    transition_lambda_temp = np.insert(lam,0,transition_lambda0)
    transition_lambda_temp = np.insert(transition_lambda_temp,-1,transition_lambda2)
    transition_lambda = np.sqrt(transition_lambda_temp[0:len(transition_lambda_temp)-1] * transition_lambda_temp[1:])
    
    # Define disk wavelength array
    dlnlambda = 1./settings.specrdisk
    lnlambdamax = np.log(settings.lambdamax*1.01) # we want lambda_disk to be slightly wider than lam for later
    lnlambdamin = np.log(settings.lambdamin*0.99)
    nlambda_disk = int(np.ceil((lnlambdamax-lnlambdamin)/dlnlambda))
    lambda_disk = np.linspace(lnlambdamin,lnlambdamax,nlambda_disk)
    lambda_disk = np.exp(lambda_disk)
    
    # Constants
    nstars = len(stars)
    ncomponents = len(disks[0])
    nplanets = len(planets[0])
    
    for i in range(0,len(stars)):
        
        # ----- SOME BASIC PARAMETERS OF THE SYSTEM -----
        # First, some parameters specific to the star/system we will use later
        s = stars.iloc[i]
        
        pixscale_AU = settings.pixscale_mas/1000. * s['dist']
        dstarAU = s['dist'] * 206265. # distance to star in AU
        GM = grav_const * s['mass']
        cosPA = np.cos(s['PA'] * np.pi/180.)
        sinPA = np.sin(s['PA'] * np.pi/180.)
        cosinc = np.cos(s['I'] * np.pi/180.)
        sininc = np.sin(s['I'] * np.pi/180.)
        
        print('i = {0:d}'.format(i))
        print('starID = {0:d}'.format(s['ID']))
        
        hiptag = ''
        if s['HIP'] > 0: hiptag = s['HIP']
        if settings.output_dir[0] == '/': tempdir = settings.output_dir
        else: tempdir = './' + settings.output_dir # if output_dir is relative or absolute
        if 'TYC2' in stars.columns:
            fits_filename = tempdir + str(s['ID']) + '-HIP_' + str(hiptag) + '-TYC_' + str(s['TYC2']) + '-mv_{0:4.2f}-L_{1:4.2f}-d_{2:4.2f}-Teff_{3:4.2f}.fits'.format(s['Vmag'],s['Lstar'],s['dist'],s['Teff'])
        else:
            fits_filename = tempdir + str(s['ID']) + '-HIP_' + str(hiptag) + '-mv_{0:4.2f}-L_{1:4.2f}-d_{2:4.2f}-Teff_{3:4.2f}.fits'.format(s['Vmag'],s['Lstar'],s['dist'],s['Teff'])
        
        if path.exists(fits_filename): continue # if the output file already exists, we skip this star
        
        # ----- START OF DISK IMAGING ----
        # Now we image the disk
        print('Creating disk image...')
        disk_contrast_image = np.zeros((settings.npix,settings.npix,nlambda_disk))
        
        # Determine range of grain sizes for this star
        settings.minsize = settings.rdust_blowout
        settings.maxsize = min(100*np.max(lam),np.max(master_rdust_boundaries)) # 100x wavelength if possible
        j = np.array([jj for jj in range(0,len(master_rdust)) if settings.minsize < master_rdust[jj] and master_rdust[jj] < settings.maxsize])
        if len(j)==0:
            print('ERROR: no dust grain sizes meet your criteria.')
            break
        nsizes = len(j)
        rdust = master_rdust[j]
        drdust = master_drdust[j]
        rdust_boundaries = master_rdust_boundaries[np.min(j):np.max(j)+2] # one unit longer to include lower and upper boundaries
        
        # We need Qsca vs lambda for all grain sizes
        # This has been previously calculated
        # So we call a routine that loads this from a grid
        composition = compositions[i] # currently we assume the same composition for all components
        
        Qsca = load_lqsca(lqq_dir, composition, rdust, rdust_boundaries, lambda_disk)
        # Qsca is now nsizes x nlambda_disk array
        
        # Now distribute a grid of points and calculate disk flux at all points
        if not settings.diskoff: disk_contrast_image = distribute_diskpoints(s, disks[i], rdust, drdust, Qsca)
        print('...done')
        # ----- END OF DISK IMAGING -----
        
        
        # ----- START OF PLANET CALCULATIONS -----
        # Now the planets
        print('Planet flux/integration...')
        p = planets[i]
        j = np.where(p[:,pllabel.index('R')] > 0)
        tempnplanets = len(j[0])+1  # add 1 to include the star
        planet_data = np.zeros((tempnplanets,16+nlambda,ntimes)) # this will hold all star + planet data vs time & wavelength
        
        # if it's just the star, there's no integration to do
        # just fill in planet_data with the time and flux
        hires, fstarhires = get_stellar_flux(s,lam,path=exovistapath)
        nhires = len(hires)
        pnu = (2.998e14)/hires                       # convert wavelength in microns to frequency in Hz
        pnu = np.append(pnu,pnu[-1]*pnu[-1]/pnu[-2]) # add on an extra one at the end for interpolation
        pdnu = np.abs(pnu[1:]-pnu[0:len(pnu)-1])     # will use this for integration
        fstar = np.zeros((ntimes,nlambda))
        fstarline = np.zeros(nlambda)
        
        if tempnplanets > 1: # there are some valid planets, do the integration stuff
            tempp = p[j]
            Mplanet = tempp[:,pllabel.index('M')]/332948.0 # mass of planets in solar masses
            Mvector = np.insert(Mplanet,0,s['mass'])  # all in solar masses

            Rplanet = tempp[:,pllabel.index('R')]
            starrad = np.sqrt(s['Lstar'] / (s['Teff']/5778.)**4) * 109.2
            Rvector = np.insert(Rplanet,0,starrad)
            Rvector *= 6.371e8/1.496e13
            
            GMplanet = np.array(grav_const * Mplanet)
            GMvector = np.array(grav_const * Mvector)
            GMplanetpstar = GMplanet + grav_const * s['mass']

            kepcoords = []            
            for cd in keplist:
                if cd in ['a','e']: kepcoords.append(tempp[:,pllabel.index(cd)])
                else: kepcoords.append(tempp[:,pllabel.index(cd)]*np.pi/180.)
                
            # Convert orbital elements to heliocentric coordinates
            cart = coord.cartesian(GMplanetpstar, kepcoords)
            
            # Append star's heliocentric coordinates to the front
            x0 = np.insert(np.array(cart[0]),0,0)
            y0 = np.insert(np.array(cart[1]),0,0)
            z0 = np.insert(np.array(cart[2]),0,0)
            vx0 = np.insert(np.array(cart[3]),0,0)
            vy0 = np.insert(np.array(cart[4]),0,0)
            vz0 = np.insert(np.array(cart[5]),0,0)
            # Center of mass coordinates:
            cmx = np.sum(x0*Mvector)/np.sum(Mvector)
            cmy = np.sum(y0*Mvector)/np.sum(Mvector)
            cmz = np.sum(z0*Mvector)/np.sum(Mvector)
            cmvx = np.sum(vx0*Mvector)/np.sum(Mvector)
            cmvy = np.sum(vy0*Mvector)/np.sum(Mvector)
            cmvz = np.sum(vz0*Mvector)/np.sum(Mvector)
        
            # Convert to barycentric coordinates before integrating
            x0 -= cmx
            y0 -= cmy
            z0 -= cmz
            vx0 -= cmvx
            vy0 -= cmvy
            vz0 -= cmvz
            # x,y,z,vx,vy,vz are now vectors w/ star + planet's baryc. coords
            # note that they are relative to the system's midplane
            
            # Later we need the latitude of the planet pointed
            # toward the observer. To do this, we start with unit vectors
            # in the z direction, and rotate it using the same conventions
            # as cartesian.pro. Note we assume zero obliquity.
            # NOTE: THIS SECTION IS USED ONLY FOR LAT-LON RESOLVED ALBEDOS.
            # IT MAY ALSO BE MISSING AN ARGPERI ROTATION.
            x_uv0 = 0.       # a unit vector in the z direction
            y_uv0 = 0.
            z_uv0 = 1.
            # rotate by inclination about x axis
            cospi = np.cos(tempp[:,pllabel.index('i')]*np.pi/180.)
            sinpi = np.sin(tempp[:,pllabel.index('i')]*np.pi/180.)
            x_uv = x_uv0
            y_uv = y_uv0 * cospi - z_uv0 * sinpi
            z_uv = y_uv0 * sinpi + z_uv0 * cospi
            # rotate by longitude of node about z axis
            cospnode = np.cos(tempp[:,pllabel.index('longnode')]*np.pi/180.)
            sinpnode = np.sin(tempp[:,pllabel.index('longnode')]*np.pi/180.)
            x_uv0 = x_uv * cospnode - y_uv * sinpnode
            y_uv0 = x_uv * sinpnode + y_uv * cospnode
            z_uv0 = z_uv
            # Now we have our unit vectors if the system midplane
            # were in the x-y plane. But the system is inclined.
            # The latitude pointed toward the observer does not change
            # as the planet orbits. So we can go ahead and calculate
            # this now.
            # rotate unit vectors about x-axis by system inclination
            x_uv = x_uv0
            y_uv = y_uv0 * cosinc - z_uv0 * sininc
            z_uv = y_uv0 * sininc + z_uv0 * cosinc
            # Finally, rotate unit vectors about z axis by system PA
            x_uv0 = x_uv * cosPA - y_uv * sinPA    # AU
            y_uv0 = x_uv * sinPA + y_uv * cosPA    # AU
            z_uv0 = z_uv                           # AU
            # xyz_uv0 now is the orbital axis of each planet projected
            # onto the sky. We'll use this again later.
            # The dot product between our rotated unit vector and
            # a unit vector pointed toward the observer now gives us
            # the latitude.
            platitude = np.arccos(z_uv0)/np.pi*180. # this sets latitude = 0 when pole-on...
            platitude = 90.-platitude # now we have -90 to +90, with +90 = North

            # ----- INTEGRATE -----
            # Step through all time steps to do this.
            # Make sure variables are expected type for call to C below
            curr_time = 0.
            desired_time = 0.

            cart0 = [x0,y0,z0,vx0,vy0,vz0]
            print('Integrating over time...')

            tlistfull = []
            translistfull = []
            for ip in range(0,tempnplanets-1): translistfull.append([])
                        
            for it in range(0,ntimes):
                if it%16==0: print(it,ntimes)
                # Integrate forward in time
                desired_time = it * settings.dt # dt is in years
                cart_end,tlist,transmaster = nbody.nbody(cart0,GMvector,Rvector,s['I'],curr_time,desired_time)
                
                x0 = cart_end[0]
                y0 = cart_end[1]
                z0 = cart_end[2]
                vx0 = cart_end[3]
                vy0 = cart_end[4]
                vz0 = cart_end[5]
                
                tlistfull = tlistfull + tlist
                for ip in range(0,tempnplanets-1): translistfull[ip] = translistfull[ip] + transmaster[ip]
                
                curr_time = desired_time # update the time
                
                # output coordinates are relative to the system plane
                # Determine the new orbital parameters of the planets
                tempx = x0[1:tempnplanets] - x0[0]  # heliocentric coordinates of just the planets, not the star
                tempy = y0[1:tempnplanets] - y0[0]
                tempz = z0[1:tempnplanets] - z0[0]
                tempvx = vx0[1:tempnplanets] - vx0[0]
                tempvy = vy0[1:tempnplanets] - vy0[0]
                tempvz = vz0[1:tempnplanets] - vz0[0]

                cartcoords = [tempx, tempy, tempz, tempvx, tempvy, tempvz]
                orb = coord.keplerian(GMplanetpstar, cartcoords)
                
                tempa = orb[0]
                tempe = orb[1]
                tempi = orb[2]/np.pi*180. # convert to degrees
                templongnode = orb[3]/np.pi*180.
                tempargperi = orb[4]/np.pi*180.
                tempmeananom = orb[5]/np.pi*180.
                
                # Convert output coordinates (relative to system midplane) to on-sky coordinates
                # rotate about x-axis by system inclination
                xtemp = x0
                ytemp = y0 * cosinc - z0 * sininc
                ztemp = y0 * sininc + z0 * cosinc
                vxtemp = vx0
                vytemp = vy0 * cosinc - vz0 * sininc
                vztemp = vy0 * sininc + vz0 * cosinc
                # rotate about z-axis by PA
                x = xtemp * cosPA - ytemp * sinPA    # AU
                y = xtemp * sinPA + ytemp * cosPA    # AU
                z = ztemp                            # AU
                vx = vxtemp * cosPA - vytemp * sinPA # AU
                vy = vxtemp * sinPA + vytemp * cosPA # AU
                vz = vztemp
                
                # Record the coordinates and orbital parameters for each body
                for ip in range(0,tempnplanets):
                    planet_data[ip,0,it] = curr_time
                    planet_data[ip,9:15,it] = [x[ip],y[ip],z[ip],vx[ip],vy[ip],vz[ip]] # barycentric coordinates
                    
                    if ip==0: # if it's the star...                
                        planet_data[ip,3:9,it] = 0. # don't record a heliocentric orbit 
                    else: # if it's a planet... 
                        planet_data[ip,3:9,it] = [tempa[ip-1],tempe[ip-1],tempi[ip-1],templongnode[ip-1],tempargperi[ip-1],tempmeananom[ip-1]] # record the heliocentric orbit
            #----- END OF INTEGRATION -----

            #----- COMPILE LIST OF TRANSITS AND ECLIPSES -----
            tlistfinal = []
            plistfinal = []
            translist0 = []
            translist1 = []
            for ip in range(1,tempnplanets):
                if translistfull[ip-1][0] != 0:
                    tlistfinal.append(tlistfull[0]*365.25)
                    plistfinal.append(ip)
                    translist0.append(translistfull[ip-1][0])
                    translist1.append(translistfull[ip-1][0])
            
            for it in range(1,len(tlistfull)):
                for ip in range(1,tempnplanets):
                    if translistfull[ip-1][it] != translistfull[ip-1][it-1]:
                        tlistfinal.append(tlistfull[it]*365.25)
                        plistfinal.append(ip)
                        translist0.append(translistfull[ip-1][it-1])
                        translist1.append(translistfull[ip-1][it])
                        
            for ip in range(1,tempnplanets):
                if translistfull[ip-1][-1] != 0:
                    tlistfinal.append(tlistfull[-1]*365.25)
                    plistfinal.append(ip)
                    translist0.append(translistfull[ip-1][-1])
                    translist1.append(translistfull[ip-1][-1])
                    
            if len(tlistfinal)==0:
                tlistfinal = [0.]
                plistfinal = [0]
                translist0 = [0]
                translist1 = [0]
            transits = np.array([np.array(tlistfinal),np.array(plistfinal),np.array(translist0),np.array(translist1)])
            
            #----- CALCULATE FLUX -----
            '''
            The integration above was done by looping over
            time, and then planets within that. That's slow
            For calculating planet flux b/c we have to load
            a different planet albedo for each planet. So here
            we loop over planets and then treat time inside of
            that loop.
            Stellar flux vs time has already been recorded
            '''
            for ip in range(1,tempnplanets):
                # Calculate phase angle
                xstar = planet_data[0,9,:]  # barycentric coordinates vs time
                ystar = planet_data[0,10,:]
                zstar = planet_data[0,11,:]              
                x = planet_data[ip,9,:]     # barycentric coordinates vs time
                y = planet_data[ip,10,:]
                z = planet_data[ip,11,:]
                x_helio = x - xstar         # heliocentric coordinates vs. time
                y_helio = y - ystar
                z_helio = z - zstar
                xx_helio = x_helio*x_helio
                yy_helio = y_helio*y_helio
                zz_helio = z_helio*z_helio
                r_helio = np.sqrt(xx_helio+yy_helio+zz_helio)
                # we assume telescope always points exactly at star
                dz = dstarAU - z_helio                                                          # z-component of vector joining point to telescope
                l = np.sqrt(xx_helio + yy_helio + dz*dz)                                        # distance from point to telescope in AU
                cosbeta = (-z_helio*(dstarAU-z_helio) + xx_helio + yy_helio) / (r_helio * l)    # cos(phase angle) calculated from a dot product
                beta = np.arccos(cosbeta)                                                       # phase angle in radians
                xpix = x_helio / pixscale_AU + settings.npix/2. # convert to pixel coordinates
                ypix = y_helio / pixscale_AU + settings.npix/2.
                
                planet_data[ip,1,:] = xpix # pixel coordinates
                planet_data[ip,2,:] = ypix # pixel coordinates
                
                planet_data[ip,15,:] = beta*180./np.pi # phase angle
                
                # Read the albedo file
                # Upon return, phi could be undefined, phase angle, or differential longitude
                # Latitude could be undefined or latitude
                albedo_file = 'geometric_albedo_files/' + albedos[i][ip-1] + '.txt'
                plambda0, phi, lat, data0 = read_albedo_file(albedo_file)
                ll = 0
                for jj in range(0,len(plambda0)):
                    if plambda0[jj]==ll: print(jj,plambda0[jj])
                    ll = plambda0[jj]
                    
                # Get planet_radius/circumstellar_distance and reformat
                Rpor = (tempp[ip-1,pllabel.index('R')] * 4.2635e-5) / r_helio # planet radius in AU divided by circumstellar distance in AU
                Rpor2 = Rpor * Rpor

                # output array
                fp = np.zeros((ntimes,nhires))
                
                # EXOPLANET INPUT MODEL CASE 1: We need to apply Lambertian phase function
                if len(data0.shape)==1: # in this case, data0 is geometric albedo
                    #print('CASE 1: albedo vs. wavelength only.')
                    
                    # Evaluate the phase function at all times
                    phasefunc = lambertian(beta)

                    # interpolate data0 to the hires array
                    # extrapolate because it was giving errors precisely at the endpoints.
                    f = interp1d(plambda0,data0,kind='cubic',bounds_error=False,fill_value='extrapolate')
                    data = f(hires)
                    
                    # calculate planet flux (modulo some unimportant factors)
                    # data is geometric albedo
                    for ii in range(0,ntimes): fp[ii] = data * Rpor2[ii] * phasefunc[ii] * fstarhires
                    
                
                # EXOPLANET INPUT MODEL CASE 2: The model is phase-resolved
                if len(data0.shape)==2: # in this case, data0 is geometric albedo times the phase function (gI)
                    #print('CASE 2: albedo vs. wavelength and phase.')
                    
                    # First, get the indices of all desired phases:
                    f = interp1d(phi,np.arange(len(phi)),kind='cubic')
                    phase_indices = f(beta*180./np.pi)
                    phase_indices = np.maximum(phase_indices, 0)
                    phase_indices = np.minimum(phase_indices, len(phi)-1)
                    phase_arg = np.argsort(np.argsort(phase_indices))
                    
                    # Now interpolate to the desired phases and wavelengths
                    # result is phase x plambda array of gI values                    
                    f2 = interp2d(np.arange(len(data0[0])),plambda0,data0,kind='cubic')
                    data = f2(phase_indices,hires)
                    data = np.maximum(data,0)

                    for di in range(0,len(data)): data[di] = data[di][phase_arg]
                    data = np.transpose(data)
                    
                    # calculate planet flux (modulo some unimportant factors)
                    # data is geometric albedo times phase function
                    for ii in range(0,ntimes): fp[ii] = data[ii] * Rpor2[ii] * fstarhires
                    
        
                # EXOPLANET INPUT MODEL CASE 3: The model is resolved by 
                # differential longitude and latitude
                if len(data0.shape)==3: # in this case, data0 is geometric albedo times the phase function (gI)
                    #print('CASE 3: albedo vs. wavelength, longitude, and latitude.')
                    
                    # Below we calculate dlong, the difference between the
                    # planet's longitude pointing toward the star and the
                    # planet's longitude pointing toward the observer.
                    # A derivation of this appears in the exoVista STG meeting
                    # notes from 09 Jan 2022. This implicitly assumes zero obliquity
                    # and the inclination of the planet's orbit doesn't
                    # change over time.
                    
                    d_x = 0.
                    d_y = 0.
                    d_z = dstarAU
                    ddotuv = d_x * x_uv0[ip-1] + d_y * y_uv0[ip-1] + d_z * z_uv0[ip-1]
                    ddotuv_x = ddotuv * x_uv0[ip-1]
                    ddotuv_y = ddotuv * y_uv0[ip-1]
                    ddotuv_z = ddotuv * z_uv0[ip-1]
                    s_x = d_x - ddotuv_x
                    s_y = d_y - ddotuv_y
                    s_z = d_z - ddotuv_z
                    p_x = -x_helio
                    p_y = -y_helio
                    p_z = -z_helio
                    p_mag = r_helio
                    o_x = p_x + s_x
                    o_y = p_y + s_y
                    o_z = p_z + s_z
                    o_mag = np.sqrt(o_x * o_x + o_y * o_y + o_z * o_z)
                    dlong = np.arccos((p_x*o_x + p_y*o_y + p_z*o_z)/(p_mag*o_mag)) # differential longitude (observer's planet longitude - stellar illumination longitude)
                    dlong *= (180./np.pi)
                    
                    # First, get the latitude index
                    f = interp1d(lat,np.arange(len(lat)),kind='cubic')
                    lat_indices = f(platitude[ip-1])
                    lat_indices = max(lat_indices,0)
                    lat_indices = min(lat_indices,len(lat)-1) # this is a single index for 1 planet
                    
                    # Now get the indices of all desired longitudes:
                    f2 = interp1d(phi,np.arange(len(phi)),kind='cubic')
                    lon_indices = f2(dlong)
                    lon_indices = np.maximum(lon_indices,0)
                    lon_indices = np.minimum(lon_indices,len(phi)-1)
                    lon_arg = np.argsort(np.argsort(lon_indices))

                    # I don't know how I got this bit working. Change it at your own risk.
                    dataslice = []                    
                    wlist = np.arange(len(data0))      # placeholder index to keep the wavelength dimension while flattening the latitude
                    llist = np.arange(len(data0[0,0]))
                    
                    for ii in range(0,len(data0[0,0])):
                        dslice = data0[:,ii,:]
                        fi = interp2d(llist,wlist,dslice,kind='cubic')
                        tempgrid = np.transpose(fi([lat_indices],wlist)[:,0])
                        dataslice.append(tempgrid)

                    # This array is (should be) longitude x wavelength.
                    dataslice = np.transpose(np.array(dataslice))
                    f3 = interp2d(np.arange(len(dataslice[0])),plambda0,dataslice,kind='cubic')
                    
                    data = f3(lon_indices,hires)
                    data = np.maximum(data,0)
                    for di in range(0,len(data)): data[di] = data[di][lon_arg]
                    data = np.transpose(data)
                    
                    # calculate planet flux (modulo some unimportant factors)
                    # data is geometric albedo times phase function
                    for ii in range(0,ntimes): fp[ii] = data[ii] * Rpor2[ii] * fstarhires
                    
                # Now that we have fp, we can bin it by wavelength
                # and divide by the observer's fstar to get contrast
                integratedFp = np.zeros((nlambda,ntimes))
                
                nl = 0
                jlist = []
                for j in range(0,len(hires)):
                    if hires[j] >= transition_lambda[nl]:
                        nl+=1
                        jlist.append(j)
                        if nl==len(transition_lambda): break
                
                for ilambda in range(0,len(lam)):
                    minj = jlist[ilambda]
                    maxj = jlist[ilambda+1]
                    # Bin the stellar spectrum on the first pass
                    for jl in range(0,ntimes):
                        planet_data[0,0,jl] = jl * settings.dt # stellar coordinates
                        planet_data[0,1:3,jl] = [settings.npix/2.,settings.npix/2.]
                        integratedFp[ilambda,jl] = np.sum(fp[jl,minj:maxj+1] * pdnu[minj:maxj+1]) / np.sum(pdnu[minj:maxj+1]) # integrate over spectral channel and divide by total dnu
                    if ip==1:
                        fstarline[ilambda] = np.sum(fstarhires[minj:maxj+1] * pdnu[minj:maxj+1]) / np.sum(pdnu[minj:maxj+1])
                if ip==1:
                    fstar[:] = fstarline
                    fstar = np.transpose(fstar)
                    planet_data[0,16:16+nlambda,:] = fstar
                cp = integratedFp / fstar             # convert to contrast by dividing by recorded fstar
                planet_data[ip,16:16+nlambda,:] = cp  # contrast
            
        # ----- END OF FLUX CALCULATION -----
                
        print('generate_scene() done for Star {0:d}.'.format(s['ID']))
        # ----- END OF PLANET FLUX/INTEGRATION -----
        
        
        # ----- CREATE A FITS FILE FOR EACH STAR -----

        # PRIMARY HDU: WAVELENGTHS FOR STAR AND PLANET SPECTRUM
        hdr0 = fits.Header()

        j = np.where(p[:,pllabel.index('R')]>0)[0]
        tempnp = len(j)

        hdr0['DATE'] = str(datetime.now())
        hdr0.comments['DATE'] = 'Date and time created'
        hdr0['VERSION'] = 2.2
        hdr0.comments['VERSION'] = 'Version of code used; used for post-processing scripts.'
        hdr0['N_EXT'] = 3+tempnp
        hdr0.comments['N_EXT'] = 'Last extension' # need to add this to the main header to enable extensions
        hdr0['SPECRES'] = settings.specres
        hdr0.comments['SPECRES'] = 'Spectral resolution for star and planet'
        hdr0['LAMMIN'] = settings.lambdamin
        hdr0.comments['LAMMIN'] = 'Minimum wavelength (microns)'
        hdr0['LAMMAX'] = settings.lambdamax
        hdr0.comments['LAMMAX'] = 'Maximum wavelength (microns)'
        hdr0['COMMENT'] = 'Stellar wavelength vector'

        hdu0 = fits.PrimaryHDU(lam, header=hdr0)

        # HDU 1: WAVELENGTHS FOR DISK SPECTRUM
        hdr1 = fits.Header()
        hdr1['SPECRES'] = settings.specrdisk
        hdr1.comments['SPECRES'] = 'Spectral resolution of disk wavelengths'
        hdr1['LAMMIN'] = settings.lambdamin
        hdr1.comments['LAMMIN'] = 'Minimum wavelength (microns)'
        hdr1['LAMMAX'] = settings.lambdamax
        hdr1.comments['LAMMAX'] = 'Maximum wavelength (microns)'
        hdr1['COMMENT'] = 'Disk wavelength vector'
        
        hdu1 = fits.ImageHDU(lambda_disk, header=hdr1)

        # HDU 2: DISK CONTRAST CUBE
        hdr2 = fits.Header()
        hdr2['NCOMP'] = ncomponents
        hdr2.comments['NCOMP'] = 'Number of disk components'
        hdr2['PXSCLMAS'] = settings.pixscale_mas
        hdr2.comments['PXSCLMAS'] = 'Pixel scale (mas)'

        for icomp in range(0,ncomponents):
            d = disks[i,icomp]
            for ih in range(1,len(dlabel)):
                if dlabel[ih]=='longnode': tag = 'LNGND' + '-' + str(icomp)
                else: tag = dlabel[ih].upper() + '-' + str(icomp)
                hdr2[tag] = d[ih]
                hdr2.comments[tag] = 'Disk comp. {0:d} {1:s}'.format(icomp,dcomments[dlabel[ih]])
                
        hdr2['MINSIZE'] = np.min(rdust_boundaries)
        hdr2.comments['MINSIZE'] = 'Minimum grain size (microns)'
        hdr2['MAXSIZE'] = np.max(rdust_boundaries)
        hdr2.comments['MAXSIZE'] = 'Minimum grain size (microns)'
        hdr2['COMPOSIT'] = composition
        hdr2.comments['COMPOSIT'] = 'Dust composition'
        hdr2['COMMENT'] = 'Disk contrast cube'

        # arrays must be transposed for backward compatibility with IDL outputs
        hdu2 = fits.ImageHDU(disk_contrast_image.T, header=hdr2)

        # HDU T: LIST OF TRANSITS AND ECLIPSES
        hdrt = fits.Header()
        hdrt['BASELINE'] = settings.timemax*365.25
        hdrt.comments['BASELINE'] = 'Length of baseline (days)'
        
        hdut = fits.ImageHDU(transits, header=hdrt)
        
        # HDU 3: STAR PROPERTIES
        hdr3 = fits.Header()
        ic = 0
        for ih in ['PA','I','ID','HIP','TYC2','dist','M_V','Vmag','Bmag','Umag','Rmag','Imag','Jmag','Hmag','Kmag','Type','Lstar','logg','Teff','angdiam','mass','rstar']:
            if ih not in s.index: hdr3[ih] = 'Invalid tag'
            elif pd.isnull(s[ih]): hdr3[ih] = 'NaN'
            else: hdr3[ih] = s[ih]
            if ih in scomments: hdr3.comments[ih] = scomments[ih]
        hdr3['PXSCLMAS'] = settings.pixscale_mas
        hdr3.comments['PXSCLMAS'] = 'Pixel scale (mas)'
        hdr3['COMMENT'] = 'Star data array'
        
        hdu3 = fits.ImageHDU(planet_data[0].T, header=hdr3)
        
        hdul = fits.HDUList([hdu0, hdu1, hdu2, hdut, hdu3])

        # HDRN: PLANET PROPERTIES
        for ip in range(1,len(planet_data)):
            hdrn = fits.Header()
            
            for ih in range(0,len(pllabel)):
                tag = pllabel[ih].upper()
                hdrn[tag] = p[ip-1,ih]
                hdrn.comments[tag] = pcomments[pllabel[ih]]
            hdrn['ALBEDO_F'] = albedos[i][ip-1]
            hdrn.comments['ALBEDO_F'] = 'Geometric albedo file'
            hdrn['PXSCLMAS'] = settings.pixscale_mas
            hdrn.comments['PXSCLMAS'] = 'Pixel scale (mas)'
            hdrn['COMMENT'] = 'Planet data array'
            
            hdun = fits.ImageHDU(planet_data[ip].T, header=hdrn)

            hdul.append(hdun)
            
        hdul.writeto(fits_filename)

    print('generate_scene done.')
    return


def distribute_diskpoints(s, disk, rdust, drdust, Qsca, xcen=0, ycen=0, xwidth=settings.npix, ywidth=settings.npix):
    # Tolerance is set to give noise level at 1/3rd of faintest
    # detectable point source (dmag=26.5) at V band for a 4 m
    # telescope.
    # NOTE: xcen and ycen are the offset of the sample area along and perpendicular to the disk's axis of inclination, respectively
    
    tol = 0.05/np.max(disk[:,dlabel.index('nzodis')])
    tol = max(tol,0.0005) # if nzodis>100, we're not going to be looking for dmag=26.5 objects anyway
    tol = min(tol,0.05)   # at least hit 5% precision even if nzodis is < 1
    tol0 = tol

    # NOTE: use pixscale_mas = 4.0 used for calibrating the exozodi.
    
    pixscale_arcsec = settings.pixscale_mas / 1000.
    auperpix = pixscale_arcsec * s['dist']
    cosinc = np.cos(-s['I']*np.pi/180) # we want to rotate in opposite direction of inclination
    sininc = np.sin(-s['I']*np.pi/180)
    cosPA = np.cos(-s['PA']*np.pi/180) # we want to rotate in opposite direction of position angle
    sinPA = np.sin(-s['PA']*np.pi/180)
    distAU = s['dist'] * 206265.0 # distance to star in AU
    nlambda = len(Qsca[0])
    
    # set bounds of output data cube
    rcen = np.sqrt(xcen**2 + ycen**2)
    xcenrot = rcen * np.cos(s['PA']*np.pi/180)
    ycenrot = rcen * np.sin(s['PA']*np.pi/180)

    xcenpix = int(xcenrot/auperpix + settings.npix/2)
    ycenpix = int(ycenrot/auperpix + settings.npix/2)
    xstart = xcenpix - int(np.ceil(xwidth/2))
    xend = xstart + xwidth
    ystart = ycenpix - int(np.ceil(ywidth/2))
    yend = ystart + ywidth
    
    masterimg = np.zeros((xwidth,ywidth,nlambda+1))
    precision = np.zeros((xwidth,ywidth))
    
    # go far enough that expected SB drops to tol of peak
    
    zmax_radial = np.max(disk[:,dlabel.index('r')])*tol**(-2./7.) # based on n ~ r^-1.5 and 1/r^2 illumination factor
    rmax = np.sqrt(2.)*(settings.npix/2.0)*auperpix
    zmax_height = np.sqrt(-2 * (np.max(disk[:,dlabel.index('hor')]) * rmax)**2 * np.log(tol)) # exponential variation with z

    # We use zmax_height (z limit based on tol and opening angle of disk) 
    # to estimate inc0, which is the inclination at which we start looking
    # through the radial extent of the disk.

    zmax = zmax_radial # *abs(sininc) > zmax_height*abs(cosinc)
    
    r,x0,y0 = rgen(settings.npix,settings.npix)
    x0 -= 0.5 # shift to 00LL convention
    y0 -= 0.5 # shift to 00LL convention
    x0 *= auperpix # convert to AUs
    y0 *= auperpix # convert to AUs
    
    r_arcsec = r * pixscale_arcsec
    
    # We do this on a pixel-by-pixel basis
    # because otherwise it takes too much RAM
    dx0 = auperpix               # in AUs
    dy0 = auperpix               # in AUs
    prevnxy = 0
    prevnz = 0
          
    nsubarraypix = 1
    bksp = ''
    
    if xstart<0 or xend>settings.npix or ystart<0 or yend>settings.npix:
        print('Error: disk sample area out of aperture.')
        exit()
    
    diskimage = wrapImage.PyImage()
    diskimage.SetupImage(s['rstar'], s['Teff'], settings.rdust_blowout, settings.tsublimate, disk, rdust, drdust, Qsca)
    
    for ix in range(xstart,xend):
        maxnxy_used = 0
        maxnz_used = 0
        for iy in range(ystart,yend):
            # reduce tolerance on pixels interior to iwa
            # this significantly improves run time
            # and memory usage, as most time is spent
            # on the central pixels
            if r_arcsec[ix,iy] < settings.iwa: tol = max(settings.iwa_tol, tol0)
            else: tol = tol0
            
            nxy = 1     # (prevnxy-2) > long(1)
            nz = 128    # (prevnz/4) > long(100)
            maxnxy = 32 # the highest nxy can ever get
            maxnz = 32768     # the higest nz can ever get
            prevnxy = 0
            prevnz = 0
            iteration = 0
            finished_this_pixel = 0
            limit_flag = 0
            
            while finished_this_pixel != 2:
                # update the dvolume values based on # of sub-points
                dx = dx0/nxy
                dy = dy0/nxy

                # create 3D array
                subpixelrows = np.arange(0,nsubarraypix)*nxy

                xref = (np.arange(0,nxy*nsubarraypix)/nxy + 0.5/nxy) * auperpix
                
                x = np.zeros((nxy*nsubarraypix,nxy*nsubarraypix,nz))
                y = np.zeros((nxy*nsubarraypix,nxy*nsubarraypix,nz))
                
                for i in range(0,len(x)):
                    x[i,:,:] = xref[i] + x0[ix]
                    y[:,i,:] = xref[i] + y0[iy]
                    
                # Linear dz takes too long          
                # dz ~ z^zpower method is faster and provides better resolution
                # near the disk midplane typically
                zpower = 2.0
                dz = 2./nz
                zvals_vector = np.arange(0,nz+1)*dz - 1 # array of endpoints of volumes
                zvals_vector *= zmax**(1./(zpower+1))
                zvals_vector = zvals_vector * np.abs(zvals_vector)**zpower
                dzvals_vector = np.abs(zvals_vector[1:]-zvals_vector[0:nz]) 
                zvals_vector = (zvals_vector[1:]+zvals_vector[0:nz])/2. # midpoints
                z = np.zeros((nxy*nsubarraypix,nxy*nsubarraypix,nz))
                dz = np.zeros((nxy*nsubarraypix,nxy*nsubarraypix,nz))
                
                for i in range(0,len(z)):
                    for j in range(0,len(z[0])):
                        z[i,j] = zvals_vector
                        dz[i,j] = dzvals_vector
                dv = dx*dy*dz
                
                # now we have nxy x nxy lines of sight in our pixel
                # resolved into nz steps along the line of sight
                
                # calculate scattering angle
                r = np.sqrt(x*x+y*y+z*z)
                l = np.sqrt(x*x+y*y+distAU*distAU)
                cosscattang = ((-x*x - y*y + z*distAU) / (r * l))
           
                '''
                Transform sky-plane xyz coords to disk-plane coords.
                Normally we do the reverse. So here we have to           
                perform the rotations in the opposite order
                and opposite direction as normal. Note that 
                sininc/cosinc and sinPA/cosPA have already been 
                calculated w/ the opposite sign included.
                rotate about z-axis by -PA
                '''
                xtemp = x * cosPA - y * sinPA  # AU
                
                ytemp = x * sinPA + y * cosPA  # AU
                ztemp = z                      # AU
                # using the opposite signs
                x = xtemp
                y = ytemp * cosinc - ztemp * sininc
                z = ytemp * sininc + ztemp * cosinc

                '''
                Approximate SB variation at these sub-points
                we ignore scattering phase function effects
                img = disk_imager_simplified(s, x, y, z, dv, limit_flag)          
                disk_imager_simplified is much faster for optimization,
                but in the end you have the run the full calculation
                at the highest resolution anyway, so there isn't much
                savings.
                disk_imager takes our grid of points, calculates the disk
                density, illuminates with starlight, and calculates a disk
                image cube for the subarray of pixels
                '''
                sdata = []
                img = diskimage.disk_imager(x, y, z, r, dv, cosscattang)
                img = np.array(img)
                
                # Now integrate all sub-pixels
                tempimg = img[:,subpixelrows,:]
                for inxy in range(1,nxy): tempimg += img[:,subpixelrows+inxy,:]
                img = tempimg[subpixelrows,:,:]
                for inxy in range(1,nxy): img += tempimg[subpixelrows+inxy,:,:]
                
                if iteration > 0:
                    diff = np.abs(img-previmg)/previmg
                    precision[ix-xstart:ix+nsubarraypix-xstart,iy-ystart:iy+nsubarraypix-ystart] = np.max(diff)
                    if np.max(diff) < tol: finished_this_pixel += 1
                    if finished_this_pixel == 0 and nz >= maxnz: finished_this_pixel += 1
                    if finished_this_pixel == 1 and nxy >= maxnxy: finished_this_pixel += 1
                previmg = img
                prevnxy = nxy
                prevnz = nz
                
                if nz > maxnz_used: maxnz_used = nz
                if nxy > maxnxy_used: maxnxy_used = nxy
                
                # If the tolerance has not been met in the z direction, double
                # the number of points in the z direction
                if finished_this_pixel == 0: nz *= 2
                # If the tolerance has not been met in the xy direction, double
                # the number of points in the xy direction
                if finished_this_pixel == 1: nxy *= 2
                
                iteration += 1
            for jx in range(0,nsubarraypix):
                for jy in range(0,nsubarraypix):
                    masterimg[ix+jx-xstart,iy+jy-ystart,0:nlambda] = img
                    masterimg[ix+jx-xstart,iy+jy-ystart,-1] = precision[ix+jx-xstart,iy+jy-ystart]
        if ix%16==0: print(ix,maxnxy_used,maxnz_used,tol0)
        
    return masterimg

def rgen(numx, numy=0):
    # Generates a 2D distribution of radii.  The origin is assumed
    # to be in the center of the matrix.  Radii are calculated
    # assuming 1 pixel = 1 unit.  The number of pixels can be odd
    # or even.
    
    # numx = # of pixels in x direction
    # numy = # of pixels in y direction
    # x: the x coordinates
    # y: the y coordinates

    if numy==0: numy = numx
    
    x = np.linspace(0,numx-1,numx) - float((numx-1)/2.)
    y = np.linspace(0,numy-1,numy) - float((numy-1)/2.)
    
    r = np.zeros((numx,numy))
    for i in range(0,numx): r[i] = np.sqrt(x[i]*x[i] + y*y)

    return r,x,y


######### WARNING: THIS ROUTINE IS LIKELY TO BE REWRITTEN #########
def get_stellar_flux(s, lam, path='./'):
    # INPUTS
    # s: star structure
    # lam: desired wavelengths (microns)
    
    # OUTPUTS
    # fstar: stellar flux in Jy
    print('Computing stellar spectrum...')
    
    Tstar = s['Teff']
    logg = s['logg']
    metallicity = 0.0
    Rstar = s['rstar']
    dstar = s['dist']
    RstarAU = Rstar * 0.00465047  # radius of star in AU
    dstarAU = s.dist * 206265.    # distance to star in AU
    
    # Load the appropriate Kurucz stellar atmosphere model
    # klambda in units of nm
    # kBnu is in units of erg s^-1 cm^-2 Hz^-1 steradian^-1

    klambda, kBnu = getkurucz(Tstar, logg, metallicity)
    klambda *= 0.001              # convert to microns
    
    # Interpolate to high spectral resolution
    specres = 10000            # internal spectral resolution

    # The transition (bin edge) wavelengths have to be recalculated here
    # to account for the different wavelength points in the Kurucz spectra.
    dlnlambda = 1./specres
    lnlambdamax = np.log(np.max(klambda))
    lnlambdamin = np.log(np.min(klambda))
    nlambda = int(np.ceil((lnlambdamax-lnlambdamin)/dlnlambda))
    interplambda = np.linspace(lnlambdamin,lnlambdamax,nlambda)
    interplambda = np.exp(interplambda)
    transition_lambda0 = np.sqrt(lam[0:len(lam)-1]*lam[1:])
    transition_lambda_temp = np.insert(transition_lambda0,0,transition_lambda0[0]*lam[0]/lam[1])
    transition_lambda = np.insert(transition_lambda_temp,-1,transition_lambda0[-1]*lam[-1]/lam[-2])
    
    interplambda = np.concatenate((interplambda,lam,transition_lambda)) # add these to the array
    # I included lambda above to make sure there was at least 1 
    # wavelength between all transition wavelengths
    interplambda = np.unique(sorted(interplambda))        
    f = interp1d(klambda, kBnu, kind='cubic',bounds_error=False,fill_value='extrapolate')  # spline option was commented out
    interpBnu = f(interplambda)
     
    # Convert to stellar flux in units of erg s^-1 cm^-2 Hz^-1 
    fstar = np.pi * interpBnu * (RstarAU / dstarAU)**2
    
    # Convert to Jy
    fstar *= 1e23
        
    return interplambda, fstar


######### WARNING: THIS ROUTINE IS LIKELY TO BE REWRITTEN #########
def getkurucz(teff, logg, metallicity=0.0):
    '''
    For a given Teff, log(g), and metallicity, retrieves the
    Kurucz & Castelli ATLAST9 stellar atmosphere models.
    Ideally it would interpolate these, but right now it just
    picks the closest match.
    As recommended by the STScI web site discussing synphot,
    we use the models with vturb=2 km/s and Delta(log(tau_ross))=0.125
    Currently we only read in the zero metallicity file, though more files
    could be added to select from in the future...
    
    INPUTS
    Teff: stellar effective temp (K)
    logg: log of surface gravity
    metallicity: currently only zero is allowed
    
    OUTPUTS
    lambda: wavelength vector (nm)
    Bnu: surface brightness (erg s^-1 cm^-2 Hz^-1 steradian^-1)
    
    OPTIONAL OUTPUTS
    Bnucont: continuum surface brightness w/out line abosrptions (erg s^-1 cm^-2 Hz^-1 steradian^-1)
    logg_model: best matching value of logg in the grid
    teff_model: best matching value of teff in the grid
    
    Vectors of logg and teff in the Kurucz model grid...
    '''
    
    gvec = [0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]
    tvec = [3500,3750,4000,4250,4500,4750,5000,5250,5500,5750,6000,6250,6500,6750,7000,7250,7500,7750,8000,8250,8500,8750,9000,9250,9500,9750,10000,10250,10500,10750,11000,11250,11500,11750,12000,12250,12500,12750,13000,14000,15000,16000,17000,18000,19000,20000,21000,22000,23000,24000,25000,26000,27000,28000,29000,30000,31000,32000,33000,34000,35000,36000,37000,38000,39000,40000,41000,42000,43000,44000,45000,46000,47000,48000,49000,50000]
    
    # First, some error handling...
    if teff > np.max(tvec): print('Teff too large.')
    if teff < np.min(tvec): print('Teff too small.')
    if logg > np.max(gvec): print('logg too large.')
    if logg < np.min(gvec): print('logg too small.')
    if metallicity != 0:
        print('Only [+0.0] metallicity currently allowed.')
        return
    
    # Open model grid file
    fin = open('fp00k2odfnew.pck','r')
    for i in range(0,22): fin.readline()

    # Read the wavelengths (in nm)
    lam = []
    line = fin.readline()
    while line.split()[0] != 'TEFF':
        entries = line.split()
        for l in entries: lam.append(float(l))
        line = fin.readline()
    lam = np.array(lam)
    nl = len(lam)
    
    # Must keep logg_model from going out of bounds
    teff_model = teff
    logg_model = logg
    if teff_model >= 6250: logg_model = max(logg_model, 0.5)
    if teff_model >= 7750: logg_model = max(logg_model, 1.0)
    if teff_model >= 8500: logg_model = max(logg_model, 1.5)
    if teff_model >= 9250: logg_model = max(logg_model, 2.0)
    if teff_model >= 12000: logg_model = max(logg_model, 2.5)
    if teff_model >= 20000: logg_model = max(logg_model, 3.0)
    if teff_model >= 27000: logg_model = max(logg_model, 3.5)
    if teff_model >= 32000: logg_model = max(logg_model, 4.0)
    if teff_model >= 40000: logg_model = max(logg_model, 4.5)
    if teff_model >= 50000: logg_model = max(logg_model, 5.0)
  
    # Get the fractional indices of the desired logg and teff
    if logg_model >= gvec[-1]: gindex = len(gvec)-2
    else: gindex = np.where(gvec >= logg_model)[0][0]-1
    if gindex<0: gindex = 0
    if gindex>=len(gvec)-1: gindex = len(gvec)-2
    fg = (logg_model - gvec[gindex]) / (gvec[gindex+1] - gvec[gindex])

    if teff_model >= tvec[-1]: tindex = len(tvec)-2
    else: tindex = np.where(tvec >= teff_model)[0][0]-1
    if tindex<0: tindex = 0
    if tindex>=len(tvec)-1: tindex = len(tvec)-2
    ft = (teff_model - tvec[tindex]) / (tvec[tindex+1] - tvec[tindex])
    
    Bnu00 = np.zeros(len(lam))
    Bnu01 = np.zeros(len(lam))
    Bnu10 = np.zeros(len(lam))
    Bnu11 = np.zeros(len(lam))
    
    while len(line)>0:
        if line.split()[0] == 'TEFF':
            tempt = float(line.split()[1])
            tempg = float(line.split()[3])

            if tempt==tvec[tindex] and tempg==gvec[gindex]:
                jl = 0
                while jl < nl:
                    line = fin.readline()
                    for jc in range(0,int(len(line)/10)):
                        Bnu00[jl] = float(line[jc*10:(jc+1)*10])
                        jl+=1
            
            if tempt==tvec[tindex] and tempg==gvec[gindex+1]:
                jl = 0
                while jl < nl:
                    line = fin.readline()
                    for jc in range(0,int(len(line)/10)):
                        Bnu01[jl] = float(line[jc*10:(jc+1)*10])
                        jl+=1
            
            if tempt==tvec[tindex+1] and tempg==gvec[gindex]:
                jl = 0
                while jl < nl:
                    line = fin.readline()
                    for jc in range(0,int(len(line)/10)):
                        Bnu10[jl] = float(line[jc*10:(jc+1)*10])
                        jl+=1
            
            if tempt==tvec[tindex+1] and tempg==gvec[gindex+1]:
                jl = 0
                while jl < nl:
                    line = fin.readline()
                    for jc in range(0,int(len(line)/10)):
                        Bnu11[jl] = float(line[jc*10:(jc+1)*10])
                        jl+=1
            else:
                try: line = fin.readline()
                except: break
        else:
            try: line = fin.readline()
            except: break
            
    Bnu0 = Bnu00 * (1.-fg) + Bnu01 * fg
    Bnu1 = Bnu10 * (1.-fg) + Bnu11 * fg
    Bnu = Bnu0 * (1.-ft) + Bnu1 * ft
    Bnu = Bnu * 4.
    
    '''
    IMPORTANT NOTE: The Kurucz flux data is in units of erg s^-1 cm^-2 Hz^-1 sr^-1.  This
    is the same units as surface brightness, so we might think that Kurucz is providing Bnu, where
    Fnu = pi * Bnu * (Rstar / d)^2 is the flux Fnu for a star of radius Rstar at distance d.  But in
    fact, for whatever reason, Kurucz is providing the quantity (pi * Bnu) / (4 * pi).  This is the
    surface flux per unit solid angle--this must be important in stellar atmospheres or something.
    Anyway, to calculate Bnu, which is what this routine returns, we multiply Kurucz's numbers by 4.
    If this sounds crazy to you, I note that this agrees with the SYNPHOT Data User's Guide, which
    can be found at http://www.stsci.edu/hst/HST_overview/documents/synphot/AppA_Catalogs9.html
    '''

    return lam,Bnu


# ----- PENDING GETTING THE FORMATTING SORTED OUT -----
def read_albedo_file(filename):
    '''
    Reads an ASCII file that contains the
    wavelengths (lambda in microns) and the 
    product of the geometric albedo and the 
    phase function (gI)
    '''
    
    # Parse the header...
    fin = open(filename,'r')
    lines = fin.readlines()
    
    nheaderlines = 0
    nphi = 1
    nlat = 1
    line = ''

    i = 0
    while lines[i][0] == '#':
        nheaderlines+=1
        line = lines[i]
        
        if str(line[1:6]).upper()=='PHASE' or str(line[1:10]).upper()=='LONGITUDE':
            j = line.find(':')
            line = line[j+1:].split()
            phi = np.zeros(len(line))
            nphi = len(phi)
            for j in range(0,nphi): phi[j] = float(line[j])
            
            # perform a check on the file:
            dphi = phi[1:nphi]-phi[0:nphi-1]
            if min(dphi) < 0:
                print('Error in file {0:s} -- phase/longitude values are not monitonically increasing.'.format(filename))
                exit()
        
        if str(line[1:9]).upper()=='LATITUDE':
            j = line.find(':')
            line = line[j+1:].split()
            lat = np.zeros(len(line))
            nlat = len(lat)
            for j in range(0,nlat): lat[j] = float(line[j])
            
            # perform a check on the file:
            dlat = lat[1:nlat]-lat[0:nlat-1]
            if min(dlat) < 0:
                print('Error in file {0:s} -- latitude values are not monitonically increasing.'.format(filename))
                exit()
        
        i += 1
        
    # Now read the data and header
    lines = lines[nheaderlines:]
    nlambda = len(lines)
    lam = np.zeros(nlambda)

    if nphi==1 and nlat==1:
        gI = np.zeros(nlambda)
        phi = [0.]
        lat = [0.]
    
        for i in range(0,nlambda):
            line = lines[i].split()
            lam[i] = float(line[0])
            gI[i] = float(line[1])

    elif nlat==1:
        gI = np.zeros((nlambda,nphi))
        lat = [0.]
        
        for i in range(0,nlambda):
            line = lines[i].split()
            lam[i] = float(line[0])
            for j in range(0,nphi):
                gI[i,j] = float(line[j+1])

    else:
        gI = np.zeros((nlambda,nphi,nlat))
    
        for i in range(0,nlambda):
            line = lines[i].split()
            lam[i] = float(line[0])
            for j in range(0,nlat):
                for k in range(0,nphi):
                    gI[i,k,j] = float(line[j*nphi+k+1])
        
    # perform a check on the file:
    dlambda = lam[1:]-lam[0:len(lam)-1]
    if np.min(dlambda) < 0:
        print('Error in file {0:d}. Lambda values are not monitonically increasing.'.format(filename))
        return
    return lam,phi,lat,gI
