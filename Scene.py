import astropy.io.fits
import numpy as np
import scipy.interpolate

# Try to find a better way to feed it a filename--star name? There might be duplicates.

class Scene():
    def __init__(self, inputfile):
        self.inputfile = inputfile
        try: self.scene = astropy.io.fits.open(inputfile)
        except:
            print('Error: input file not found.')
            exit()
            
        # Constants defined by the FITS file structure
        self.n_ext = len(self.scene)
        if self.n_ext < 4:
            print('Error missing extensions in FITS file.')
            exit()

        self.version = 0.
        if 'VERSION' in self.scene[0].header:
            self.version = self.scene[0].header['VERSION']            
            
        if self.version <= 2.1:
            self.specstart = 15
            self.hstar = 3
        else:
            self.specstart = 16
            self.hstar = 4
        self.star_ext = self.hstar
        self.planet_ext = self.hstar+1 #first planet extension
    
        self.nplanets = len(self.scene)-self.planet_ext
        if(abs(self.scene[-1].header['A']/np.sqrt(self.scene[self.hstar].header['LSTAR'])-1)<1.e-5):
            self.nplanets -= 1 # remove the extra Earth twin if it is present

        # Define extension numbers
        self.lam_ext = 0
        self.disklam_ext = 1
        self.disk_ext = 2

        # Wavelength arrays
        self.lambdas = astropy.io.fits.getdata(self.inputfile, ext=self.lam_ext) #read wavelength extension
        self.nlambda = len(self.lambdas)
        
        self.lambdas_disk = astropy.io.fits.getdata(self.inputfile, ext=self.disklam_ext) # disk wavelengths
        self.nlambda_disk = len(self.lambdas_disk)

    def getStarLambda(self):
        return self.lambdas
        
    def getDiskLambda(self):
        return self.lambdas_disk

    def getAngDiam(self):
        d, h = astropy.io.fits.getdata(self.inputfile, ext=self.star_ext, header=True)
        return h['ANGDIAM']

    def getPixScale(self):
        d, h = astropy.io.fits.getdata(self.inputfile, ext=self.star_ext, header=True)
        return h['PXSCLMAS']
            
    def getXYstar(self,time=0):
        # Works for an array of times (but is constant anyway).
        d, h = astropy.io.fits.getdata(self.inputfile, ext=self.star_ext, header=True)
        if d.ndim == 1: d = np.expand_dims(d, 0)
        
        t = d[:,0] # time vector
        x = d[:,1] # heliocentric x location vector (pix)
        y = d[:,2] # heliocentric y location vector (pix)
        xstar = x[0] # pick the first entry by default
        ystar = y[0]

        # If the fits file contains a vector of times, interpolate...
        if len(t) > 1: 
            x_interp = scipy.interpolate.interp1d(t, x, kind='quadratic')
            y_interp = scipy.interpolate.interp1d(t, y, kind='quadratic')
            xstar = x_interp(time)
            ystar = y_interp(time)
            
        return np.transpose(np.array([xstar,ystar]))
            
    def getXYstarBary(self,time=0):
        # Works for an array of times (but is constant anyway).
        d, h = astropy.io.fits.getdata(self.inputfile, ext=self.star_ext, header=True)
        if d.ndim == 1: d = np.expand_dims(d, 0)
        
        t = d[:,0] # time vector
        x = d[:,9] # heliocentric x location vector (pix)
        y = d[:,10] # heliocentric y location vector (pix)
        xstar = x[0] # pick the first entry by default
        ystar = y[0]

        # If the fits file contains a vector of times, interpolate...
        if len(t) > 1: 
            x_interp = scipy.interpolate.interp1d(t, x, kind='quadratic')
            y_interp = scipy.interpolate.interp1d(t, y, kind='quadratic')
            xstar = x_interp(time)
            ystar = y_interp(time)
            
        return np.transpose(np.array([xstar,ystar]))
    
    def getRV(self,time=0):
        # Works for an array of times
        d, h = astropy.io.fits.getdata(self.inputfile, ext=self.star_ext, header=True)
        if d.ndim == 1: d = np.expand_dims(d, 0)
        
        t = d[:,0] # time vector
        v = d[:,11] # barycentric z-velocity
        vstar = v[0] # pick the first entry by default

        # If the fits file contains a vector of times, interpolate...
        if len(t) > 1: 
            v_interp = scipy.interpolate.interp1d(t, v, kind='quadratic')
            vstar = v_interp(time)

        rv = np.transpose(np.array(vstar))
        rv *= 1.496e13/31557600 # convert AU/yr to cm/s
            
        return rv

    def getStarSpec(self,time=0):
        d, h = astropy.io.fits.getdata(self.inputfile, ext=self.star_ext, header=True)
        if d.ndim == 1: d = np.expand_dims(d, 0)
        
        t = d[:,0] # time vector
        fstar = d[0,self.specstart:self.specstart+self.nlambda] # grab the stellar flux of first time entry
        
        # If the fits file contains a vector of times, interpolate...
        if len(t) > 1:
            for ii in range(self.nlambda):
                fstar_interp = scipy.interpolate.interp1d(t, d[:,self.specstart+ii], kind='quadratic')                
                fstar[ii] = fstar_interp(time)
                
        return fstar
            
    def getXYplanets(self,time=0):
        # Does not work for an array input.
        xyplanet = np.zeros((self.nplanets,2))
        for ip in range(1,self.nplanets): #loop over all planets
            d, h = astropy.io.fits.getdata(self.inputfile, ext=self.hstar+ip, header=True)
            if d.ndim == 1: d = np.expand_dims(d, 0)
            t = d[:,0] # time vector
            x = d[:,1] # heliocentric x position vector (pix)
            y = d[:,2] # heliocentric y position vector (pix)
            
            xyplanet[ip,0] = x[0] # pick the first entry by default 
            xyplanet[ip,1] = y[0]
            
            if len(t) > 1:
                x_interp = scipy.interpolate.interp1d(t, x, kind='quadratic')
                y_interp = scipy.interpolate.interp1d(t, y, kind='quadratic')
                xyplanet[ip,0] = x_interp(time)
                xyplanet[ip,1] = y_interp(time)

        return xyplanet

    def getPlanetSpec(self,time=0):
        fplanet = np.zeros((self.nplanets,self.nlambda))
        fstar = self.getStarSpec()
        for ip in range(1,self.nplanets): #loop over all planets
            d, h = astropy.io.fits.getdata(self.inputfile, ext=self.hstar+ip, header=True)
            if d.ndim == 1: d = np.expand_dims(d, 0)
            t = d[:,0] # time vector
            
            contrast = d[0,self.specstart:self.specstart+self.nlambda]
            fplanet[ip,:] = contrast * fstar #convert to flux
                
            if len(t) > 1:
                for ii in range(self.nlambda):
                    contrast_interp = scipy.interpolate.interp1d(t, d[:,self.specstart+ii], kind='quadratic')
                    contrast = contrast_interp(time)
                    fplanet[ip, ii] = contrast * fstar[ii]

        return fplanet
    
    def getOrbits(self):
        # orbits: (a, e, i, longnode, argperi, meananom)
        
        orbits = np.zeros((self.nplanets,6))
        for i in range(0,self.nplanets):
            d, h = astropy.io.fits.getdata(self.inputfile, ext=self.planet_ext+i, header=True)
            orbits[i,0] = float(h['A'])
            orbits[i,1] = float(h['E'])
            orbits[i,2] = float(h['I'])
            orbits[i,3] = float(h['LONGNODE'])
            orbits[i,4] = float(h['ARGPERI'])
            orbits[i,5] = float(h['MEANANOM'])

        return orbits

    def countPlanets(self):
        return self.nplanets

    def countPlanetTypes(self):
        fbound = open('planetbins.dat','r')
        line = fbound.readline().strip()
        if line!='Radius':
            print('Error: wrong format in planetbins.dat')
            exit()
        rbound = np.array(fbound.readline().split()).astype(float)
        rlen = len(rbound)-1
        
        atype = fbound.readline().strip()
        alines = fbound.readlines()
        if len(alines)!=len(rbound)-1:
            print('Error: wrong format in planetbins.dat')
            exit()
        alen = len(alines[0].split())
        abound = np.zeros((len(alines),alen))
        for i in range(0,len(alines)):
            abound[i] = np.array(alines[i].split()).astype(float)
        if atype=='Flux': abound = 1./np.sqrt(abound)

        bins = np.zeros((rlen,alen-1))
        d, h = astropy.io.fits.getdata(self.inputfile, ext=self.star_ext, header=True)
        lstar = float(h['LSTAR'])
        
        for i in range(0,self.nplanets):
            d, h = astropy.io.fits.getdata(self.inputfile, ext=self.planet_ext+i, header=True)
            r = float(h['R'])
            a = float(h['A'])/np.sqrt(lstar)
            
            rbin = np.where(rbound>=r)[0][0]-1
            if 0<=rbin<rlen:
                abin = np.where(abound[rbin]>=a)[0][0]-1
                if 0<=abin<alen:
                    bins[rbin,abin] += 1
                
        return bins
    
    def countEECs(self):
        d, h = astropy.io.fits.getdata(self.inputfile, ext=self.star_ext, header=True)
        lstar = float(h['LSTAR'])

        eecs = 0
        
        for i in range(0,self.nplanets):
            d, h = astropy.io.fits.getdata(self.inputfile, ext=self.planet_ext+i, header=True)
            r = float(h['R'])
            a = float(h['A'])/np.sqrt(lstar)
            if 0.95 < a < 1.67 and 0.8/np.sqrt(a) < r < 1.4: eecs += 1
            
        return eecs

    def getTransits(self):
        # Returns a list with four sublists:
        # 1: Time of event in days
        # 2: ID of planet
        # 3: State before event
        # 4: State after event
        # States: 1: in transit, -1: in eclipse, 0: out of transit or eclipse
        # Includes start/end time when a transit/eclipse is in progress
        
        if self.version <= 2.1:
            print('Warning: attempt to retrieve transit data, but transit data are not available.')
            return None
        transits, h = astropy.io.fits.getdata(self.inputfile, ext=3, header=True)
        
        return transits

    def getTransitingPlanets(self):
        
        if self.version <= 2.1:
            print('Warning: attempt to retrieve transiting planets, but transiting planets are not available.')
            return None
        d, h = astropy.io.fits.getdata(self.inputfile, ext=3, header=True)
        plist = np.unique(np.array(d[1]))
        
        if plist == [0]:
            print('No transiting planets found.')
            return None
        else:
            return plist
        
    def getTransitTimes(self, planet):
        # Returns two lists of transit and eclipse times (if they exist) for the given planet ID
        # Each entry within a list is a pair consisting of the ingress and egress time
        # Or the start/end time of the simulation if the event is in progress
        
        if self.version <= 2.1:
            print('Warning: attempt to retrieve transit times, but transit times are not available.')
            return None, None
        d, h = astropy.io.fits.getdata(self.inputfile, ext=3, header=True)
        
        if planet not in d[1]:
            print('Planet {0:d} does not transit.'.format(planet))
            return None, None
        else:
            transits = []
            eclipses = []
            i = np.where(np.array(d[1])==planet)[0]
            
            if d[0][i[0]]==0.:
                if d[0][i[0]]==1: transits.append([0.])
                else: eclipses.append([0.])
                
            for j in i:
                if d[3][j]==d[2][j]:
                    if d[2][j]==1: transits[-1].append(d[0][j])
                    else: eclipses[-1].append(d[j][0])
                elif d[3][j]==1: transits.append([d[0][j]])
                elif d[3][j]==-1: eclipses.append([d[0][j]])
                elif d[2][j]==1: transits[-1].append(d[0][j])
                else: eclipses[-1].append(d[0][j])
                
        return transits, eclipses

    def getDiskImage(self):
        fstar = self.getStarSpec()
        temp = astropy.io.fits.getdata(self.inputfile, ext=self.disk_ext)
        contrast = temp[0:self.nlambda_disk,:,:] #3D contrast data cube
        cprecision = temp[self.nlambda_disk,:,:] #2D contrast precision
        
        #Interpolate the disk image cube to the desired wavelength spacing
        lambda_indices = np.searchsorted(self.lambdas_disk, self.lambdas) - 1
        
        #index in log lambda space (have to add on fractional indices)
        frac_lambda_indices = (lambda_indices + 
                               (np.log(self.lambdas) - np.log(self.lambdas_disk[lambda_indices])) /
                               (np.log(self.lambdas_disk[lambda_indices + 1]) -
                                np.log(self.lambdas_disk[lambda_indices])))
    
        contrast_interp = scipy.interpolate.interp1d(np.arange(len(self.lambdas_disk)), contrast, axis=0, kind='cubic')
        diskimage = np.multiply(contrast_interp(frac_lambda_indices).T, fstar).T

        return diskimage
