import numpy as np

# Physical constants
grav_const    = 4.*np.pi**2 # in AU^3 yr^-2 solar_mass^-1
c = 2.998e8 / 1.496e11 * 365.25 * 24. * 60. * 60. # AU yr^-1
planck = 6.62608e-27

# Data structure parameters
maxnplanets   = 30
mincomponents = 1 # Minimum number of disk components, warm dust / exo-zodi required
maxcomponents = 3 # Maximum number of disk components, highest plausible is 3-4.

# File paths
exovistapath  = './'
lqq_dir       = 'lqq_files/'

# Define debris disk grain sizes
'''
If we did Mie theory on the fly using dustmap, we'd have to use
a LOT of different grain sizes to average over Mie ringing
artifacts. Instead, we have pre-calculated Qabs and Qsca
for different grain size ranges (in lqq_files folder).
Here is the master list of pre-calculated grain sizes we can
select from. I previously determined a size resolution ~ 5 was
sufficient to accurately reproduce Qsca and Qabs for any distribution
of grain sizes.
'''

# Variables used in theory to generate the dust filenames procedurally.
'''
sizeres = 5.
master_maxsize = 1000.
master_minsize = 0.1
dlnsize = 1./sizeres
master_lnmaxsize = np.log(master_maxsize)
master_lnminsize = np.log(master_minsize)
master_nsizes = int(np.ceil((master_lnmaxsize-master_lnminsize)/dlnsize))
'''

# Arrays needed to generate the exact dust filenames.
master_rdust = np.array([0.1103, 0.1342, 0.1632, 0.1986, 0.2415, 0.2938, 0.3574, 0.4348, 0.5289, 0.6434, 0.7827, 0.9522, 1.1583, 1.4091, 1.7141, 2.0852, 2.5366, 3.0858, 3.7538, 4.5664, 5.5550, 6.7575, 8.2204, 10.0000, 12.1648, 14.7983, 18.0019, 21.8991, 26.6399, 32.4070, 39.4225, 47.9569, 58.3388, 70.9682, 86.3317, 105.0211, 127.7565, 155.4138, 189.0584, 229.9864, 279.7748, 340.3416, 414.0202, 503.6490, 612.6804, 745.3159, 906.6649])

master_rdust_boundaries = np.array([0.1000, 0.1216, 0.1480, 0.1800, 0.2190, 0.2664, 0.3241, 0.3942, 0.4796, 0.5834, 0.7097, 0.8633, 1.0502, 1.2776, 1.5541, 1.8906, 2.2999, 2.7977, 3.4034, 4.1402, 5.0365, 6.1268, 7.4532, 9.0666, 11.0294, 13.4171, 16.3217, 19.8551, 24.1534, 29.3822, 35.7430, 43.4808, 52.8937, 64.3444, 78.2739, 95.2190, 115.8323, 140.9082, 171.4126, 208.5206, 253.6620, 308.5758, 375.3775, 456.6408, 555.4957, 675.7517, 822.0413, 1000.0001])

master_nsizes = len(master_rdust_boundaries)
master_drdust = master_rdust_boundaries[1:master_nsizes]-master_rdust_boundaries[0:master_nsizes-1]

# Table headings
starbase = {'ID':0, 'HIP':0, 'TYC2':'', 'dist':10., 'M_V':0., 'Vmag':0., 'Bmag':0., 'Umag':float('nan'), 'Rmag':float('nan'), 'Imag':float('nan'), 'Jmag':float('nan'), 'Hmag':float('nan'), 'Kmag':float('nan'), 'Type':'', 'Lstar':0., 'logg':0., 'Teff':0., 'angdiam':0., 'mass':0., 'rstar':0., 'PA':0., 'I':60.}
intlist = ['ID', 'HIP']
strlist = ['TYC2', 'WDS', 'Type']
keplist = ('a','e','i','longnode','argperi','meananom')
pllabel = ('M','R','a','e','i','longnode','argperi','meananom')
dlabel = ('n', 'longnode', 'i', 'nzodis', 'r', 'dror', 'rinner', 'eta', 'hor', 'g0', 'g1', 'g2', 'w0', 'w1', 'w2')

# FITS file comments
scomments = {'PA':'System midplane position angle (deg)',
             'I':'System midplane inclination (deg)',
             'ID':'Internal catalog #',
             'HIP':'Hipparcos designation',
             'TYC2':'Tycho-2 designation',
             'dist':'Distance (pc)',
             'Vmag':'Absolute V band mag',
             'Type':'Spectral type',
             'Lstar':'Bolometric luminosity (solar luminosities)',
             'Teff':'Effective temperature (K)',
             'angdiam':'Angular diameter (mas)',
             'mass':'Mass (solar masses)',
             'Rstar':'Radius (solar radii)'}

pcomments = {'M':'Mass (Earth masses)',
             'R':'Radius (Earth radii)',
             'a':'Semi-major axis (AU)',
             'e':'Eccentricity',
             'i':'Inclination (deg)',
             'longnode':'Longitude of ascending node (deg)',
             'argperi':'Argument of pericenter (deg)',
             'meananom':'Mean anomaly (deg)'}

dcomments = {'longnode':'longitude of ascending node (deg)',
             'i':'inclination (deg)',
             'nzodis':'exozodi level (zodis)',
             'r':'peak density radius (AU)',
             'dror':'Gaussian peak width',
             'rinner':'truncation radius (AU)',
             'eta':'T_PR/T_coll',
             'hor':'scale height',
             'g0':'HG scattering asym. parameter 1',
             'w0':'HG function weight 1',
             'g1':'HG scattering asym. parameter 2',
             'w1':'HG function weight 2',
             'g2':'HG scattering asym. parameter 3',
             'w2':'HG function weight 3'}
