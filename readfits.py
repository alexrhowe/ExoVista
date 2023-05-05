import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
from astropy.utils.data import get_pkg_data_filename
import os
import sys
import time

filename = '999-HIP_-TYC_SUN-mv_4.83-L_1.00-d_10.00-Teff_5778.00.fits'

planetcolors = {'Archean_Earth':(0.0,1.0,0.0),'Earth':(0.0,1.0,0.0),'Hazy_Archean_Earth':(0.0,1.0,0.0),'Proterozoic_Earth-hi_o2':(0.0,1.0,0.0),
                'Proterozoic_Earth-lo_o2':(0.0,1.0,0.0),'Venus':(1.0,1.0,0.0),'Mars':(1.0,0.0,0.0),'Early_Mars':(1.0,0.0,0.0),
                'Jupiter':(1.0,0.5,0.0),'Saturn':(0.5,1.0,0.0),'Uranus':(0.0,1.0,1.0),'Neptune':(0.0,0.0,1.0),
                'Warm_Neptune_1AU_Clouds':(0.0,0.0,1.0),'Warm_Neptune_1AU_NoClouds':(0.0,0.0,1.0),'Warm_Neptune_2AU':(0.0,0.0,1.0)}

plt.style.use('dark_background')
plt.rcParams['animation.embed_limit'] = 2**128

lambda_ref    = 0.50  # reference wavelength in microns
mirror_size   = 8.0   # mirror diameter in meters
disk_gain     = 1.0   # multiplies the brightness of the disk image
log_disk      = False # color the disk brightness on a logarithmic scale
planet_bright = False # scale planet brightnesses relative to the maxima of their phase curves
color_code    = False # color-code the planets based on type
fitsdir = './output'  # directory of FITS files

filename = ''

if len(sys.argv)>1:
    filename = sys.argv[1]
    try: hdul = fits.open(filename)
    except:
        print('Error: file not found or in wrong format.')
        filename = ''
        
if filename == '':
    print('Listing available files in FITS file directory...')
    time.sleep(2)
    fitslist = []
    for root, dirs, files in os.walk(fitsdir):
        for file0 in files:
            if file0[-5:] == '.fits': fitslist.append(fitsdir + '/' + file0)
    fitslist = sorted(fitslist)
    
    for i in range(0,len(fitslist)): print(i, fitslist[i]) # print the list with indices if you need it

    while not os.path.exists(filename):
        print('Enter file number or other file name.')
        filenum = input()
        
        try:
            filenum = int(filenum)
            filename = fitslist[filenum]
        except:
            filename = filenum

        try:
            hdul = fits.open(filename)
            hdul.info()
        except:
            print('Error: file not found or in wrong format.')
            filename = ''


# Check whether FITS file is complete.
if len(hdul) < 4:
    print('Error: missing extensions in FITS file.')
    exit()

# Check version number for backwards compatibility.
version = 0.
if 'VERSION' in hdul[0].header:
    version = hdul[0].header['VERSION']

if version <= 2.1:
    specstart = 15
    hstar = 3
else:
    specstart = 16
    hstar = 4

specstart = 16 # temp for the current batch of files
    
speclen = len(hdul[0].data)

if hdul[hstar].data.ndim > 1:
    ntimes = len(hdul[hstar].data)
    npoints = len(hdul[hstar].data[0])
    
nplanets = len(hdul)-hstar
if(abs(hdul[-1].header['A']/np.sqrt(hdul[hstar].header['LSTAR'])-1)<1.e-5): nplanets -= 1 # remove the extra Earth twin if it is present

stardata = np.zeros((ntimes,npoints))
planetdata = np.zeros((nplanets,ntimes,npoints))

stardata = hdul[hstar].data
for i in range(1,nplanets): planetdata[i-1] = hdul[i+hstar].data # note here and later that counting starts from 1 because it's counted from hstar

# Plot 1: Image of disk and planets marked by brightness relative to full phase

ild = np.where(hdul[1].data >= lambda_ref)[0][0]-1
print('Disk brightness plotted at {0:5.1f} nm'.format(hdul[1].data[ild]*1000))
disk = hdul[2].data[ild,:,:]
mindisk = np.log10(np.min(disk))
maxdisk = np.log10(np.max(disk))
if log_disk: disk = (np.log10(disk)-mindisk)/(maxdisk-mindisk) # logarithmic disk brightness for testing

# Compute lambda/D
pixscale = hdul[2].header['PXSCLMAS'] # pixel size in mas
loD = (lambda_ref / (mirror_size*1.e6) * 180./np.pi * 3.6e6) / pixscale

# Black out the central region within 1.5*lambda/D of the star
for i in range(0,250):
    for j in range(0,250):
        if (i-125)**2 + (j-125)**2 <= (loD*1.5)**2: disk[i,j] = np.min(disk)

# Apply gain
disk = np.minimum(disk*disk_gain,np.max(disk))

fig = plt.figure(figsize=(10.8,10.8))
ax = fig.add_subplot(111)
ax.set_xlim(0,250)
ax.set_ylim(0,250)
plt.subplots_adjust(left=0.001, bottom=0.001, right=0.999, top=0.999, wspace=0, hspace=0)

xc = np.arange(0,251)
yc = np.arange(0,251)
ax.pcolor(xc,yc,disk,cmap='inferno')

coords = np.zeros((nplanets,2,ntimes))

for i in range(1,nplanets):
    coords[i-1,0,:] = planetdata[i-1,:,1]
    coords[i-1,1,:] = planetdata[i-1,:,2]

line, = ax.plot([], [])
dots = []

for i in range(0,nplanets):
    lobj = ax.plot([],[],'o',markersize=10,c='w')[0]
    dots.append(lobj)

ilam = np.where(hdul[0].data >= lambda_ref)[0][0]-1
print('Planet brightness plotted at {0:5.1f} nm'.format(hdul[0].data[ilam]*1000))
starbright = np.log10(stardata[0,ilam])
minbright = starbright-15.

albedos = []
plbright = np.zeros((nplanets,ntimes))
for i in range(1,nplanets):
    plbright[i-1] = planetdata[i-1,:,ilam]
    plbright[i-1] /= np.max(plbright[i-1])  # normalize each planet's brightness to the maximum over its orbit
    albedos.append(hdul[i+hstar].header['ALBEDO_F'])

for i in range(1,nplanets):
    x = planetdata[i-1,0,1]
    y = planetdata[i-1,0,2]
    if (x-125)**2 + (y-125)**2 <= (loD*1.5)**2: continue
    bright = plbright[i-1,0]
    dots[i].set_data([x],[y])
    
    # Set planet colors
    if not planet_bright: bright = 1.0 # Sets all planets to maximum brightness
    newcolor = (bright,bright,bright) # Sets all planets to greyscale based on their brightness
    
    # Color-codes planets based on type.
    if color_code and albedos[i] in planetcolors:
        newcolor = tuple([bright*x for x in planetcolors[albedos[i]]])
    
    dots[i].set_color(newcolor)
plt.style.use('default')

# Plot 2: Disk brightness profile along the x-axis
# Recommend to use with PA=0.

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
x = (np.arange(250)-124.5)*0.0796*np.sqrt(2.)
y = np.zeros(250)
for i in range(0,250): y[i] = disk[i,i]*x[i]**2
ax2.set_xlabel('r (AU)')
ax2.set_ylabel('Disk Brightness Profile (cgs)')
ax2.plot(x,y)

# Plot 3: Planetary contrast spectra at t=0

wave = hdul[0].data
fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
for i in range(1,nplanets):
    spec = planetdata[i-1,0,specstart:]
    ax3.plot(wave,spec,label='Planet {0:d}'.format(i))
plt.legend()

# Plot 4: Contrast phase curves of the planets at the reference wavelength

ptime = stardata[:,ilam+15]
time = np.arange(len(ptime))*10

fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
for i in range(1,nplanets):
    ptime = planetdata[i-1,:,ilam+specstart]
    ax4.plot(time,ptime,label='Planet {0:d}'.format(i))
plt.legend()

# Plot 5: 2-D trajectories of the planets on the plane of the sky

fig5 = plt.figure()
ax5 = fig5.add_subplot(111)
for i in range(1,nplanets):
    xtime = planetdata[i-1,:,1]
    ytime = planetdata[i-1,:,2]
    ax5.plot(xtime,ytime,label='Planet {0:d}'.format(i))
    ax5.plot(125,125,'o',c='k',markersize=10)

# Printing a list of transit and eclipse events
    
if version >= 2.2:
    transits = hdul[3].data
    if len(transits[0])==1 and transits[0,0] == 0.:
        print('No transits or eclipses.')
    else:
        print('Transit and Eclipse Event List')
        for i in range(0,len(transits[0])):
            event = ''
            if transits[2,i]+transits[3,i]>0: event='transit'
            if transits[2,i]+transits[3,i]<0: event='eclipse'
            if transits[2,i]==0:
                print('Planet {0:d} {1:s} ingress at {2:.2f} days'.format(int(transits[1,i]),event,transits[0,i]))
            elif transits[3,i]==0:
                print('Planet {0:d} {1:s}  egress at {2:.2f} days'.format(int(transits[1,i]),event,transits[0,i]))
            elif transits[2,i]==transits[3,i]:
                print('Planet {0:d} {1:s} in progress at {2:.2f} days'.format(int(transits[1,i]),event,transits[0,i]))
            else: print('Error in event format at {2:.2f} days'.format(transits[0,i]))

'''
fig.savefig('disk_image.png')
fig2.savefig('disk_profile.png')
fig3.savefig('planet_spectra.png')
fig4.savefig('phase_curves.png')
fig5.savefig('planet_motions.png')
'''
plt.show()
