import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
from astropy.utils.data import get_pkg_data_filename
import os

filename = 'output/-1-HIP_-TYC_SUN-mv_4.83-L_1.00-d_10.00-Teff_5778.00.fits'

planetcolors = {'Archean_Earth':(0.0,1.0,0.0),'Earth':(0.0,1.0,0.0),'Hazy_Archean_Earth':(0.0,1.0,0.0),'Proterozoic_Earth-hi_o2':(0.0,1.0,0.0),
                'Proterozoic_Earth-lo_o2':(0.0,1.0,0.0),'Venus':(1.0,1.0,0.0),'Mars':(1.0,0.0,0.0),'Early_Mars':(1.0,0.0,0.0),
                'Jupiter':(1.0,0.5,0.0),'Saturn':(0.5,1.0,0.0),'Uranus':(0.0,1.0,1.0),'Neptune':(0.0,0.0,1.0),
                'Warm_Neptune_1AU_Clouds':(0.0,0.0,1.0),'Warm_Neptune_1AU_NoClouds':(0.0,0.0,1.0),'Warm_Neptune_2AU':(0.0,0.0,1.0)}

plt.style.use('dark_background')
plt.rcParams['animation.embed_limit'] = 2**128

lambda_ref = 0.50 # reference wavelength in microns

'''
# Find available stars to plot
fitsdir = './output'
fitslist = []
for root, dirs, files in os.walk(fitsdir):
    for file0 in files:
        if file0[-5:] == '.fits': fitslist.append(fitsdir + '/' + file0)
fitslist = sorted(fitslist)

for i in range(0,len(fitslist)): print(i, fitslist[i]) # print the list with indices if you need it
hdul = fits.open(fitslist[0])
'''

hdul = fits.open(filename)
hdul.info()

if len(hdul) < 4:
    print('Error missing extensions in FITS file.')
    exit()
notime = False
if hdul[3].data.ndim > 1:
    ntimes = len(hdul[3].data)
    npoints = len(hdul[3].data[0])
else:
    notime = True
    ntimes = 1
    npoints = len(hdul[3].data)
nplanets = len(hdul)-4
if(abs(hdul[-1].header['A']-1)<1.e-5): nplanets -= 1 # remove the extra Earth twin if it is present

stardata = np.zeros((ntimes,npoints))
planetdata = np.zeros((nplanets,ntimes,npoints))

if notime:
    stardata = np.expand_dims(hdul[3].data, axis=0)
    for i in range(0,nplanets): planetdata[i] = np.expand_dims(hdul[i+4].data, axis=0)
else:
    stardata = hdul[3].data
    for i in range(0,nplanets): planetdata[i] = hdul[i+4].data

# Plot 1: Image of disk and planets marked by brightness relative to full phase

ild = np.where(hdul[1].data >= lambda_ref)[0][0]-1
print('Disk brightness plotted at {0:5.1f} nm'.format(hdul[1].data[ild]*1000))
disk = hdul[2].data[ild,:,:]
mindisk = np.log10(np.min(disk))
maxdisk = np.log10(np.max(disk))
#disk = (np.log10(disk)-mindisk)/(maxdisk-mindisk) # logarithmic disk brightness for testing

# Compute lambda/D
pixscale = hdul[2].header['PXSCLMAS'] # pixel size in mas
D = 6.7 # mirror diameter in meteres
loD = (lambda_ref / (D*1.e6) * 180./np.pi * 3.6e6) / pixscale

# Black out the central region within 1.5*lambda/D of the star
for i in range(0,250):
    for j in range(0,250):
        if (i-125)**2 + (j-125)**2 <= (loD*1.5)**2: disk[i,j] = np.min(disk)

gain = 1.0 # amplify the disk at the cost of saturating the center
disk = np.minimum(disk*gain,np.max(disk))

fig = plt.figure(figsize=(10.8,10.8))
ax = fig.add_subplot(111)
ax.set_xlim(0,250)
ax.set_ylim(0,250)
plt.subplots_adjust(left=0.001, bottom=0.001, right=0.999, top=0.999, wspace=0, hspace=0)

xc = np.arange(0,250)
yc = np.arange(0,250)
ax.pcolor(xc,yc,disk,cmap='inferno')

coords = np.zeros((nplanets,2,ntimes))

for i in range(0,nplanets):
    coords[i,0,:] = planetdata[i,:,1]
    coords[i,1,:] = planetdata[i,:,2]

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
for i in range(0,nplanets):
    plbright[i] = planetdata[i,:,ilam]
    plbright[i] /= np.max(plbright[i])  # normalize each planet's brightness to the maximum over its orbit
    if i>0: albedos.append(hdul[i+4].header['ALBEDO_F'])
plbright /= np.max(plbright[1:,0])

for i in range(0,nplanets):
    x = planetdata[i,0,1]
    y = planetdata[i,0,2]
    if (x-125)**2 + (y-125)**2 <= (loD*1.5)**2: continue
    bright = plbright[i,0]
    dots[i].set_data([x],[y])

    # Set planet colors
    newcolor = (bright,bright,bright) # Sets all planets to greyscale based on their brightness
    #if ntimes > 1: newcolor=(1.0,1.0,1.0) # Sets all planets to white.
    # Color-codes planets based on type.
    '''
    if i>0:
        newcolor = tuple([bright*x for x in planetcolors[albedos[i-1]]])
        print(i,albedos[i-1])
    newcolor = tuple([bright*x for x in newcolor])
    '''
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
for i in range(0,nplanets):
    spec = planetdata[i,0,15:]
    ax3.plot(wave,spec,label='Planet {0:d}'.format(i))
plt.legend()

# Plot 4: Contrast phase curves of the planets at the reference wavelength
# Will be blank if there is only one timestep.

ptime = stardata[:,ilam+15]
time = np.arange(len(ptime))*10

fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
for i in range(0,nplanets):
    ptime = planetdata[i,:,ilam+15]
    ax4.plot(time,ptime,label='Planet {0:d}'.format(i))
plt.legend()

#fig.savefig('disk_image.png')
#fig2.savefig('disk_profile.png')
#fig3.savefig('planet_spectra.png')
#fig4.savefig('phase_curves.png')
plt.show()
