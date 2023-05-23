import numpy as np
import sys
import os
import time
from astropy.io import fits

filename = ''
fitsdir = './output'  # directory of FITS files

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

hip = hdul[hstar].header['HIP']
fout = open('HIP_{0:d}.input.dat'.format(hip),'w')

fout.write('Star\n')
fout.write('ID\t{0:s}\n'.format('999' + str(hdul[hstar].header['ID'])))
fout.write('HIP\t{0:d}\n'.format(hip))
fout.write('dist\t{0:f}\n'.format(hdul[hstar].header['DIST']))
fout.write('M_V\t{0:f}\n'.format(hdul[hstar].header['M_V']))
fout.write('BmV\t{0:f}\n'.format(hdul[hstar].header['BMAG']-hdul[hstar].header['VMAG']))
fout.write('Type\t{0:s}\n'.format(hdul[hstar].header['TYPE']))
fout.write('Lstar\t{0:f}\n'.format(hdul[hstar].header['LSTAR']))
fout.write('logg\t{0:f}\n'.format(hdul[hstar].header['LOGG']))
fout.write('Teff\t{0:f}\n'.format(hdul[hstar].header['TEFF']))
fout.write('mass\t{0:f}\n'.format(hdul[hstar].header['MASS']))
fout.write('rstar\t{0:f}\n'.format(hdul[hstar].header['RSTAR']))
fout.write('PA\t{0:f}\n'.format(hdul[hstar].header['PA']))
fout.write('I\t{0:f}\n'.format(hdul[hstar].header['I']))

fout.write('Planets\n')
fout.write('M\t\tR\t\ta\t\te\t\ti\t\tlongnode\targperi\t\tmeananom\talbedo\n')
for j in range(1,nplanets):
    M = hdul[hstar+j].header['M']
    R = hdul[hstar+j].header['R']
    a = hdul[hstar+j].header['A']
    e = hdul[hstar+j].header['E']
    i = hdul[hstar+j].header['I']
    longnode = hdul[hstar+j].header['LONGNODE']
    argperi = hdul[hstar+j].header['ARGPERI']
    meananom = hdul[hstar+j].header['MEANANOM']
    albedo = hdul[hstar+j].header['ALBEDO_F']
    fout.write('{0:f}\t{1:f}\t{2:f}\t{3:f}\t{4:f}\t{5:f}\t{6:f}\t{7:f}\t{8:s}\n'.format(M,R,a,e,i,longnode,argperi,meananom,albedo))

fout.write('Disks\n')
fout.write('nzodis\t\tr\t\tdror\t\trinner\t\thor\t\tg0\t\tg1\t\tg2\t\tw0\t\tw1\t\tw2\t\tcomposition\n')
for i in range(0,hdul[2].header['NCOMP']):
    nzodis = hdul[2].header['NZODIS-'+str(i)]
    r = hdul[2].header['R-'+str(i)]
    dror = hdul[2].header['DROR-'+str(i)]
    rinner = hdul[2].header['RINNER-'+str(i)]
    hor = hdul[2].header['HOR-'+str(i)]
    g0 = hdul[2].header['G0-'+str(i)]
    g1 = hdul[2].header['G1-'+str(i)]
    g2 = hdul[2].header['G2-'+str(i)]
    w0 = hdul[2].header['W0-'+str(i)]
    w1 = hdul[2].header['W1-'+str(i)]
    w2 = hdul[2].header['W2-'+str(i)]
    composition = hdul[2].header['COMPOSIT']
    fout.write('{0:f}\t{1:f}\t{2:f}\t{3:f}\t{4:f}\t{5:f}\t{6:f}\t{7:f}\t{8:f}\t{9:f}\t{10:f}\t{11:s}\n'.format(nzodis,r,dror,rinner,hor,g0,g1,g2,w0,w1,w2,composition))
    
fout.write('Settings\n')
fout.write('pixscale\t{0:f}\n'.format(hdul[2].header['PXSCLMAS']*0.001))
if version >= 2.3: fout.write('iwa\t\t{0:f}\n'.format(hdul[2].header['IWA']))
fout.write('npix\t\t{0:d}\n'.format(hdul[2].header['NAXIS1']))
fout.write('specres\t\t{0:f}\n'.format(hdul[0].header['SPECRES']))
fout.write('specresdisk\t{0:f}\n'.format(hdul[1].header['SPECRES']))
fout.write('lammin\t\t{0:f}\n'.format(hdul[0].header['LAMMIN']))
fout.write('lammax\t\t{0:f}\n'.format(hdul[0].header['LAMMAX']))
if version >= 2.3: fout.write('rdust_blowout\t{0:f}\n'.format(hdul[2].header['DUSTBLOW']))
if version >= 2.3: fout.write('tsublimate\t{0:f}\n'.format(hdul[2].header['TSUB']))
