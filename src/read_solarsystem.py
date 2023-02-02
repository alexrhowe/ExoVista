import numpy as np
import pandas as pd
from src import generate_disks
from src.constants import *
from src import Settings

settings = Settings.Settings()

def read_solarsystem(settings,system_file='example_system.dat'):

    settings = settings
    fin = open(system_file,'r')
    
    # Set up the star structure with default values
    s = pd.DataFrame(starbase, index=[0])

    l = fin.readline()
    if l.split()[0] != 'Star':
        print('Error: wrong input file format.')
        exit()

    BmV = 0.
    while True:
        try: l = fin.readline()
        except: break
        if l.split()[0] == 'Planets': break
        
        tag = l.split()[0]
        val = l.split()[1]

        if tag in s.columns.values:
            if tag in strlist: s[tag] = val
            elif tag in intlist: s[tag] = int(val)
            else: s[tag] = float(val)
        elif tag=='BmV': BmV = float(val)
    
    s['Vmag'] = s['M_V']-5 + 5*np.log10(s['dist']) # apparent V band mag
    if s['Bmag'].loc[0]==0.: s['Bmag'] = s['Vmag'] + BmV  # apparent B band mag
    s['angdiam'] = 0.465*s['rstar'] / (s['dist']/10.)
    
    # Add the planet structure
    i = 0
    pllist = fin.readline().split()
    maxnplanets = 30
    albedo_files = []
    planet = np.zeros((30,len(pllabel)))

    while i < maxnplanets:
        try: l = fin.readline()
        except: break
        if l.split()[0] == 'Disks': break

        l = l.split()
        for tag in pllist:
            if tag=='albedo': albedo_files.append(l[pllist.index(tag)])
            elif tag in pllabel:
                j = pllabel.index(tag)
                planet[i,j] = float(l[pllist.index(tag)])
            else: continue
        i += 1

    i = 0
    planet = [planet]
    albedo_files = [albedo_files]
    dlist = fin.readline().split()
    
    # Call generate_disks to get the disk structure added
    disk,composition = generate_disks.generate_disks(s,planet,settings)
    # reset all disk params to zero
    disk = np.zeros(disk.shape)[0]
    
    while i < len(disk):
        try: l = fin.readline()
        except: break
        if len(l)==0: break

        l = l.split()
        for tag in dlist:
            if tag=='composition': composition = [l[dlist.index(tag)]]
            elif tag in dlabel:
                j = dlabel.index(tag)
                disk[i,j] = float(l[dlist.index(tag)])
            else: continue
            
        beta = 0.5
        alpha = beta * grav_const * s['mass'] / c          # quantity related to PR drag to be used later
        tpr = disk[i,dlabel.index('r')]**2 / (4 * alpha)                      # yr
        torbit = np.sqrt(disk[i,dlabel.index('r')]**3 / s['mass'].values[0])  # yr
        if i==2: tcoll = torbit / (4*np.pi * (disk[1,dlabel.index('nzodis')] + disk[2,dlabel.index('nzodis')]) * 1e-7) # yr
        else: tcoll = torbit / (4*np.pi * disk[i,dlabel.index('nzodis')] * 1e-7)    # yr
        disk[i,dlabel.index('eta')] = tpr / tcoll
        disk[i,dlabel.index('i')] = s['I']
        
        i += 1
        
    disk = np.array([disk])
    
    return s,planet,albedo_files,disk,composition
