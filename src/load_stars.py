import math
import numpy as np
from scipy.interpolate import interp1d
import pandas as pd
from src.constants import *

def load_target_list(target_list_file):

    # Read in target list.
    
    fin = open(target_list_file,'r')
    header = fin.readline().replace(',',' ').replace('|',' ').split()
    hlen = len(header)
    
    lines = fin.readlines()
    start = 0
    for line in lines:
        if len(line.replace('|',',').split(','))<hlen:
            start+=1
        else: break
    dlen = len(lines)-start

    grid = []
    target_list = []
    finalheader = []
    headerindex = []
    
    # Clean up header from file.
    for i in range(start,hlen):
        header[i] = header[i].split('(')[0]
        header[i] = header[i].split(',')[0]
        header[i] = header[i].replace('_','')
        
        for h in alias:
            if header[i].lower().startswith(h.lower()):
                if h.lower()=='ra' and header[i].lower().startswith('rad'): continue
                if alias[h] in finalheader: continue
                finalheader.append(alias[h])
                headerindex.append(i)
                break
            else: continue

    if 'Lstar' not in finalheader:
        print('Error: \"Lstar\" missing from input file.')
        exit()
    if 'dist' not in finalheader:
        print('Error: \"dist\" missing from input file.')
        exit()
    if 'Type' not in finalheader:
        print('Error: \"Type\" missing from input file.')
        exit()
    if 'Vmag' not in finalheader:
        print('Error: \"Vmag\" missing from input file.')
        exit()
    if 'Umag' not in finalheader and 'Bmag' not in finalheader:
        print('Error: need a U or B magnitude in input file.')
        exit()
    if 'Rmag' not in finalheader and 'Rmag' not in finalheader and 'Rmag' not in finalheader and 'Rmag' not in finalheader and 'Rmag' not in finalheader:
        print('Error: need an R, I, J, H, or K magnitude in input file.')
        exit()        
            
    for i in range(0,dlen):
        grid.append([])
        line = lines[i+start].replace('|',',').split(',')
        
        for j in range(0,len(headerindex)):
            if finalheader[j] in intlist:
                try: grid[i].append(int(line[headerindex[j]].strip()))
                except:
                    try: grid[i].append(int(float(line[headerindex[j]].strip())))
                    except: grid[i].append(math.nan)
            elif finalheader[j] in strlist:
                grid[i].append(line[headerindex[j]].strip())
            else:
                try: grid[i].append(float(line[headerindex[j]].strip()))
                except: grid[i].append(math.nan)
                
    target_list = pd.DataFrame(grid,columns=finalheader)
    
    # Finished reading target list.
    
    # Compute derived quantities.

    Lstar = target_list['Lstar'].values
    Bmag = target_list['Bmag'].values
    Vmag = target_list['Vmag'].values
    
    BmV = Bmag - Vmag
    # Angular diameter normalized to Vmag=0. From Eq 2 and Table 1 in Boyajian et al. 2014, valid for LC IV and V stars
    logtheta_norm = 0.49612 + 1.11136*BmV - 1.18694*BmV**2. + 0.91974*BmV**3. - 0.19526*BmV**4.
    logtheta = logtheta_norm - 0.2*Vmag
    angdiam = 10.**logtheta # Actual star angular diameter in mass
    
    target_list['BmV'] = BmV
    target_list['angdiam'] = angdiam

    return target_list

def load_stars(target_list_file):

    target_list = load_target_list(target_list_file)
    
    # Filter down target list.
    
    # Eliminate targets with non-finite values for necessary quantities.
    target_list = target_list[np.isfinite(target_list['Lstar']) & np.isfinite(target_list['dist']) & np.isfinite(target_list['Vmag'])]
    
    # Eliminate targets that do not have appropriate comparison magnitudes.
    target_list = target_list[target_list['Umag'].notna() | target_list['Bmag'].notna()]
    target_list = target_list[target_list['Rmag'].notna() | target_list['Imag'].notna() | target_list['Jmag'].notna() | target_list['Hmag'].notna() | target_list['Kmag'].notna()]

    # Limit Lstar to 0-100 Lsun
    target_list = target_list[(target_list['Lstar'] > 0.) & (target_list['Lstar'] < 100.)]

    # Require a Main Sequence or Subgiant spectral type. NOTE: WILL NOT WORK FOR A I/IV BINARY.
    prefixes = ('O','B','A','G','F','G','K','M','o','b','a','f','g','k','m')
    target_list = target_list[target_list['Type'].str.startswith(prefixes,na=False)]
    target_list = target_list[~target_list['Type'].str.contains('II',na=False) & ~(target_list['Type'].str.contains('I',na=False) & ~target_list['Type'].str.contains('IV',na=False))]
    
    # Use Eric Mamajek's table of MS stars to interpolate mass and effective temperature
    fin2 = open('mamajek_dwarf.txt')
    data = fin2.readlines()
    data = [d for d in data if d[0]!='#']
    slen = len(data)
    
    mamajeklogT = []
    mamajeklogL = []
    mamajekBmV = []
    mamajeklogM = []
    
    for i in range(0,slen):
       line = data[i].split()
       try: ttemp = float(line[2])
       except: ttemp = math.nan
       try: ltemp = float(line[3])
       except: ltemp = math.nan
       try: bmvtemp = float(line[7])
       except: bmvtemp = math.nan
       try: mtemp = np.log10(float(line[29]))
       except: mtemp = math.nan
    
       if np.isfinite(ttemp) and np.isfinite(mtemp) and np.isfinite(bmvtemp) and mtemp > -3.2:
           mamajeklogT.append(ttemp)
           mamajeklogL.append(ltemp)
           mamajekBmV.append(bmvtemp)
           mamajeklogM.append(mtemp)

    isort = np.argsort(mamajekBmV)
    mamajeklogT = np.array(mamajeklogT)[isort]
    mamajeklogL = np.array(mamajeklogL)[isort]
    mamajekBmV = np.array(mamajekBmV)[isort]
    mamajeklogM = np.array(mamajeklogM)[isort]
    
    mass_interp = interp1d(mamajekBmV, mamajeklogM, kind='linear', fill_value='extrapolate')
    teff_interp = interp1d(mamajekBmV, mamajeklogT, kind='linear', fill_value='extrapolate')
    # Done reading in reference table.
    
    # Computer interpolated quantities.
    BmV = target_list['BmV'].values
    
    mass = 10. ** mass_interp(BmV)
    Teff = 10. ** teff_interp(BmV)
    
    # Now calculate log(g)
    rsunAU = 0.00465047
    rstarAU = target_list['dist'].values * (target_list['angdiam'].values/1000.)
    rstar = rstarAU / rsunAU / 2.0
    logg = np.log10(mass) - 2.*np.log10(rstar) + 4.437

    target_list['mass'] = mass
    target_list['Teff'] = Teff
    target_list['rstar'] = rstar
    target_list['logg'] = logg
 
    # Dust map can't handle logg > 5.25
    target_list = target_list[target_list['logg'] <= 5.24]

    # Stellar model can't handle Teff < 3500
    target_list = target_list[target_list['Teff'] >= 3500]

    # Fix indices of target list.
    target_list.index = range(len(target_list.index))

    if 'ID' not in target_list.columns:
        target_list['ID'] = np.arange(0,len(target_list))
        
    if 'nexozodis' in target_list.columns:
        nexozodis = target_list['nexozodis'].values
    else: nexozodis = None
    
    return target_list,nexozodis
