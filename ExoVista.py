import numpy as np
import pandas as pd
import os
from src import load_stars
from src import generate_planets
from src import generate_disks
from src import generate_scene
from src import read_solarsystem

# Generates a universe of simulated planetary systems based on a stellar target list.

parallel = True
maxcores = 1000000
try: import multiprocessing
except:
    print('multiprocessing module not available. Continuing with serial processing.')
    parallel = False
    

def call_Exo(stars,planets,disks,albedos,compositions):

    kwargs = {'timemax':5.0,'output_dir':'output'} # "standard" configuration
    generate_scene.generate_scene(stars,planets,disks,albedos,compositions,**kwargs)
    

if __name__ == '__main__':
    
    target_list_file = 'master_target_list-usingDR2-50_pc.txt'
    stars, nexozodis = load_stars.load_stars(target_list_file)
    print('\n{0:d} stars in range.'.format(len(stars)))
    
    planets, albedos = generate_planets.generate_planets(stars, imin=0.01, addearth=True, usebins=True) # zero inclination causes numerical instability in the orbit integrator
    
    disks, compositions = generate_disks.generate_disks(stars, planets, nexozodis=nexozodis)
    
    if parallel:
        cores = min(maxcores,os.cpu_count())
        percore = int(np.ceil(len(stars)/cores))
        pool = multiprocessing.Pool(cores)

        if len(stars)<cores: cores = len(stars)
        
        inputs = []
        for i in range(0,cores):
            imin = i*percore
            imax = (i+1)*percore
            inputs.append([stars.iloc[imin:imax],planets[imin:imax],disks[imin:imax],albedos[imin:imax],compositions[imin:imax]])
            
        pool.starmap(call_Exo, [inputs[j][:] for j in range(0,cores)])
        pool.close()
        pool.join()
    else:
        call_Exo(stars,planets,disks,albedos,compositions)

    print('Done')
