import numpy as np
import pandas as pd
import os
from src import load_stars
from src import generate_planets
from src import generate_disks
from src import generate_scene
from src import read_solarsystem
from src import Settings

# ExoVista v2.4

# Generates a universe of simulated planetary systems based on a stellar target list.

parallel = True
maxcores = 1000000
try: import multiprocessing
except:
    print('multiprocessing module not available. Continuing with serial processing.')
    parallel = False    

if __name__ == '__main__':

    settings = Settings.Settings(timemax=10.0, output_dir='output') # "standard" configuration
    
    #target_list_file = 'master_target_list-usingDR2-50_pc.txt'
    target_list_file = 'target_list8.txt'
    stars, nexozodis = load_stars.load_stars(target_list_file)
    print('\n{0:d} stars in model ranges.'.format(len(stars)))
    
    planets, albedos = generate_planets.generate_planets(stars, settings, addearth=True)
    disks, compositions = generate_disks.generate_disks(stars, planets, settings, nexozodis=nexozodis)
    print('Generating scenes. (This may take a while.)')
    
    if parallel:
        cores = min(maxcores,os.cpu_count())
        cores = min(cores,len(stars))
        percore = int(np.ceil(len(stars)/cores))
        if percore > 1 and percore*(cores-1) >= len(stars): percore -= 1
        
        pool = multiprocessing.Pool(cores)
        
        inputs = []
        for i in range(0,cores):
            imin = i*percore
            imax = (i+1)*percore
            inputs.append([stars.iloc[imin:imax],planets[imin:imax],disks[imin:imax],albedos[imin:imax],compositions[imin:imax],settings])
            
        pool.starmap(generate_scene.generate_scene, [inputs[j][:] for j in range(0,cores)])
        pool.close()
        pool.join()
    else:
        generate_scene.generate_scene(stars,planets,disks,albedos,compositions,settings)

    print('Done')
