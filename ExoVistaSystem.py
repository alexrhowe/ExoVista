import numpy as np
import pandas as pd
import sys
from os import path
from src import read_solarsystem
from src import generate_scene
from src import Settings

# ExoVista v2.33

# Generates a single, user-defined planetary system.

filename = ''
filein = ''
ndisk = 2

if len(sys.argv)>1:
    filein = sys.argv[1]
    if not path.exists(filein):
        print('Error: input file not found.')
        filein = ''
        
while not path.exists(filename):
    if filein != '':
        filename = filein
    else:
        print('Enter input file name.')
        filename = input()
        
    if not path.exists(filename):
        print('Error: file not found.')
        filename = ''
        filein = ''
        continue

    line = ''
    
    fin = open(filename,'r')
    line = fin.readline()
    
    if line != 'Star\n':
        print('Error: wrong file format.')
        filename = ''
        filein = ''
        continue
    while line != 'Planets\n':
        try:
            line = fin.readline()
        except:
            print('Error: file contains no planets.')
            filename = ''
            filein = ''
            break
    while line != 'Disks\n':
        try:
            line = fin.readline()
        except:
            print('Error: file contains no disk components.')
            filename = ''
            filein = ''
            break
    ndisk = -1
    while line != 'Settings\n':
        try:
            line = fin.readline()
            if line != 'Settings\n': ndisk += 1
        except: break
        if len(line.split())==0:
            ndisk -= 1
            break

settings = Settings.Settings(output_dir='./output', ncomponents=ndisk, timemax=10.) # "standard" configuration
s,p,a,d,c,new_settings = read_solarsystem.read_solarsystem(settings,system_file=filename)
print('Generating scene...')
generate_scene.generate_scene(s,p,d,a,c,new_settings)
