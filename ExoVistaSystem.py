import numpy as np
import pandas as pd
from src import read_solarsystem
from src import generate_scene

# Generates a single, user-defined planetary system.

kwargs = {'timemax':1.0e-10,'output_dir':'output2','diskoff':False,'specres':300.,'configtest':3}
s,p,a,d,c = read_solarsystem.read_solarsystem('solar_system.dat')
generate_scene.generate_scene(s,p,d,a,c,**kwargs)
