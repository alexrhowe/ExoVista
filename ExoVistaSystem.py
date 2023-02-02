import numpy as np
import pandas as pd
from src import read_solarsystem
from src import generate_scene
from src import Settings

# ExoVista v2.1

# Generates a single, user-defined planetary system.

settings = Settings.Settings(output_dir='output', ncomponents=3, timemax=10.0) # "standard" configuration
s,p,a,d,c = read_solarsystem.read_solarsystem(settings,system_file='solar_system.dat')
generate_scene.generate_scene(s,p,d,a,c,settings)
