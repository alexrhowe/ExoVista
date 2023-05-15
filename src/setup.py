from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

# Build command:
# python setup.py build_ext --inplace

#import os
#os.environ['CC'] = 'gcc'
#os.environ['CXX'] = 'g++'

ext_modules1=[Extension("wrapImage",sources=["wrapImage.pyx"],language="c++")]
ext_modules2=[Extension("wrapIntegrator",sources=["wrapIntegrator.pyx"],language="c++")]

setup(name="wrapImage",ext_modules=cythonize(ext_modules1))
setup(name="wrapIntegrator",ext_modules=cythonize(ext_modules2))
