from setuptools import setup
from Cython.Build import cythonize
import numpy

import os
os.chdir("hortensia_latest")

setup(
    name='hortensia int4c module',
    ext_modules=cythonize("cython/aid4c.pyx"),
    include_dirs=[numpy.get_include()],
    zip_safe=False)

setup(
    name='hortensia calcD module',
    ext_modules=cythonize('cython/calcD.pyx'),
    include_dirs=[numpy.get_include()],
    zip_safe=False)
