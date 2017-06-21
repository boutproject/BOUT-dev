from distutils.core import setup
from Cython.Build import cythonize
import os

os.environ["CXX"] = "mpicxx"
os.environ["CC"]  = "mpicc"

setup(
    ext_modules = cythonize("field3d.pyx",
                            language="c++",             # generate C++ code
                            include_path=['../../include'],
                            )
)
