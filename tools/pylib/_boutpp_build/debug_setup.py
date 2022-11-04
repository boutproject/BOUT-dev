from distutils.core import setup
from Cython.Build import cythonize
import os

os.environ["CXX"] = "mpic++"
os.environ["CC"] = "mpic++"

setup(
    ext_modules=cythonize(
        "debug.pyx",
        language="c++",  # generate C++ code
    )
)
