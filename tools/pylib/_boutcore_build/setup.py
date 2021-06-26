from setuptools import Extension, setup
from Cython.Build import cythonize
import os
import numpy

from boutconfig import config as conf

libs = []
ldirs = []
lextra = []
for l in conf["libs"].split():
    if l.startswith("-l"):
        libs.append(l[2:])
    elif l.startswith("-L"):
        ldirs.append(l[2:])
    else:
        lextra.append(l)

include = []
flags = []
for i in conf["cflags"].split():
    if i.startswith("-I"):
        include.append(i[2:])
    else:
        flags.append(i)
include.append(numpy.get_include())

extensions = [
    Extension(
        "boutcore",
        ["boutcore.pyx"],
        include_dirs=include,
        libraries=libs,
        library_dirs=ldirs,
        extra_compile_args=flags,
        extra_link_args=lextra,
    )
]

# Setting LDSHARED sets the linker. However, setting it might drop
# flags - thus it is not set by default. If the compiler used here is
# different to the python compiler, LDSHARED might need to be set. In
# that case it might be required to check the linking commands, and
# add the apropriate flags.

setup(
    ext_modules=cythonize(
        extensions,
        language_level=3,
        compiler_directives=dict(binding=True, embedsignature=True),
    )
)
