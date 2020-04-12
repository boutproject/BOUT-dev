#!/usr/bin/env python3

import os
import sys

def error(*args):
    print(*args)
    sys.exit(1)

def info():
    print(f"""This program, {sys.argv[0]} provides
a convenient way to create modules for the environment module for
BOUT++ that sets PYTHONPATH, PATH and BOUT_TOP, so that switching
between different BOUT++ installations is easy.
{sys.argv[0]} stores the loaded modules as prerequesit, so that you
don't accidentially mix e.g. mpi versions. Feel free to edit the
crated module file, if there are non-needed dependencies.

{sys.argv[0]} can be used interactively, or the location where to
create the modfiles, the name and the BOUT++ directory can be
passed:
  {sys.argv[0]} [modulepath [name [bout-top]]]

If modulepath is not passed, the script checks whether any of the
paths in $MODULEPATH are writeable, and asks which of them to
choose. If non are writeabel, and you can or want not to run as
root, you can extend this. Create a new folder and at it to your
`.bashrc` e.g. with:
```
 mkdir ~/.module
 echo 'export MODULEPATH=$HOME/.module:$MODULEPATH'>> ~/.bashrc
 source ~/.bashrc
```

name can by any name by which you will be later able to recognice
this specific bout installation, e.g. `3d-metrics`. You can later
load the module with `module load bout/name` e.g. `module load
bout/3d-metrics`.

bout-top is the root directory of the BOUT++ installation. Note that
it only works with in-source-builds.""")

def ask_modpath():
    try:
        modpath=os.environ['MODULEPATH']
    except KeyError:
        error("MODULEPATH not set. Please ensure environment modules are loaded.")

    writeable=[]
    for mp in modpath.split(':'):
        if os.access(mp, os.W_OK):
            writeable.append(mp)

    if writeable == []:
        error("Did not find a writeable directory in $MODULEPATH.\nPlease append a writable path.")

    if len(writeable) == 1:
        return mp[0]

    print("Several paths found to store module files. Please select one:")
    for i in range(len(writeable)):
        print("%3d: %s"%(i, writeable[i]))
    select = input("Your choice: ")
    return writeable[int(select)]

def ask_name():
    print("The module will be prefixed with `bout/`")
    print("Please enter the name of the module.")
    name = input("bout/")
    return name

def ask_bouttop():
    this = os.path.realpath(__file__)
    top = "/".join(this.split("/")[:-2])
    print("Do you want to use `%s` as BOUT_TOP?"%top)
    while True:
        inp = input("[Y/n/<alt-path>]")
        if inp in ["", "y", "Y"]:
            return top
        if os.path.isdir(inp):
            return os.path.realpath(inp)
        if inp not in "nN":
            print("`%s` is not a valid directory"%inp)

def create_mod(mp, name, top):
    if not os.path.isdir(mp+"/bout"):
        os.mkdir(mp+"/bout")
    fn = mp+"/bout/"+name
    if 'LOADEDMODULES' in os.environ:
        prereq=[]
        for pr in os.environ['LOADEDMODULES'].split(":"):
            if pr.startswith("bout/"):
                continue
            prereq.append(f"prereq {pr}\n")
        prereq = "".join(prereq)
    else:
        prereq=""
    with open(fn,"w") as f:
        f.write(f"""#%Module 1.0
#
#  BOUT module for use with 'environment-modules' package
#  Created by bout-add-mod-path v0.9

# Only allow one bout module to be loaded at a time
conflict bout
# Require all modules that where loaded at generation time
{prereq}

setenv        BOUT_TOP         {top}
prepend-path  PATH             {top}/bin
prepend-path  LD_LIBRARY_PATH  {top}/lib
prepend-path  PYTHONPATH       {top}/tools/pylib
prepend-path  IDL_PATH        +{top}/tools/idllib:'<IDL_DEFAULT>'
""")
    print(f"created `{fn}`")

if __name__ == "__main__":
    for help in ["-h", "--help", "-?"]:
        if help in sys.argv:
            info()
            sys.exit(0)

    lsa = len(sys.argv)
    mp = ask_modpath() if lsa < 2 else sys.argv[1]
    name = ask_name() if lsa < 3 else sys.argv[2]
    top = ask_bouttop() if lsa < 4 else sys.argv[3]

    create_mod(mp, name, top)
