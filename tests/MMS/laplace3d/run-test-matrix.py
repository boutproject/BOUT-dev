#!/usr/bin/env python3

from boutdata import collect
from boututils.run_wrapper import shell, shell_safe, launch_safe, getmpirun
from matplotlib import pyplot
from sys import exit

MPIRUN = getmpirun()

nproc = 8

print("Making matrix check")
shell_safe("make -f makefile-test-matrix> make.log")

# Delete old data
shell("rm data/BOUT.dmp.*.nc")

args = ""

# Command to run
cmd = "./test-matrix "+args
# Launch using MPI
s, out = launch_safe(cmd, runcmd=MPIRUN, nproc=nproc, pipe=True)

print(out)

x = collect('x', path='data', xguards=True, yguards=True, info=False)
#rhs = collect('rhs', path='data', xguards=True, yguards=True, info=False)
bout_rhs = collect('bout_rhs', path='data', xguards=True, yguards=True, info=False)
petsc_rhs = collect('petsc_rhs', path='data', xguards=True, yguards=True, info=False)

error = petsc_rhs - bout_rhs;

for j in range(1, x.shape[1]-1):
    pyplot.figure('j = '+str(j))

    pyplot.subplot(221)
    pyplot.pcolor(x[1:-1, j, :], cmap='viridis')
    pyplot.colorbar()
    pyplot.title('x')

    pyplot.subplot(222)
    pyplot.pcolor(bout_rhs[1:-1, j, :], cmap='viridis')
    pyplot.colorbar()
    pyplot.title('bout_rhs')

    pyplot.subplot(223)
    pyplot.pcolor(petsc_rhs[1:-1, j, :], cmap='viridis')
    pyplot.colorbar()
    pyplot.title('petsc_rhs')

    pyplot.subplot(224)
    pyplot.pcolor(error[1:-1, j, :], cmap='viridis')
    pyplot.colorbar()
    pyplot.title('error')

    #pyplot.show()
    pyplot.savefig('test-matrix-result.pdf')

exit(0)
