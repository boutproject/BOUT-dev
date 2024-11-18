from __future__ import print_function
from boututils.datafile import DataFile
from boututils.volume_integral import volume_integral
from boutdata.collect import collect

# Integrate over a volume


path = "./data/"

gfile = "./cbm18_dens8.grid_nx68ny64.nc"

with DataFile(gfile) as f:
    g = {v: f.read(v) for v in f.keys()}

var = collect("P", path=path)

sol = volume_integral(var, g)

print(sol)
