from __future__ import print_function
from boututils.file_import import file_import
from boututils.volume_integral import volume_integral
from boutdata.collect import collect
# Integrate over a volume


path='./data/'

gfile='./cbm18_dens8.grid_nx68ny64.nc'

g = file_import(gfile)

var=collect("P", path=path)

sol=volume_integral(var, g)

print(sol)
