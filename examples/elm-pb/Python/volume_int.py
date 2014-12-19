from boututils import file_import, volume_integral
from boutdata import collect
# Integrate over a volume
 

path='./data/'
 
gfile='./cbm18_dens8.grid_nx68ny64.nc'
 
g = file_import(gfile)
   
var=collect("P", path=path)   
   
sol=volume_integral(var, g)

print sol
