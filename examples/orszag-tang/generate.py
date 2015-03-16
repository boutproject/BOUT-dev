from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
# the Scientific Python netCDF 3 interface
# http://dirac.cnrs-orleans.fr/ScientificPython/
#from Scientific.IO.NetCDF import NetCDFFile as Dataset
# the 'classic' version of the netCDF4 python interface
# http://code.google.com/p/netcdf4-python/
import numpy as np
from netCDF4 import Dataset
from numpy import dtype # array module from http://numpy.scipy.org
"""
This example writes some surface pressure and temperatures
The companion program sfc_pres_temp_rd.py shows how to read the netCDF
data file created by this program.

This example demonstrates the netCDF Python API.
It will work either with the Scientific Python NetCDF version 3 interface
(http://dirac.cnrs-orleans.fr/ScientificPython/)
of the 'classic' version of the netCDF4 interface. 
(http://netcdf4-python.googlecode.com/svn/trunk/docs/netCDF4_classic-module.html)
To switch from one to another, just comment/uncomment the appropriate
import statements at the beginning of this file.

Jeff Whitaker <jeffrey.s.whitaker@noaa.gov> 20070202
"""
# Adapted for BOUT++ by
# George Breyiannis, JAEA, Nov 2013
#

# the output array to write will be nx x ny
ny = 100; nx = ny + 4

# dy of grid
dy = old_div(1.0, np.float(ny))
dx = dy
# create grid
dxarr=np.zeros((nx,ny),dtype='float32')+dx
dyarr=np.zeros((nx,ny),dtype='float32')+dy

xarr=np.arange(0.,np.float(nx),1.,dtype='float32')*dx
yarr=np.arange(0.,np.float(ny),1.,dtype='float32')*dy

# compute initial variables

rho=np.zeros((nx,ny),dtype='float32')+old_div(25.,(36.*np.pi))
p=np.zeros((nx,ny),dtype='float32')+old_div(5.,(12.*np.pi))

rho=1.
p=old_div(rho,3.)

v_x=np.zeros((nx,ny),dtype='float32')
Bx=np.zeros((nx,ny),dtype='float32')

for y in range(ny):
    v_x[:,y]=-np.sin(2.*np.pi*yarr[y])
    Bx[:,y]=-np.sin(2.*np.pi*yarr[y])
 
#Bx=Bx/np.sqrt(4.*np.pi)


v_y=np.zeros((nx,ny),dtype='float32')
By=np.zeros((nx,ny),dtype='float32')

for x in range(nx):
    v_y[x,:]=np.sin(2.*np.pi*xarr[x])
    By[x,:]=np.sin(4.*np.pi*xarr[x])
    
#By=By/np.sqrt(4.*np.pi)


# Domain inside core (periodic)

ixseps1 = nx
ixseps2 = nx

# open a new netCDF file for writing.
ncfile = Dataset('otv.grd.128.nc','w', format='NETCDF3_CLASSIC') 

# output data.

# create the nx and ny dimensions.
ncfile.createDimension('x',nx)
ncfile.createDimension('y',ny)
ncfile.createDimension('single',1)

# create and write nx,ny variables
nxx=ncfile.createVariable('nx','i4',('single'))
nyy=ncfile.createVariable('ny','i4',('single'))

nxx[:]=nx
nyy[:]=ny

# Define the coordinate variables. They will hold the coordinate
# information, that is, xarr,yarr
dx = ncfile.createVariable('dx',dtype('float32').char,('x','y'))
dy = ncfile.createVariable('dy',dtype('float32').char,('x','y',))

# write data to coordinate vars.
dx[:,:] = dxarr
dy[:,:] = dyarr

# create and write ixseps* dimensions.
ix1=ncfile.createVariable('ixseps1','i4',('single'))
ix2=ncfile.createVariable('ixseps2','i4',('single'))

ix1[:]=ixseps1
ix2[:]=ixseps2

# create the corresponding variables 
rho0 = ncfile.createVariable('rho0',dtype('float32').char,('x','y'))
p0 = ncfile.createVariable('p0',dtype('float32').char,('x','y'))
v0_x = ncfile.createVariable('v0_x',dtype('float32').char,('x','y'))
v0_y = ncfile.createVariable('v0_y',dtype('float32').char,('x','y'))
B0x = ncfile.createVariable('B0x',dtype('float32').char,('x','y'))
B0y = ncfile.createVariable('B0y',dtype('float32').char,('x','y'))

# write data to variables.

rho0[:,:]=rho
p0[:,:]=p
v0_x[:,:]=v_x
v0_y[:,:]=v_y
B0x[:,:]=Bx
B0y[:,:]=By


ncfile.close()
print('*** SUCCESS writing file otv.grd.py.nc!')
