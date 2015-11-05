from builtins import range
####################################
# computes polslice for time series
####################################

import numpy as np
from boutdata.collect import collect
from boututils.plotpolslice import plotpolslice

###########################
# Specify parameters


path='./data/'

variable="P"

p = collect(variable, path=path)

period=15

grid='../cbm18_dens8.grid_nx68ny64.nc'

########################################################
# Call plotpolslice once to get extended poloidal grid

r,z,fun=plotpolslice(p[0,:,:,:],grid,period=period,rz=1)

nx=r.shape[0] # number of points in r
ny=r.shape[1] # number of points in z
nt=p.shape[0] # time intervals


fm=np.zeros((nt,nx,ny)) # array to store the time sequence of the poloidal cross section

#Compute all time frames

for k in range(nt):
    fm[k,:,:]=plotpolslice(p[k,:,:,:],grid,period=period,rz=0)

np.savez('pslice',fm=fm, z=z, r=r)
