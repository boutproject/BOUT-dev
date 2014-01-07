#import matplotlib
#matplotlib.use('Qt4Agg')
#from pylab import *

import numpy
from bunch import Bunch
import create_grid
from process_grid import process_grid
from analyse_equil_3 import analyse_equil
from saveobject import saveobject
import cPickle as pickle 
from pylab import figure, show, draw, plot
from read_geqdsk import read_geqdsk

datafile='geqdsk_circular131101.efit'
datafile='g118898.03400'

g=read_geqdsk(datafile)



# plot and check the field
fig=figure(num=0, figsize=(10, 12))
ax = fig.add_subplot(111)
nlev = 100
minf = numpy.min(g.psi)
maxf = numpy.max(g.psi)
levels = numpy.arange(numpy.float(nlev))*(maxf-minf)/numpy.float(nlev-1) + minf
ax.contour(g.r,g.z,g.psi, levels=levels)
ax.set_aspect('equal')
   
show(block=False)

plot(g.rbdry,g.zbdry,'g-')

draw()


npsigrid=numpy.arange(numpy.size(g.pres)).astype(float)/(numpy.size(g.pres)-1)

#fpsi = numpy.zeros((2, numpy.size(fpol)), numpy.float64)
#fpsi[0,:] = (simagx + npsigrid * ( sibdry -simagx ))
#fpsi[1,:] = fpol


rz_grid = Bunch(nr=g.nx, nz=g.ny,   # Number of grid points
                   r=g.r[:,0], z=g.z[0,:], # R and Z as 1D arrays
                   simagx=g.simagx, sibdry=g.sibdry, # Range of psi
                   psi=g.psi, # Poloidal flux in Weber/rad on grid points
                   npsigrid=npsigrid, # Normalised psi grid for fpol, pres and qpsi
                   fpol=g.fpol, # Poloidal current function on uniform flux grid
                   pres=g.pres, # Plasma pressure in nt/m^2 on uniform flux grid
                   qpsi=g.qpsi, # q values on uniform flux grid
                   nlim=g.nlim, rlim=g.xlim, zlim=g.ylim) # Wall boundary




critical = analyse_equil(g.psi,g.r[:,0],g.z[0,:])


settings = Bunch(psi_inner=0.7, 
            psi_outer=1.1, 
            nrad=68, 
            npol=64, 
            rad_peaking=[0.0], 
            pol_peaking=[0.0], 
            parweight=0.0)


boundary = numpy.array([rz_grid.rlim, rz_grid.zlim])

mesh=create_grid.create_grid( g.psi, g.r[:,0], g.z[0,:], settings, 
                                    critical=critical, boundary=boundary,
                    iter=0, fpsi = None, fast='fast')


#save mesh object for faster re-iterations of process_grid (optional)
saveobject(mesh, 'mesh')
#
#read mesh object for faster re-iterations of process_grid (optional)
with open('mesh', 'rb') as input:
    mesh = pickle.load(input)


process_grid( rz_grid, mesh )
