from __future__ import division
from future import standard_library
standard_library.install_aliases()
from past.utils import old_div
import numpy as np
from boututils.bunch import Bunch
from create_grid  import create_grid
from process_grid import process_grid
from analyse_equil_3 import analyse_equil
from saveobject import saveobject
import pickle as pickle 
from pylab import figure, show, draw, plot
from read_geqdsk import read_geqdsk


def grid(g,psi_in=None,psi_out=None,nrd=None,npl=None):

  
  
  
  # plot and check the field
  fig=figure(num=0, figsize=(10, 12))
  ax = fig.add_subplot(111)
  nlev = 100
  minf = np.min(g.psi)
  maxf = np.max(g.psi)
  levels = np.arange(np.float(nlev))*(maxf-minf)/np.float(nlev-1) + minf
  ax.contour(g.r,g.z,g.psi, levels=levels)
  ax.set_aspect('equal')
     
  show(block=False)
  
  plot(g.xlim,g.ylim,'g-')
  
  draw()
  
  
  npsigrid=old_div(np.arange(np.size(g.pres)).astype(float),(np.size(g.pres)-1))
  
  #fpsi = numpy.zeros((2, np.size(fpol)), np.float64)
  #fpsi[0,:] = (simagx + npsigrid * ( sibdry -simagx ))
  #fpsi[1,:] = fpol
  
  rz_grid = Bunch(nr=g.nx,          # Number of grid points in R direction           [1]
                  nz=g.ny,          # Number of grid points in Z direction           [1]
                  r=g.r[:,0],       # 1D array of R position                         [m]   
                  z=g.z[0,:],       # 1D array of Z position                         [m]
                  simagx=g.simagx,  # poloidal flux function at the axis             [Wb/rad]
                  sibdry=g.sibdry,  # poloidal flux function at the plasma boundary  [Wb/rad]
                  psi=g.psi,        # 2D array of poloidal flux on RZ plane          [Wb/rad] 
                  npsigrid=npsigrid,# uniform flux grid normalised in [0,1] 
                  fpol=g.fpol,      # poloidal current function on uniform flux grid [T*m]
                  pres=g.pres,      # plasma pressure in nt/m^2 on uniform flux grid [Pa]
                  qpsi=g.qpsi,      # q values on uniform flux grid                  [1]
                  nlim=g.nlim,      # Number of wall boundary points                 [1]
                  rlim=g.xlim,      # 1D array of wall bounday points in R           [m]   
                  zlim=g.ylim)      # 1D array of wall bounday points in Z           [m]   
    
  critical = analyse_equil(g.psi,g.r[:,0],g.z[0,:]) # find X/O-points
  
  #settings = Bunch(psi_inner=psi_in,  # inner radial boundary defined by npsi
  #                 psi_outer=psi_out, # outer radial boundary defined by npsi
  #                 nrad=nrd,          # number of radial   grid points 
  #                 npol=npl,          # number of poloidal grid points 
  #                 rad_peaking=[0.0], # radial separatrix peaking factor
  #                 pol_peaking=[0.0], # poloidal separatrix peaking factor
  #                 parweight=0.0)     # for separatrix 
  
  settings = Bunch(psi_inner=psi_in,  # inner radial boundary defined by npsi
                   psi_outer=psi_out, # outer radial boundary defined by npsi
                   nrad=[nrd],        # number of radial   grid points as 1D array 
                   npol=[npl],        # number of poloidal grid points as 1D array
                   rad_peaking=[0.0], # radial separatrix peaking factor
                   pol_peaking=[0.0], # poloidal separatrix peaking factor
                   parweight=0.0)     # for separatrix  (not used yet)  

  boundary = np.array([rz_grid.rlim, rz_grid.zlim])

  mesh= create_grid( g.psi,             # 2D array of poloidal flux [Wb/rad] 
                     g.r[:,0],          # 1D array of R position    [m]   
                     g.z[0,:],          # 1D array of R position    [m]   
                     settings,          # settings for grid generation    
                     critical=critical, # information of X/O-points
                     boundary=boundary, # position of PFCs 
                     iter=0,            # number of iteration
                     fpsi = None,       # 2D array of current function [T*m]
                     fast='fast')       # use fast option

  #save mesh object for faster re-iterations of process_grid (optional)
  saveobject(mesh, 'mesh')
  
  #read mesh object for faster re-iterations of process_grid (optional)
  with open('mesh', mode='rb') as input: mesh = pickle.load(input)

  process_grid( rz_grid, mesh )

  return


if __name__ == '__main__':


  datafile='g118898.03400'
  #datafile='g014220.00200'
  g=read_geqdsk(datafile)

  grid(g,0.6,0.8,64,128)
