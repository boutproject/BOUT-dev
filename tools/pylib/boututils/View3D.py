"""
View a 3D rendering of the magnetic field lines and the streamlines of the rational surfaces.
The quality of the later can be used as an indicator of the quality of the grid. The magnetic field
is computed from efit_analyzed.py. The script can be used as a template to show additional properties of the field

based on enthought's example by Gael Varoquaux <gael.varoquaux@normalesup.org>
https://docs.enthought.com/mayavi/mayavi/auto/example_magnetic_field.html#example-magnetic-field

"""
from __future__ import absolute_import
from __future__ import division
from builtins import range
from past.utils import old_div


from boutdata.collect import collect
import numpy as np

import sys

if sys.version_info[0]>=3:
    message = "View3D uses the VTK library through mayavi, which"+\
              " is currently only available in python 2"
    raise ImportError(message)
else:
    from mayavi import mlab

from .read_geqdsk import read_geqdsk
from boututils.View2D import View2D
from scipy import interpolate
from .boutgrid import *


def View3D(g,path=None, gb=None):
  ##############################################################################
  # Resolution

  n=51

  #compute Bxy
  [Br,Bz,x,y,q]=View2D(g,option=1)


  rd=g.r.max()+.5
  zd=g.z.max()+.5
  ##############################################################################
  # The grid of points on which we want to evaluate the field
  X, Y, Z = np.mgrid[-rd:rd:n*1j, -rd:rd:n*1j, -zd:zd:n*1j]
  ## Avoid rounding issues :
  #f = 1e4  # this gives the precision we are interested by :
  #X = np.round(X * f) / f
  #Y = np.round(Y * f) / f
  #Z = np.round(Z * f) / f

  r = np.c_[X.ravel(), Y.ravel(), Z.ravel()]

  ##############################################################################
  # Calculate field
  # First initialize a container matrix for the field vector :
  B = np.empty_like(r)


  #Compute Toroidal field
  # fpol is given between simagx (psi on the axis) and sibdry (
  # psi on limiter or separatrix). So the toroidal field (fpol/R) and the q profile are within these boundaries
  # For each r,z we have psi thus we get fpol if (r,z) is within the boundary (limiter or separatrix) and fpol=fpol(outer_boundary) for outside

  #The range of psi is g.psi.max(), g.psi.min() but we have f(psi) up to the limit. Thus we use a new extended variable padded up to max psi
  # set points between psi_limit and psi_max

  add_psi=np.linspace(g.sibdry,g.psi.max(),10)

  # define the x (psi) array
  xf=np.arange(np.float(g.qpsi.size))*(g.sibdry-g.simagx)/np.float(g.qpsi.size-1) + g.simagx

  # pad the extra values excluding the 1st value

  xf=np.concatenate((xf, add_psi[1::]), axis=0)

  # pad fpol with corresponding points

  fp=np.lib.pad(g.fpol, (0,9), 'edge')

  # create interpolating function

  f = interpolate.interp1d(xf, fp)

  #calculate Toroidal field

  Btrz = old_div(f(g.psi), g.r)


  rmin=g.r[:,0].min()
  rmax=g.r[:,0].max()
  zmin=g.z[0,:].min()
  zmax=g.z[0,:].max()


  B1p,B2p,B3p,B1t,B2t,B3t = magnetic_field(g,X,Y,Z,rmin,rmax,zmin,zmax, Br,Bz,Btrz)

  bpnorm = np.sqrt(B1p**2 + B2p**2 + B3p**2)
  btnorm = np.sqrt(B1t**2 + B2t**2 + B3t**2)

  BBx=B1p+B1t
  BBy=B2p+B2t
  BBz=B3p+B3t
  btotal = np.sqrt(BBx**2 + BBy**2 + BBz**2)

  Psi = psi_field(g,X,Y,Z,rmin,rmax,zmin,zmax)

  ##############################################################################
  # Visualization

  # We threshold the data ourselves, as the threshold filter produce a
  # data structure inefficient with IsoSurface
  #bmax = bnorm.max()
  #
  #B1[B > bmax] = 0
  #B2[B > bmax] = 0
  #B3[B > bmax] = 0
  #bnorm[bnorm > bmax] = bmax

  mlab.figure(1, size=(1080,1080))#, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5))

  mlab.clf()

  fieldp = mlab.pipeline.vector_field(X, Y, Z, B1p, B2p, B3p,
                                    scalars=bpnorm, name='Bp field')

  fieldt = mlab.pipeline.vector_field(X, Y, Z, B1t, B2t, B3t,
                                    scalars=btnorm, name='Bt field')

  field = mlab.pipeline.vector_field(X, Y, Z, BBx, BBy, BBz,
                                    scalars=btotal, name='B field')



  field2 = mlab.pipeline.scalar_field(X, Y, Z, Psi, name='Psi field')

  #vectors = mlab.pipeline.vectors(field,
  #                      scale_factor=1,#(X[1, 0, 0] - X[0, 0, 0]),
  #                      )

  #vcp1 = mlab.pipeline.vector_cut_plane(fieldp,
  #                                        scale_factor=1,
  #                                        colormap='jet',
  #                                        plane_orientation='y_axes')
  ##
  #vcp2 = mlab.pipeline.vector_cut_plane(fieldt,
  #                                        scale_factor=1,
  #                                        colormap='jet',
  #                                        plane_orientation='x_axes')


  # Mask random points, to have a lighter visualization.
  #vectors.glyph.mask_input_points = True
  #vectors.glyph.mask_points.on_ratio = 6

  #vcp = mlab.pipeline.vector_cut_plane(field1)
  #vcp.glyph.glyph.scale_factor=5*(X[1, 0, 0] - X[0, 0, 0])
  # For prettier picture:
  #vcp1.implicit_plane.widget.enabled = False
  #vcp2.implicit_plane.widget.enabled = False

  iso = mlab.pipeline.iso_surface(field2,
                                  contours=[Psi.min()+.01],
                                  opacity=0.4,
                                  colormap='bone')

  for i in range(q.size):
      iso.contour.contours[i+1:i+2]=[q[i]]

  iso.compute_normals = True
  #

  #mlab.pipeline.image_plane_widget(field2,
  #                            plane_orientation='x_axes',
  #                            #slice_index=10,
  #                            extent=[-rd, rd, -rd, rd, -zd,zd]
  #                        )
  #mlab.pipeline.image_plane_widget(field2,
  #                            plane_orientation='y_axes',
  #                           # slice_index=10,
  #                            extent=[-rd, rd, -rd,rd, -zd,zd]
  #                        )



  #scp = mlab.pipeline.scalar_cut_plane(field2,
  #                                        colormap='jet',
  #                                        plane_orientation='x_axes')
  # For prettier picture and with 2D streamlines:
  #scp.implicit_plane.widget.enabled = False
  #scp.enable_contours = True
  #scp.contour.number_of_contours = 20

  #

  # Magnetic Axis

  s=mlab.pipeline.streamline(field)
  s.streamline_type = 'line'
  s.seed.widget = s.seed.widget_list[3]
  s.seed.widget.position=[g.rmagx,0.,g.zmagx]
  s.seed.widget.enabled = False


  #  q=i surfaces

  for i in range(np.shape(x)[0]):

      s=mlab.pipeline.streamline(field)
      s.streamline_type = 'line'
  ##s.seed.widget = s.seed.widget_list[0]
  ##s.seed.widget.center = 0.0, 0.0, 0.0
  ##s.seed.widget.radius = 1.725
  ##s.seed.widget.phi_resolution = 16
  ##s.seed.widget.handle_direction =[ 1.,  0.,  0.]
  ##s.seed.widget.enabled = False
  ##s.seed.widget.enabled = True
  ##s.seed.widget.enabled = False
  #
      if x[i].size>1 :
          s.seed.widget = s.seed.widget_list[3]
          s.seed.widget.position=[x[i][0],0.,y[i][0]]
          s.seed.widget.enabled = False


  # A trick to make transparency look better: cull the front face
  iso.actor.property.frontface_culling = True

  #mlab.view(39, 74, 0.59, [.008, .0007, -.005])
  out=mlab.outline(extent=[-rd, rd, -rd, rd, -zd, zd], line_width=.5 )
  out.outline_mode = 'cornered'
  out.outline_filter.corner_factor = 0.0897222


  w = mlab.gcf()
  w.scene.camera.position = [13.296429046581462, 13.296429046581462, 12.979811259697154]
  w.scene.camera.focal_point = [0.0, 0.0, -0.31661778688430786]
  w.scene.camera.view_angle = 30.0
  w.scene.camera.view_up = [0.0, 0.0, 1.0]
  w.scene.camera.clipping_range = [13.220595435695394, 35.020427055647517]
  w.scene.camera.compute_view_plane_normal()
  w.scene.render()
  w.scene.show_axes = True

  mlab.show()

  if(path is not None):
  #BOUT data
  #path='../Aiba/'
  #
  #gb = file_import(path+'aiba.bout.grd.nc')
  #gb = file_import("../cbm18_8_y064_x516_090309.nc")
  #gb = file_import("cbm18_dens8.grid_nx68ny64.nc")
  #gb = file_import("/home/ben/run4/reduced_y064_x256.nc")

    data = collect('P', path=path)
    data = data[50,:,:,:]
  #data0=collect("P0", path=path)
  #data=data+data0[:,:,None]

    s = np.shape(data)
    nz = s[2]


    sgrid = create_grid(gb, data, 1)

  # OVERPLOT the GRID
  #mlab.pipeline.add_dataset(sgrid)
  #gr=mlab.pipeline.grid_plane(sgrid)
  #gr.grid_plane.axis='x'


  ## pressure scalar cut plane from bout
    scpb = mlab.pipeline.scalar_cut_plane(sgrid,
                                          colormap='jet',
                                          plane_orientation='x_axes')

    scpb.implicit_plane.widget.enabled = False
    scpb.enable_contours = True
    scpb.contour.filled_contours=True
  #
    scpb.contour.number_of_contours = 20
  #
  #
  #loc=sgrid.points
  #p=sgrid.point_data.scalars

  # compute pressure from scatter points interpolation
  #pint=interpolate.griddata(loc, p, (X, Y, Z), method='linear')
  #dpint=np.ma.masked_array(pint,np.isnan(pint)).filled(0.)
  #
  #p2 = mlab.pipeline.scalar_field(X, Y, Z, dpint, name='P field')
  #
  #scp2 = mlab.pipeline.scalar_cut_plane(p2,
  #                                        colormap='jet',
  #                                        plane_orientation='y_axes')
  #
  #scp2.implicit_plane.widget.enabled = False
  #scp2.enable_contours = True
  #scp2.contour.filled_contours=True
  #scp2.contour.number_of_contours = 20
  #scp2.contour.minimum_contour=.001



  # CHECK grid orientation
  #fieldr = mlab.pipeline.vector_field(X, Y, Z, -BBx, BBy, BBz,
  #                                  scalars=btotal, name='B field')
  #
  #sg=mlab.pipeline.streamline(fieldr)
  #sg.streamline_type = 'tube'
  #sg.seed.widget = sg.seed.widget_list[3]
  #sg.seed.widget.position=loc[0]
  #sg.seed.widget.enabled = False



  #OUTPUT grid

  #ww = tvtk.XMLStructuredGridWriter(input=sgrid, file_name='sgrid.vts')
  #ww.write()

  return

def magnetic_field(g,X,Y,Z,rmin,rmax,zmin,zmax,Br,Bz,Btrz):

    rho = np.sqrt(X**2 + Y**2)
    phi=np.arctan2(Y,X)

    br=np.zeros(np.shape(X))
    bz=np.zeros(np.shape(X))
    bt=np.zeros(np.shape(X))

    nx,ny,nz=np.shape(X)

    mask = (rho >= rmin) & (rho <= rmax) & (Z >= zmin) & (Z <= zmax)
    k=np.argwhere(mask==True)

    fr=interpolate.interp2d(g.r[:,0], g.z[0,:], Br.T)
    fz=interpolate.interp2d(g.r[:,0], g.z[0,:], Bz.T)
    ft=interpolate.interp2d(g.r[:,0], g.z[0,:], Btrz.T)

    for i in range(len(k)):
        br[k[i,0],k[i,1],k[i,2]]=fr(rho[k[i,0],k[i,1],k[i,2]],Z[k[i,0],k[i,1],k[i,2]])
        bz[k[i,0],k[i,1],k[i,2]]=fz(rho[k[i,0],k[i,1],k[i,2]],Z[k[i,0],k[i,1],k[i,2]])
        bt[k[i,0],k[i,1],k[i,2]]=ft(rho[k[i,0],k[i,1],k[i,2]],Z[k[i,0],k[i,1],k[i,2]])

 # Toroidal component
    B1t=-bt*np.sin(phi)
    B2t=bt*np.cos(phi)
    B3t=0*bz

 # Poloidal component
    B1p=br*np.cos(phi)
    B2p=br*np.sin(phi)
    B3p=bz


    # Rotate the field back in the lab's frame
    return B1p,B2p,B3p,B1t,B2t,B3t


def psi_field(g,X,Y,Z,rmin,rmax,zmin,zmax):

    rho = np.sqrt(X**2 + Y**2)

    psi=np.zeros(np.shape(X))

    nx,ny,nz=np.shape(X)

    mask = (rho >= rmin) & (rho <= rmax) & (Z >= zmin) & (Z <= zmax)
    k=np.argwhere(mask==True)

    f=interpolate.interp2d(g.r[:,0], g.z[0,:], g.psi.T)

    for i in range(len(k)):
        psi[k[i,0],k[i,1],k[i,2]]=f(rho[k[i,0],k[i,1],k[i,2]],Z[k[i,0],k[i,1],k[i,2]])

    # Rotate the field back in the lab's frame
    return psi


if __name__ == '__main__':
     path='../../tokamak_grids/pyGridGen/'
     g=read_geqdsk(path+"g118898.03400")
     View3D(g)
     mlab.show()
