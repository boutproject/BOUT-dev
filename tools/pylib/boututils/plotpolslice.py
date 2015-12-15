from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div
import numpy as np
from boututils.file_import import file_import
import sys

if sys.version_info[0]>=3:
    message = "polplotslice uses the VTK library through mayavi, which"+\
              " is currently only available in python 2"
    raise ImportError(message)
else:
    from mayavi import mlab
    from tvtk.tools import visual



def zinterp( v, zind):
  #v = REFORM(v)
  nz = np.size(v)
  z0 = np.round(zind)

  p = zind - float(z0)          # between -0.5 and 0.5

  if p < 0.0 :
      z0 = z0 - 1
      p = p + 1.0


  z0 = ((z0 % (nz-1)) + (nz-1)) % (nz-1)

  # for now 3-point interpolation

  zp = (z0 + 1) % (nz - 1)
  zm = (z0 - 1 + (nz-1)) % (nz - 1)

  result = 0.5*p*(p-1.0)*v[zm] \
    + (1.0 - p*p)*v[z0] \
    + 0.5*p*(p+1.0)*v[zp]

  return result


def plotpolslice(var3d,gridfile,period=1,zangle=0.0, rz=1, fig=0):
    """ data2d = plotpolslice(data3d, 'gridfile' , period=1, zangle=0.0, rz:return (r,z) grid also=1, fig: to do the graph, set to 1 ) """

    g=file_import(gridfile)

    nx=var3d.shape[0]
    ny=var3d.shape[1]
    nz=var3d.shape[2]


    zShift=g.get('zShift')
    rxy=g.get('Rxy')
    zxy=g.get('Zxy')

    dz = 2.0*np.pi / float(period*nz)

    ny2=ny
    nskip=np.zeros(ny-1)
    for i in range(ny-1):
        ip=(i+1)%ny
        nskip[i]=0
        for x in range(nx):
            ns=old_div(np.max(np.abs(zShift[x,ip]-zShift[x,i])),dz)-1
            if ns > nskip[i] : nskip[i] = ns

    nskip = np.int_(np.round(nskip))
    ny2 = np.int_(ny2 + np.sum(nskip))

    print("Number of poloidal points in output:" + str(ny2))

    var2d = np.zeros((nx, ny2))
    r = np.zeros((nx, ny2))
    z = np.zeros((nx, ny2))

    ypos = 0
    for y in range (ny-1) :
      # put in the original points
        for x in range (nx):
            zind = old_div((zangle - zShift[x,y]),dz)
            var2d[x,ypos] = zinterp(var3d[x,y,:], zind)
     #     IF KEYWORD_SET(profile) THEN var2d[x,ypos] = var2d[x,ypos] + profile[x,y]
            r[x,ypos] = rxy[x,y]
            z[x,ypos] = zxy[x,y]

        ypos = ypos + 1

        print((y, ypos))

      # and the extra points

        for x in range (nx):
            zi0 = old_div((zangle - zShift[x,y]),dz)
            zip1 = old_div((zangle - zShift[x,y+1]),dz)

            dzi = old_div((zip1 - zi0), (nskip[y] + 1))

            for i in range (nskip[y]):
                zi = zi0 + float(i+1)*dzi # zindex
                w = old_div(float(i+1),float(nskip[y]+1)) # weighting

                var2d[x,ypos+i] = w*zinterp(var3d[x,y+1,:], zi) + (1.0-w)*zinterp(var3d[x,y,:], zi)
            #  IF KEYWORD_SET(profile) THEN var2d[x,ypos+i] = var2d[x,ypos+i] + w*profile[x,y+1] + (1.0-w)*profile[x,y]
                r[x,ypos+i] = w*rxy[x,y+1] + (1.0-w)*rxy[x,y]
                z[x,ypos+i] = w*zxy[x,y+1] + (1.0-w)*zxy[x,y]



        ypos = ypos + nskip[y]


  # FINAL POINT

    for x in range(nx):
        zind = old_div((zangle - zShift[x,ny-1]),dz)
        var2d[x,ypos] = zinterp(var3d[x,ny-1,:], zind)
     # IF KEYWORD_SET(profile) THEN var2d[x,ypos] = var2d[x,ypos] + profile[x,ny-1]
        r[x,ypos] = rxy[x,ny-1]
        z[x,ypos] = zxy[x,ny-1]


    if(fig==1):

        f = mlab.figure(size=(600,600))
        # Tell visual to use this as the viewer.
        visual.set_viewer(f)


        s = mlab.mesh(r,z,var2d, colormap='PuOr')#, wrap_scale='true')#, representation='wireframe')
        s.enable_contours=True
        s.contour.filled_contours=True
        mlab.view(0,0)

    else:
        # return according to opt
        if rz==1 :
            return r,z,var2d
        else:
            return var2d
