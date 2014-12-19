####################################
# Replicates IDL's showdata and more
####################################

import numpy as np
import os.path
from boutdata import collect
from boututils import plotpolslice

from tvtk.tools import visual
try:
    from enthought.mayavi import mlab
except ImportError:
    try: from mayavi import mlab
    except ImportError:
        print "No mlab available"
     
from boututils import anim
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


if os.path.isfile('pslice.npy') :
        print 'Reading from file'
	fm=np.load('pslice.npy')
else:
	print 'Compute all time frames'

	for k in xrange(nt):
    		fm[k,:,:]=plotpolslice(p[k,:,:,:],grid,period=period,rz=0)
    
	np.save('pslice',fm)


########################################################
# Set up the window

f = mlab.figure(size=(800,600))
# Tell visual to use this as the viewer.
visual.set_viewer(f)

########################################################
# Do the appropriate graph

#s = mlab.contour_surf(r,z,fun, contours=30, line_width=.5, transparent=True)
#s=mlab.surf(r,z,fun, colormap='Spectral')
s = mlab.mesh(r,z,fm[90,:,:], colormap='PuOr')#, wrap_scale='true')#, representation='wireframe')
s.enable_contours=True
s.contour.filled_contours=True

############ NOTE ############
# In order to get the 3D effect plot only one time frame e.g. mlab.mesh(r,z,fm[10,:,:], colormap='PuOr') and don't do anim()
########################################################


# Define perspective and optional attributes. You can also implement from the window afterwards

mlab.view(0,0)
mlab.draw(f)
mlab.colorbar(orientation="vertical")
#mlab.axes()
#mlab.outline()


########################################################
# mlab animation 

#anim.anim(s,fm)
