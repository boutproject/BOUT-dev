from __future__ import print_function
####################################
# Replicates IDL's showdata and more
####################################

import numpy as np

from tvtk.tools import visual
try:
    from enthought.mayavi import mlab
except ImportError:
    try: from mayavi import mlab
    except ImportError:
        print("No mlab available")

from boututils.anim import anim
###########################
# Read polslice array

npzfile=np.load('pslice.npz')
r=npzfile['r']
z=npzfile['z']
fm=npzfile['fm']

########################################################
# Set up the window

f = mlab.figure(size=(800,600))
# Tell visual to use this as the viewer.
visual.set_viewer(f)

########################################################
# Do the appropriate graph

#s = mlab.contour_surf(r,z,fun, contours=30, line_width=.5, transparent=True)
#s=mlab.surf(r,z,fun, colormap='Spectral')
s = mlab.mesh(r,z,fm[0,:,:], scalars=fm[0,:,:], colormap='PuOr')#, wrap_scale='true')#, representation='wireframe')
s.enable_contours=True
s.contour.filled_contours=True

# Define perspective and optional attributes. You can also implement from the window afterwards

mlab.view(0,0)
#mlab.view(-94.159958841373324,
# 53.777002382688906,
# 8.2311808018087582)
mlab.draw(f)
mlab.colorbar(orientation="vertical")
#mlab.axes()
#mlab.outline()


########################################################
# mlab animation

anim(s,fm, save=True)
