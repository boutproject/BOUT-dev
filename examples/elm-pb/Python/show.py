####################################
# Replicates IDL's showdata and more
####################################

import numpy as np
from boutdata import collect
from boututils import plotpolslice
from enthought.mayavi import mlab
from tvtk.tools import visual

###########################
# Specify parameters


path='../data/'

variable="P"

p = collect(variable, path=path)

period=15

grid='../cbm18_dens8.grid_nx68ny64.nc'

########################################################
# Call plotpolslice once to get extended poloidal grid 

r,z,fun=plotpolslice(p[0,:,:,:],grid,period=period,rz=1)

nx=r.shape[0] # number of points in r
ny=r.shape[1] # number of points in z
nt=10#fun.shape[0] # time intervals


fm=np.zeros((nt,nx,ny)) # array to store the time sequence of the poloidal cross section

#Compute all time frames

for k in xrange(1,nt):
    fm[k,:,:]=plotpolslice(p[k,:,:,:],grid,period=period,rz=0)

########################################################
# Set up the window

f = mlab.figure(size=(800,600))
# Tell visual to use this as the viewer.
visual.set_viewer(f)

########################################################
# Do the appropriate graph

#s = mlab.contour_surf(r,z,fun, contours=30, line_width=.5, transparent=True)
#s=mlab.surf(r,z,fun, colormap='Spectral')
s = mlab.mesh(r,z,fun, colormap='PuOr')#, wrap_scale='true')#, representation='wireframe')
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


@mlab.show
@mlab.animate(delay=250)
def anim():

    for i in range(1,nt):#nt+1):
        s.mlab_source.scalars = fm[i,:,:]
        title="t="+np.string0(i)
        mlab.title(title,height=1.1, size=0.26)
    #    mlab.savefig('../Movie/anim%d.png'%i)  # uncomment to save images for movie creation
        yield

 #Run the animation.
anim()
