from builtins import range
from boutdata import collect

path='../data/'

data = collect("P", path=path)

nt=data.shape[0]

ns=data.shape[1]
ne=data.shape[2]
nz=data.shape[3]

    
from enthought.mayavi import mlab
from tvtk.tools import visual
from scipy import mgrid
from enthought.mayavi.mlab import *
f = mlab.figure(size=(600,600))
# Tell visual to use this as the viewer.
visual.set_viewer(f)

#x, y, z= mgrid[0:ns:1, 0:ne:1, 0:nz:1]

#First way

s1 = contour_surf(data[0,:,:,10], contours=30, line_width=.5, transparent=True)
s = surf(data[0,:,:,10], colormap='Spectral')#, wrap_scale='true')#, representation='wireframe')


# second way

#x, y= mgrid[0:ns:1, 0:ne:1]
#s = mesh(x,y,data[0,:,:,10], colormap='Spectral')#, wrap_scale='true')#, representation='wireframe')
#s.enable_contours=True
#s.contour.filled_contours=True


#p=plot3d(x,y,z,data[10,:,:,:], tube_radius=0.025, colormap='Spectral')
#p=points3d(x,y,z,data[10,:,:,:], colormap='Spectral')

#s=contour3d(x,y,z,data[1,:,:,:], contours=4, transparent=True)

#mlab.view(0.,0.)
colorbar()
#axes()
#outline()

@mlab.show
@mlab.animate(delay=250)
def anim():

    for i in range(1,nt):
        s1.mlab_source.scalars = data[i,:,:,10]
        s.mlab_source.scalars = data[i,:,:,10]
   #    mlab.savefig('Movie/anim%d.png'%i)
        yield

# Run the animation.
anim()
