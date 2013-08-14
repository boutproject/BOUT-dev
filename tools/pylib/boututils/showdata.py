"""
A Function to animate time dependent data from BOUT++
Requires numpy, mpl_toolkits, matplotlib, boutdata libaries.  

Movies require FFmpeg to be installed.

2D arrays will produce an animated lineplot
3D arrays will produce an animated filled contour plot by default
Other options for 3D data are surface plots and polar plots.  

The 'tslice' variable is used to control the time value that is printed on each
frame of the animation.  If the input data matches the time values found within
BOUT++'s dmp data files, then these time values will be used.  Otherwise, an
integer counter is used.  

Written by Luke Easy
le590@york.ac.uk
Last Updated 13/8/2013
"""

#import numpy as np
from mpl_toolkits.mplot3d import axes3d
from matplotlib import pyplot as plt
from matplotlib import animation
from numpy import linspace, meshgrid, array, min, max, floor, pi
from boutdata import collect
import sys

def showdata(var, surf = 0, polar = 0, tslice = 0, movie = 0):
    plt.ioff()

    # Determine shape of array
    dims = array(var.shape)
    ndims = dims.shape[0]

    # Perform check to make sure that data is either 2D or 3D
    if (ndims < 2):
        print 'Error, data must be either 2 or 3 dimensional.  Exiting'
        sys.exit()
        
    if (ndims > 3):
        print 'Error, data must be either 2 or 3 dimensional.  Exiting'
        sys.exit()

    if ((ndims == 2) & (polar != 0)):
        print 'Error, data must be  3 dimensional (time, r, theta) for polar plots.  Exiting'
        sys.exit()
        
    # Collect time data from file
    if (tslice == 0):           # Only wish to collect time data if it matches 
        t = collect('t_array')
        
    # Obtain size of data arrays
    Nt = var.shape[0]
    Nx = var.shape[1]
    if (ndims == 3):            # Only need Ny if we have third dimension
        Ny = var.shape[2]

    # Generate grid for plotting
    x = linspace(0, Nx, Nx)
    if (ndims == 3):
        y = linspace(0, Ny, Ny)
        
    # Determine range of data.  Used to ensure constant colour map and
    # to set y scale of line plot.
    fmax = max(var)
    fmin = min(var)

    # Set contour colour levels
    clevels = linspace(fmin, fmax, 25)

    # Create figures for animation plotting

    if (polar == 0):
        if (ndims == 2):
            fig = plt.figure()
            ax = plt.axes(xlim=(0,Nx), ylim=(fmin,fmax))
            line, = ax.plot([], [], lw=2)
            plt.xlabel(r'x')
            plt.ylabel(r'f')

        if (ndims == 3):
            if (surf == 0): 
                fig = plt.figure()
                ax = plt.axes(xlim=(0, Nx), ylim=(0, Ny))
                plt.contourf(x,y,var[0,:,:].T, 25, lw=0, levels=clevels)
                plt.colorbar()
                plt.xlabel(r'x')
                plt.ylabel(r'y')
            else:
                fig = plt.figure()
                x,y = meshgrid(x,y)
                ax = fig.add_subplot(111, projection='3d')
                xstride=int(floor(Nx/20))
                ystride=int(floor(Ny/20))
                
    if (polar != 0):
        r = linspace(0,Nx-1,Nx) 
        theta = linspace(0,2*pi,Ny) 
        r,theta = meshgrid(r,theta)
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        f = var[0,:,:].T
        graph = ax.contourf(theta, r, f, levels=clevels)
        plt.colorbar(graph)
        ax.set_rmax(Nx-1)

    # animation functions
    def lineplot1(i):
        f = var[i,:]
        line.set_data(x,f)
        plt.title(r't = %1.2e' % t[i] )
        return line,
        
    def lineplot2(i):
        f = var[i,:]
        line.set_data(x,f)
        plt.title(r't = %i' % i)
        return line,
    
    def contour1(i): 
        f = var[i,:,:].T
        graph = plt.contourf(x, y, f, lw=0, levels=clevels)
        plt.title(r't = %1.2e' % t[i] )
        return graph

    def contour2(i): 
        f = var[i,:,:].T
        graph = plt.contourf(x, y, f, lw=0, levels=clevels)
        plt.title(r't = %i' % i)
        return graph
    
    def surface1(i): 
        f = var[i,:,:].T
        ax = fig.add_subplot(111, projection='3d')
        graph = ax.plot_wireframe(x, y, f, rstride=ystride, cstride=xstride)
        ax.set_zlim(fmin,fmax)
        plt.title(r't = %1.2e' % t[i] )
        return graph

    def surface2(i): 
        f = var[i,:,:].T
        ax = fig.add_subplot(111, projection='3d')
        graph = ax.plot_wireframe(x, y, f, rstride=ystride, cstride=xstride)
        ax.set_zlim(fmin,fmax)
        plt.title(r't = %i' % i)
        return graph
    
    def polar1(i):
        f = var[i,:,:].T        
        graph = ax.contourf(theta, r, f, levels=clevels)
        ax.set_rmax(Nx-1)
        plt.title(r't = %1.2e' % t[i])
        return graph
        
    def polar2(i):
        f = var[i,:,:]
        graph = ax.contourf(theta, r, f, levels=clevels) 
        ax.set_rmax(Nx-1)
        plt.title(r't = %i' % i)
        return graph
        

    # Call animation functions
    if (polar ==0):
        if (ndims == 2):
            if (tslice ==0):
                anim = animation.FuncAnimation(fig, lineplot1, frames=Nt)
            else:
                anim = animation.FuncAnimation(fig, lineplot2, frames=Nt)
        
        elif (ndims == 3):
                
            if (surf == 0): 
                if (tslice ==0):
                    anim = animation.FuncAnimation(fig, contour1, frames=Nt)
                else:
                    anim = animation.FuncAnimation(fig, contour2, frames=Nt)
            else:
                if (tslice ==0):
                    anim = animation.FuncAnimation(fig, surface1, frames=Nt)
                else:
                    anim = animation.FuncAnimation(fig, surface2, frames=Nt)

    if (polar !=0):
        if (tslice ==0):
            anim = animation.FuncAnimation(fig, polar1, frames=Nt)
        else:
            anim = animation.FuncAnimation(fig, polar1, frames=Nt)
                
    # Save movie with given name
    if ((isinstance(movie,basestring)==1)):
        anim.save(movie+'.mp4')

    # Save movie with default name
    if ((isinstance(movie,basestring)==0)):
        if (movie != 0):
            anim.save('animation.mp4')

    # Show animation
    if (movie == 0):
        plt.show()

"""
To do list
1. Options for number of contour levels
2. Options to output less often
3. Speed up animations
"""


