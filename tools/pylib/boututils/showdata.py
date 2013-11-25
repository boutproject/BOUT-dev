"""
A Function to animate time dependent data from BOUT++
Requires numpy, mpl_toolkits, matplotlib, boutdata libaries.  

Can now animate more than one plot at once.  To do so, enter your data arrays in a list.  E.g. showdata([var1, var2, var3])

Movies require FFmpeg to be installed.

2D arrays will produce an animated lineplot
3D arrays will produce an animated filled contour plot by default
Other options for 3D data are surface plots and polar plots.  

The tslice variable is used to control the time value that is printed on each
frame of the animation.  If the input data matches the time values found within
BOUT++'s dmp data files, then these time values will be used.  Otherwise, an
integer counter is used.  

Written by Luke Easy
le590@york.ac.uk
Last Updated 25/11/2013
"""

#import numpy as np
from mpl_toolkits.mplot3d import axes3d
from matplotlib import pyplot as plt
from matplotlib import animation
from numpy import linspace, meshgrid, array, min, max, floor, pi
from boutdata import collect
import sys

def showdatatest(vars, titles=[], surf = 0, polar = 0, tslice = 0, movie = 0, intv = 1, Ncolors = 25 ):
    plt.ioff()

    # Check to make sure polar and surface plots aren't both true
    if (surf ==1 & polar ==1):
        print 'What do you want?!? I cant do polar plots AND surface plots at the same time!  Exiting'
        sys.exit()

    # Check to see whether vars is a list or not.
    if isinstance(vars, list):
        Nvar = len(vars)
    else:
        vars = [vars]
        Nvar = len(vars)

    # Sort out titles

    if len(titles) == 0:
        for i in range(0,Nvar):
            titles.append(('Var' + str(i+1)))
        
    elif len(titles) != Nvar:
        print 'The length of the titles input list must match the length of the vars list.  Exiting'
        sys.exit()
    
    # Determine shapes of arrays
    dims = []
    Ndims = []
    for i in range(0,Nvar):
        dims.append(array((vars[i].shape)))
        Ndims.append(dims[i].shape[0])     

        # Perform check to make sure that data is either 2D or 3D
        if (Ndims[i] < 2):
            print 'Error, data must be either 2 or 3 dimensional.  Exiting'
            sys.exit()

        if (Ndims[i] > 3):
            print 'Error, data must be either 2 or 3 dimensional.  Exiting'
            sys.exit()

        if ((Ndims[i] == 2) & (polar != 0)):
            print 'Error, data must be  3 dimensional (time, r, theta) for polar plots.  Exiting'
            sys.exit()

        if ((Ndims[i] != Ndims[0])):
            print 'Error, Number of dimensions must be the same for all variables.  Exiting'
            sys.exit()
    
    # Obtain size of data arrays
    Nt = []
    Nx = []
    Ny = []
    for i in range(0, Nvar):
        Nt.append(vars[i].shape[0])
        Nx.append(vars[i].shape[1])

        if (Nt[i] != Nt[0]):
            print 'Error, Dimensions must be the same for all variables.  Exiting'
            sys.exit()

        if (Nx[i] != Nx[0]):
            print 'Error, Dimensions must be the same for all variables.  Exiting'
            sys.exit()   
        
        if (Ndims[i]== 3):
            Ny.append(vars[i].shape[2])
            if (Ny[i] != Ny[0]):
                print 'Error, Dimensions must be the same for all variables.  Exiting'
                sys.exit()  
            

    # Collect time data from file
    if (tslice == 0):           # Only wish to collect time data if it matches 
        t = collect('t_array')
        
    # Obtain number of frames
    Nframes = Nt[0]/intv

    # Generate grid for plotting
    x = linspace(0,Nx[0]-1, Nx[0])
    if (Ndims[0] == 3):
        y = linspace(0, Ny[0]-1, Ny[0])

    # Determine range of data.  Used to ensure constant colour map and
    # to set y scale of line plot.
    fmax = []
    fmin = []
    clevels = []
    for i in range(0,Nvar):
        fmax.append(max(vars[i]))
        fmin.append(min(vars[i]))
        clevels.append(linspace(fmin[i], fmax[i], Ncolors))
    
    # Create figures for animation plotting
    if (Nvar < 2):
        row = 1
        col = 1
        h = 6.0
        w = 8.0
    elif (Nvar <3):
        row = 1
        col = 2
        h = 6.0
        w = 12.0
    elif (Nvar < 5):
        row = 2
        col = 2
        h = 8.0
        w = 12.0
        
    elif (Nvar < 7):
        row = 2
        col = 3
        h = 8.0
        w = 14.0
        
    elif (Nvar < 10) :
        row = 3
        col = 3
        h = 12.0
        w = 14.0
    else:
        print 'too many variables, for now ;)'
        sys.exit()

    fig = plt.figure(figsize=(w,h))
    ax = []
    if (polar == 0):
        if (Ndims[0] == 2):
            lines = []
            for i in range(0, Nvar):
                ax.append(fig.add_subplot(row,col,i+1))
                ax[i].set_xlim((0,Nx[0]-1))
                ax[i].set_ylim((fmin[i], fmax[i]))
                lines.append(ax[i].plot([],[],lw=2)[0]) 
                #Need the [0] to 'unpack' the line object from tuple.  Alternatively:
                #lines[i], = lines[i]
                ax[i].set_xlabel(r'x')
                ax[i].set_ylabel(titles[i])
                #ax[i].set_ylabel(r'f%i' % (i+1))
                #ax[i].set_title(titles[i])
            title = fig.suptitle(r't = 0', fontsize=14  )

        if (Ndims[0] == 3):
            if (surf == 0): 
                plots = []
                cbars = []
                for i in range(0,Nvar):
                    ax.append(fig.add_subplot(row,col,i+1))
                    ax[i].set_xlim((0,Nx[0]-1))
                    ax[i].set_ylim((0,Ny[0]-1))
                    ax[i].set_xlabel(r'x')
                    ax[i].set_ylabel(r'y')
                    ax[i].set_title(titles[i])
                    plots.append(ax[i].contourf(x,y,vars[i][0,:,:].T, Ncolors, lw=0, levels=clevels[i] ))
                    plt.axes(ax[i])
                    cbars.append(fig.colorbar(plots[i], format='%1.1e', ))
                title = fig.suptitle(r't = 0', fontsize=14)
    
            else:             #Surface plots
                x,y = meshgrid(x,y)
                plots = []
                xstride=int(floor(Nx[0]/20))
                ystride=int(floor(Ny[0]/20))
                for i in range(0,Nvar):
                    ax.append(fig.add_subplot(row,col,i+1, projection='3d'))
                    plots.append(ax[i].plot_wireframe(x, y, vars[i][0,:,:].T, rstride=ystride, cstride=xstride))
                    title = fig.suptitle(r'', fontsize=14 )
                    ax[i].set_xlabel(r'x')
                    ax[i].set_ylabel(r'y')
                    ax[i].set_zlabel(titles[i])

    if (polar != 0):
        r = linspace(1,Nx[0],Nx[0])
        theta = linspace(0,2*pi,Ny[0])
        r, theta = meshgrid(r,theta)
        plots = []
        cbars = []
        for i in range(0,Nvar):
            ax.append(fig.add_subplot(row,col,i+1, projection='polar'))
            plots.append(ax[i].contourf(theta, r, vars[i][0,:,:].T, levels=clevels[i]))
            plt.axes(ax[i])
            cbars.append(fig.colorbar(plots[i], format='%1.1e'))
            ax[i].set_rmax(Nx[0]-1)
            ax[i].set_title(titles[i])
        title = fig.suptitle(r't = 0', fontsize=14  )

    # Animation functions
    def surface1(i):
        index = i*intv
        for j in range(0,Nvar):
            ax[j] = fig.add_subplot(row,col,j+1, projection='3d')
            plots[j] = ax[j].plot_wireframe(x, y, vars[j][i,:,:].T, rstride=ystride, cstride=xstride)
            ax[j].set_zlim(fmin[j],fmax[j])
            ax[j].set_xlabel(r'x')
            ax[j].set_ylabel(r'y')
            ax[j].set_title(titles[j])
        title.set_text('t = %1.2e' % t[index])
        return plots
    
    def surface2(i):
        index = i*intv
        for j in range(0,Nvar):
            ax[j] = fig.add_subplot(row,col,j+1, projection='3d')
            plots[j] = ax[j].plot_wireframe(x, y, vars[j][i,:,:].T, rstride=ystride, cstride=xstride)
            ax[j].set_zlim(fmin[j],fmax[j])
            ax[j].set_xlabel(r'x')
            ax[j].set_ylabel(r'y')
            ax[j].set_title(titles[j])
        title.set_text('t = %1.2e' % t[index])
        return plots
    
    def lineplot1(i):
        index = i*intv
        for j in range(0,Nvar):
            lines[j].set_data(x, vars[j][index,:])
        title.set_text('t = %1.2e' % t[index])
        return lines

    def lineplot2(i):
        index = i*intv
        for j in range(0,Nvar):
            lines[j].set_data(x, vars[j][index,:])
        title.set_text('t = %i' % index)
        return lines

    def contour1(i):
        index = i*intv
        for j in range(0,Nvar):
            plots[j] = ax[j].contourf(x,y,vars[j][index,:,:].T, Ncolors, lw=0, levels=clevels[j])
        title.set_text('t = %1.2e' % t[index])
        return plots

    def contour2(i):
        index = i*intv
        for j in range(0,Nvar):
            plots[j] = ax[j].contourf(x,y,vars[j][index,:,:].T, Ncolors, lw=0, levels=clevels[j])
        title.set_text('t = %i' % index)
        return plots

    def polar1(i):
        index = i*intv
        for j in range(0,Nvar):
            plots[j] = ax[j].contourf(theta, r, vars[j][i,:,:].T, levels=clevels[j])
            ax[j].set_rmax(Nx[0]-1)
        title.set_text('t = %1.2e' % t[index])
        return plots
        
    def polar2(i):
        index = i*intv
        for j in range(0,Nvar):
            plots[j] = ax[j].contourf(theta, r, vars[j][i,:,:].T, levels=clevels[j])
            ax[j].set_rmax(Nx[0]-1)
        title.set_text('t = %i' % index)
        return plots

    # Call Animation functions
    if (polar ==0):
        if (Ndims[0] == 2):
            if (tslice ==0):
                anim = animation.FuncAnimation(fig, lineplot1, frames=Nframes)
            else:
                anim = animation.FuncAnimation(fig, lineplot2, frames=Nframes)

        elif (Ndims[0] == 3):

            if (surf == 0): 
                if (tslice ==0):
                    anim = animation.FuncAnimation(fig, contour1, frames=Nframes)
                else:
                    anim = animation.FuncAnimation(fig, contour2, frames=Nframes)
            else:
                if (tslice ==0):
                    anim = animation.FuncAnimation(fig, surface1, frames=Nframes)
                else:
                    anim = animation.FuncAnimation(fig, surface2, frames=Nframes)

    if (polar !=0):
        if (tslice ==0):
            anim = animation.FuncAnimation(fig, polar1, frames=Nframes)
        else:
            anim = animation.FuncAnimation(fig, polar1, frames=Nframes)

    # Save movie with given name
    if ((isinstance(movie,basestring)==1)):
        anim.save(movie+'.mp4')

    # Save movie with default name
    if ((isinstance(movie,basestring)==0)):
        if (movie != 0):
            anim.save('animation.mp4')

    # Show animation
    if (movie == 0):
        plt.show()   # Save movie with given name
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
1. Speed up animations ????
2. Look at theta in polar plots - perioidic?!?
"""

