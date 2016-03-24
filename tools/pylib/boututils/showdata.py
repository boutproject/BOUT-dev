"""
Visualisation and animation routines

Written by Luke Easy
le590@york.ac.uk
Last Updated 19/3/2015
Additional functionality by George Breyiannis 26/12/2014

"""
from __future__ import print_function
from __future__ import division
try:
    from builtins import str
    from builtins import chr
    from builtins import range
except:
    pass

#import numpy as np
from mpl_toolkits.mplot3d import axes3d
from matplotlib import pyplot as plt
from matplotlib import animation
from numpy import linspace, meshgrid, array, min, max, floor, pi
from boutdata.collect import collect


####################################################################
# Specify manually ffmpeg path
#plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

FFwriter = animation.FFMpegWriter()
####################################################################


###################
#http://stackoverflow.com/questions/16732379/stop-start-pause-in-python-matplotlib-animation
#
j=-2
pause = False
###################


def showdata(vars, titles=[], legendlabels = [], surf = [], polar = [], tslice = 0, movie = 0, intv = 1, Ncolors = 25, x = [], y = [], global_colors = False, symmetric_colors = False,hold_aspect=False):
    """
    A Function to animate time dependent data from BOUT++
    Requires numpy, mpl_toolkits, matplotlib, boutdata libaries.

    To animate multiple variables on different axes:
    showdata([var1, var2, var3])

    To animate more than one line on a single axes:
    showdata([[var1, var2, var3]])

    The default graph types are:
    2D (time + 1 spatial dimension) arrays = animated line plot
    3D (time + 2 spatial dimensions) arrays = animated contour plot.

    To use surface or polar plots:
    showdata(var, surf = 1)
    showdata(var, polar = 1)

    Can plot different graph types on different axes.  Default graph types will be used depending on the dimensions of the input arrays.  To specify polar/surface plots on different axes:
    showdata([var1,var2], surf = [1,0], polar = [0,1])

    Movies require FFmpeg to be installed.

    The tslice variable is used to control the time value that is printed on each
    frame of the animation.  If the input data matches the time values found within
    BOUT++'s dmp data files, then these time values will be used.  Otherwise, an
    integer counter is used.

    During animation click once to stop in the current frame. Click again to continue.

    global_colors = True: if "vars" is a list the colorlevels are determined from the mximum of the maxima and and the minimum of the  minima in all fields in vars.

    symmetric_colors = True: colorlevels are symmetric.
    """
    plt.ioff()

    # Check to see whether vars is a list or not.
    if isinstance(vars, list):
        Nvar = len(vars)
    else:
        vars = [vars]
        Nvar = len(vars)

    if Nvar < 1:
        raise ValueError("No data supplied")

    # Check to see whether each variable is a list - used for line plots only
    Nlines = []
    for i in range(0, Nvar):
        if isinstance(vars[i], list):
            Nlines.append(len(vars[i]))
        else:
            Nlines.append(1)
            vars[i] = [vars[i]]

    # Sort out titles
    if len(titles) == 0:
        for i in range(0,Nvar):
            titles.append(('Var' + str(i+1)))
    elif len(titles) != Nvar:
        raise ValueError('The length of the titles input list must match the length of the vars list.')

    # Sort out legend labels
    if len(legendlabels) == 0:
        for i in range(0,Nvar):
            legendlabels.append([])
            for j in range(0,Nlines[i]):
                legendlabels[i].append(chr(97+j))
    elif (isinstance(legendlabels[0], list) != 1):
        if Nvar != 1:
            check = 0
            for i in range(0,Nvar):
                if len(legendlabels) != Nlines[i]:
                    check = check+1
            if check == 0:
                print("Warning, the legendlabels list does not contain a sublist for each variable, but it's length matches the number of lines on each plot. Will apply labels to each plot")
                legendlabelsdummy = []
                for i in range(0, Nvar):
                    legendlabelsdummy.append([])
                    for j in range(0,Nlines[i]):
                        legendlabelsdummy[i].append(legendlabels[j])
                legendlabels = legendlabelsdummy
            else:
                print("Warning, the legendlabels list does not contain a sublist for each variable, and it's length does not match the number of lines on each plot. Will default apply labels to each plot")
                legendlabels = []
                for i in range(0,Nvar):
                    legendlabels.append([])
                    for j in range(0,Nlines[i]):
                        legendlabels[i].append(chr(97+j))
        else:
            if (Nlines[0] == len(legendlabels)):
                legendlabels = [legendlabels]
    elif len(legendlabels) != Nvar:
        print("Warning, the length of the legendlabels list does not match the length of the vars list, will continue with default values")
        legendlabels = []
        for i in range(0,Nvar):
            legendlabels.append([])
            for j in range(0,Nlines[i]):
                legendlabels[i].append(chr(97+j))
    else:
        for i in range(0,Nvar):
            if isinstance(legendlabels[i], list):
                if len(legendlabels[i]) != Nlines[i]:
                    print('Warning, the length of the legendlabel (sub)list for each plot does not match the number of datasets for each plot. Will continue with default values')
                legendlabels[i] = []
                for j in range(0,Nlines[i]):
                    legendlabels[i].append(chr(97+j))
            else:
                legendlabels[i] = [legendlabels[i]]
            if len(legendlabels[i]) != Nlines[i]:
                print('Warning, the length of the legendlabel (sub)list for each plot does not match the number of datasets for each plot.  Will continue with default values')
                legendlabels[i] = []
                for j in range(0,Nlines[i]):
                    legendlabels[i].append(chr(97+j))


    # Sort out surf list
    if isinstance(surf, list):
        if (len(surf) == Nvar):
            for i in range(0, Nvar):
                if surf[i] >= 1:
                    surf[i] = 1
                else:
                    surf[i] = 0
        elif (len(surf) == 1):
            if surf[0] >= 1:
                surf[0] = 1
            else:
                surf[0] = 0
            if (Nvar > 1):
                for i in range(1,Nvar):
                    surf.append(surf[0])
        elif (len(surf) == 0):
            for i in range(0,Nvar):
                surf.append(0)
        else:
            print('Warning, length of surf list does not match number of variables.  Will default to no polar plots')
            for i in range(0,Nvar):
                surf.append(0)

    else:
        surf = [surf]
        if surf[0] >= 1:
            surf[0] = 1
        else:
            surf[0] = 0
        if (Nvar > 1):
            for i in range(1,Nvar):
                surf.append(surf[0])

    # Sort out polar list
    if isinstance(polar, list):
        if (len(polar) == Nvar):
            for i in range(0, Nvar):
                if polar[i] >= 1:
                    polar[i] = 1
                else:
                    polar[i] = 0
        elif (len(polar) == 1):
            if polar[0] >= 1:
                polar[0] = 1
            else:
                polar[0] = 0
            if (Nvar > 1):
                for i in range(1,Nvar):
                    polar.append(polar[0])
        elif (len(polar) == 0):
            for i in range(0,Nvar):
                polar.append(0)
        else:
            print('Warning, length of polar list does not match number of variables.  Will default to no polar plots')
            for i in range(0,Nvar):
                polar.append(0)
    else:
        polar = [polar]
        if polar[0] >= 1:
            polar[0] = 1
        else:
            polar[0] = 0
        if (Nvar > 1):
            for i in range(1,Nvar):
                polar.append(polar[0])

    # Determine shapes of arrays
    dims = []
    Ndims = []
    lineplot = []
    contour = []
    for i in range(0,Nvar):
        dims.append([])
        Ndims.append([])
        for j in range(0, Nlines[i]):
            dims[i].append(array((vars[i][j].shape)))
            Ndims[i].append(dims[i][j].shape[0])
            # Perform check to make sure that data is either 2D or 3D
            if (Ndims[i][j] < 2):
                raise ValueError('data must be either 2 or 3 dimensional.  Exiting')

            if (Ndims[i][j] > 3):
                raise ValueError('data must be either 2 or 3 dimensional.  Exiting')

            if ((Ndims[i][j] == 2) & (polar[i] != 0)):
                print('Warning, data must be  3 dimensional (time, r, theta) for polar plots.  Will plot lineplot instead')

            if ((Ndims[i][j] == 2) & (surf[i] != 0)):
                print('Warning, data must be  3 dimensional (time, x, y) for surface plots.  Will plot lineplot instead')

            if ((Ndims[i][j] == 3) & (Nlines[i] != 1)):
                raise ValueError('cannot have multiple sets of 3D (time + 2 spatial dimensions) on each subplot')


            if ((Ndims[i][j] != Ndims[i][0])):
                raise ValueError('Error, Number of dimensions must be the same for all variables on each plot.')

        if (Ndims[i][0] == 2): # Set polar and surf list entries to 0
            polar[i] = 0
            surf[i] = 0
            lineplot.append(1)
            contour.append(0)
        else:
            if ((polar[i] == 1) & (surf[i] == 1)):
                print('Warning - cannot do polar and surface plots at the same time.  Default to contour plot')
                contour.append(1)
                lineplot.append(0)
                polar[i] = 0
                surf[i] = 0
            elif (polar[i] == 1) | (surf[i] == 1):
                contour.append(0)
                lineplot.append(0)
            else:
                contour.append(1)
                lineplot.append(0)

    # Obtain size of data arrays
    Nt = []
    Nx = []
    Ny = []
    for i in range(0, Nvar):
        Nt.append([])
        Nx.append([])
        Ny.append([])
        for j in range(0, Nlines[i]):
            Nt[i].append(vars[i][j].shape[0])
            Nx[i].append(vars[i][j].shape[1])
            if (Nt[i][j] != Nt[0][0]):
                raise ValueError('time dimensions must be the same for all variables.')

            #if (Nx[i][j] != Nx[i][0]):
            #    raise ValueError('Dimensions must be the same for all variables on each plot.')

            if (Ndims[i][j] == 3):
                Ny[i].append(vars[i][j].shape[2])
                #if (Ny[i][j] != Ny[i][0]):
                #    raise ValueError('Dimensions must be the same for all variables.')

    # Collect time data from file
    if (tslice == 0):           # Only wish to collect time data if it matches
        try:
            t = collect('t_array')
            if t == None:
                raise ValueError("t_array is None")
            if len(t) != Nt[0][0]:
                raise ValueError("t_array is wrong size")
        except:
            t = linspace(0,Nt[0][0], Nt[0][0])

    # Obtain number of frames
    Nframes = int(Nt[0][0]/intv)

    # Generate grids for plotting
    # Try to use provided grids where possible 
    xnew = []
    ynew = []
    for i in range(0,Nvar):
        xnew.append([])
        try:
            xnew[i].append(x[i])
        except:	
            for j in range(0, Nlines[i]):
                xnew[i].append(linspace(0,Nx[i][j]-1, Nx[i][j]))

        #x.append(linspace(0,Nx[i][0]-1, Nx[i][0]))
        
        if (Ndims[i][0] == 3):
            try:
                ynew.append(y[i])
            except:        
                ynew.append(linspace(0, Ny[i][0]-1, Ny[i][0]))
        else:
            ynew.append(0)
    x = xnew
    y = ynew
    # Determine range of data.  Used to ensure constant colour map and
    # to set y scale of line plot.
    fmax = []
    fmin = []
    xmax = []
    dummymax = []
    dummymin = []
    clevels = []

    for i in range(0,Nvar):

        dummymax.append([])
        dummymin.append([])
        for j in range(0,Nlines[i]):
            dummymax[i].append(max(vars[i][j]))
            dummymin[i].append(min(vars[i][j]))

        fmax.append(max(dummymax[i]))
        fmin.append(min(dummymin[i]))

        if(symmetric_colors):
            absmax = max(abs(fmax[i]),abs(fmin[i]))
            fmax[i] = absmax
            fmin[i] = -absmax

        for j in range(0,Nlines[i]):
            dummymax[i][j] = max(x[i][j])
        xmax.append(max(dummymax[i]))


        if not (global_colors):
            clevels.append(linspace(fmin[i], fmax[i], Ncolors))
    if(global_colors):
        fmaxglobal = max(fmax)
        fminglobal = min(fmin)
        for i in range(0,Nvar):
            fmax[i]  = fmaxglobal
            fmin[i]  = fminglobal
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
        raise ValueError('too many variables...')


    fig = plt.figure(figsize=(w,h))
    title = fig.suptitle(r' ', fontsize=14  )

    # Initiate all list variables required for plotting here
    ax = []
    lines = []
    plots = []
    cbars = []
    xstride = []
    ystride = []
    r = []
    theta = []

    # Initiate figure frame
    for i in range(0,Nvar):
        lines.append([])
        if (lineplot[i] == 1):
            ax.append(fig.add_subplot(row,col,i+1))
            ax[i].set_xlim((0,xmax[i]))
            ax[i].set_ylim((fmin[i], fmax[i]))
            for j in range(0,Nlines[i]):
                lines[i].append(ax[i].plot([],[],lw=2, label = legendlabels[i][j])[0])
                #Need the [0] to 'unpack' the line object from tuple.  Alternatively:
                #lines[i], = lines[i]
            ax[i].set_xlabel(r'x')
            ax[i].set_ylabel(titles[i])
            if (Nlines[i] != 1):
                legendneeded = 1
                for k in range(0,i):
                    if (Nlines[i] == Nlines[k]):
                        legendneeded = 0
                if (legendneeded == 1):
                    plt.axes(ax[i])
                    plt.legend(loc = 0)
            # Pad out unused list variables with zeros
            plots.append(0)
            cbars.append(0)
            xstride.append(0)
            ystride.append(0)
            r.append(0)
            theta.append(0)

        elif (contour[i] == 1):
            ax.append(fig.add_subplot(row,col,i+1))
            #ax[i].set_xlim((0,Nx[i][0]-1))
            #ax[i].set_ylim((0,Ny[i][0]-1))
            ax[i].set_xlim(min(x[i]),max(x[i]))
            ax[i].set_ylim(min(y[i]),max(y[i]))
            ax[i].set_xlabel(r'x')
            ax[i].set_ylabel(r'y')
            ax[i].set_title(titles[i])
            if hold_aspect:
                ax[i].set_aspect('equal')
            plots.append(ax[i].contourf(x[i][0],y[i],vars[i][0][0,:,:].T, Ncolors, lw=0, levels=clevels[i] ))
            plt.axes(ax[i])
            cbars.append(fig.colorbar(plots[i], format='%1.1e'))
            # Pad out unused list variables with zeros
            lines[i].append(0)
            xstride.append(0)
            ystride.append(0)
            r.append(0)
            theta.append(0)

        elif (surf[i] == 1):
            x[i][0],y[i] = meshgrid(x[i][0],y[i])
            if (Nx[i][0]<= 20):
                xstride.append(1)
            else:
                xstride.append(int(floor(Nx[i][0]/20)))
            if (Ny[i][0]<=20):
                ystride.append(1)
            else:
                ystride.append(int(floor(Ny[i][0]/20)))
            ax.append(fig.add_subplot(row,col,i+1, projection='3d'))
            plots.append(ax[i].plot_wireframe(x[i][0], y[i], vars[i][0][0,:,:].T, rstride=ystride[i], cstride=xstride[i]))
            title = fig.suptitle(r'', fontsize=14 )
            ax[i].set_xlabel(r'x')
            ax[i].set_ylabel(r'y')
            ax[i].set_zlabel(titles[i])
            # Pad out unused list variables with zeros
            lines[i].append(0)
            cbars.append(0)
            r.append(0)
            theta.append(0)

        elif (polar[i] == 1):
            r.append(linspace(1,Nx[i][0], Nx[i][0]))
            theta.append(linspace(0,2*pi, Ny[i][0]))
            r[i],theta[i] = meshgrid(r[i], theta[i])
            ax.append(fig.add_subplot(row,col,i+1, projection='polar'))
            plots.append(ax[i].contourf(theta[i], r[i], vars[i][0][0,:,:].T, levels=clevels[i]))
            plt.axes(ax[i])
            cbars.append(fig.colorbar(plots[i], format='%1.1e'))
            ax[i].set_rmax(Nx[i][0]-1)
            ax[i].set_title(titles[i])
            # Pad out unused list variables with zeros
            lines[i].append(0)
            xstride.append(0)
            ystride.append(0)



    def onClick(event):
        global pause
        pause ^= True


    def control():
        global j, pause
        if j == Nframes-1 : j = -1
        if not pause:
            j=j+1

        return j


    # Animation function
    def animate(i):
        j=control()

        index = j*intv

        for j in range(0,Nvar):
                if (lineplot[j] == 1):
                    for k in range(0,Nlines[j]):
                        lines[j][k].set_data(x[j][k], vars[j][k][index,:])
                elif (contour[j] == 1):
                    plots[j] = ax[j].contourf(x[j][0],y[j],vars[j][0][index,:,:].T, Ncolors, lw=0, levels=clevels[j])
                elif (surf[j] == 1):
                    ax[j] = fig.add_subplot(row,col,j+1, projection='3d')
                    plots[j] = ax[j].plot_wireframe(x[j][0], y[j], vars[j][0][index,:,:].T, rstride=ystride[j], cstride=xstride[j])
                    ax[j].set_zlim(fmin[j],fmax[j])
                    ax[j].set_xlabel(r'x')
                    ax[j].set_ylabel(r'y')
                    ax[j].set_title(titles[j])
                elif (polar[j] == 1):
                    plots[j] = ax[j].contourf(theta[j], r[j], vars[j][0][index,:,:].T, levels=clevels[j])
                    ax[j].set_rmax(Nx[j][0]-1)

        if (tslice == 0):
            title.set_text('t = %1.2e' % t[index])
        else:
            title.set_text('t = %i' % index)
        return plots

    def init():
        global j, pause
        j=-2
        pause = False
        return animate(0)






    # Call Animation function

    fig.canvas.mpl_connect('button_press_event', onClick)
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=Nframes)

    # Save movie with given name
    if ((isinstance(movie,str)==1)):
        try:
            anim.save(movie+'.mp4',writer = FFwriter, fps=30, extra_args=['-vcodec', 'libx264'])
        except Exception:
            print("Save failed: Check ffmpeg path")

    # Save movie with default name
    if ((isinstance(movie,str)==0)):
        if (movie != 0):
            try:
                anim.save('animation.mp4',writer = FFwriter, fps=28, extra_args=['-vcodec', 'libx264'])
            except Exception:
                print("Save failed: Check ffmpeg path")

    # Show animation
    if (movie == 0):
        plt.show()


"""
To do list
1. Speed up animations ????
2. Look at theta in polar plots - perioidic?!?
3. Log axes, colorbars
4. Figureplot
"""
