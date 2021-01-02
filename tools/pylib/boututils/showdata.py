"""
Visualisation and animation routines

Written by Luke Easy
le590@york.ac.uk
Last Updated 19/3/2015
Additional functionality by George Breyiannis 26/12/2014

"""
from __future__ import print_function
from __future__ import division
from builtins import str, chr, range

from matplotlib import pyplot as plt
from matplotlib import animation
from numpy import linspace, meshgrid, array, min, max, abs, floor, pi, isclose
from boutdata.collect import collect
from boututils.boutwarnings import alwayswarn

####################################################################
# Specify manually ffmpeg path
#plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

FFwriter = animation.FFMpegWriter()
####################################################################


###################
#https://stackoverflow.com/questions/16732379/stop-start-pause-in-python-matplotlib-animation
#
j=-2
pause = False
###################


def showdata(vars, titles=[], legendlabels=[], surf=[], polar=[], tslice=0, t_array=None,
             movie=0, fps=28, dpi=200, intv=1, Ncolors=25, x=[], y=[],
             global_colors=False, symmetric_colors=False, hold_aspect=False,
             cmap=None, clear_between_frames=None, return_animation=False, window_title=""):
    """A Function to animate time dependent data from BOUT++

    To animate multiple variables on different axes:

    >>> showdata([var1, var2, var3])

    To animate more than one line on a single axes:

    >>> showdata([[var1, var2, var3]])

    The default graph types are:
    2D (time + 1 spatial dimension) arrays = animated line plot
    3D (time + 2 spatial dimensions) arrays = animated contour plot.

    To use surface or polar plots:

    >>> showdata(var, surf=1)
    >>> showdata(var, polar=1)

    Can plot different graph types on different axes.  Default graph types will
    be used depending on the dimensions of the input arrays.  To specify
    polar/surface plots on different axes:

    >>> showdata([var1, var2], surf=[1, 0], polar=[0, 1])

    Movies require FFmpeg (for .mp4) and/or ImageMagick (for .gif) to be
    installed.  The 'movie' option can be set to 1 (which will produce an mp4
    called 'animation.mp4'), to a name with no extension (which will produce an
    mp4 called '<name>.mp4')

    The `tslice` variable is used to control the time value that is printed on
    each frame of the animation. If the input data matches the time values
    found within BOUT++'s dmp data files, then these time values will be used.
    Otherwise, an integer counter is used.

    The `cmap` variable (if specified) will set the colormap used in the plot
    cmap must be a matplotlib colormap instance, or the name of a registered
    matplotlib colormap

    During animation click once to stop in the current frame. Click again to
    continue.

    Parameters
    ----------
    vars : array_like or list of array_like
        Variable or list of variables to plot
    titles : str or list of str, optional
        Title or list of titles for each axis
    legendlabels : str or list of str, optional
        Legend or list of legends for each variable
    surf : list of int
        Which axes to plot as a surface plot
    polar : list of int
        Which axes to plot as a polar plot
    tslice : list of int
        Use these time values from a dump file (see above)
    t_array : array
        Pass in t_array using this argument to use the simulation time in plot
        titles. Otherwise, just use the t-index.
    movie : int
        If 1, save the animation to file
    fps : int
        Number of frames per second to use when saving animation
    dpi : int
        Dots per inch to use when saving animation
    intv : int
        ???
    Ncolors : int
        Number of levels in contour plots
    x, y : array_like, list of array_like
        X, Y coordinates
    global_colors : bool
        If "vars" is a list the colorlevels are determined from the
        maximum of the maxima and and the minimum of the minima in all
        fields in vars
    symmetric_colors : bool
        Colour levels are symmetric
    hold_aspect : bool
        Use equal aspect ratio in plots
    cmap : colormap
        A matplotlib colormap instance to use
    clear_between_frames : bool, optional
        - Default (None) - all plots except line plots will clear between frames
        - True - all plots will clear between frames
        - False - no plots will clear between frames
    return_animation : bool
        Return the matplotlib animation instance
    window_title : str
        Give a title for the animation window

    TODO
    ----
    - Replace empty lists in signature with None
    - Use bools in sensible places
    - Put massive list of arguments in kwargs
    - Speed up animations ????
    - Look at theta in polar plots - periodic?!?
    - Log axes, colorbars
    - Figureplot

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
                alwayswarn("The legendlabels list does not contain a sublist for each variable, but its length matches the number of lines on each plot. Will apply labels to each plot")
                legendlabelsdummy = []
                for i in range(0, Nvar):
                    legendlabelsdummy.append([])
                    for j in range(0,Nlines[i]):
                        legendlabelsdummy[i].append(legendlabels[j])
                legendlabels = legendlabelsdummy
            else:
                alwayswarn("The legendlabels list does not contain a sublist for each variable, and it's length does not match the number of lines on each plot. Will default apply labels to each plot")
                legendlabels = []
                for i in range(0,Nvar):
                    legendlabels.append([])
                    for j in range(0,Nlines[i]):
                        legendlabels[i].append(chr(97+j))
        else:
            if (Nlines[0] == len(legendlabels)):
                legendlabels = [legendlabels]
    elif len(legendlabels) != Nvar:
        alwayswarn("The length of the legendlabels list does not match the length of the vars list, will continue with default values")
        legendlabels = []
        for i in range(0,Nvar):
            legendlabels.append([])
            for j in range(0,Nlines[i]):
                legendlabels[i].append(chr(97+j))
    else:
        for i in range(0,Nvar):
            if isinstance(legendlabels[i], list):
                if len(legendlabels[i]) != Nlines[i]:
                    alwayswarn('The length of the legendlabel (sub)list for each plot does not match the number of datasets for each plot. Will continue with default values')
                legendlabels[i] = []
                for j in range(0,Nlines[i]):
                    legendlabels[i].append(chr(97+j))
            else:
                legendlabels[i] = [legendlabels[i]]
            if len(legendlabels[i]) != Nlines[i]:
                alwayswarn('The length of the legendlabel (sub)list for each plot does not match the number of datasets for each plot.  Will continue with default values')
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
            alwayswarn('Length of surf list does not match number of variables.  Will default to no polar plots')
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
            alwayswarn('Length of polar list does not match number of variables.  Will default to no polar plots')
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
                alwayswarn('Data must be  3 dimensional (time, r, theta) for polar plots.  Will plot lineplot instead')

            if ((Ndims[i][j] == 2) & (surf[i] != 0)):
                alwayswarn('Data must be  3 dimensional (time, x, y) for surface plots.  Will plot lineplot instead')

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
                alwayswarn('Cannot do polar and surface plots at the same time.  Default to contour plot')
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

    # Obtain number of frames
    Nframes = int(Nt[0][0]/intv)

    # Generate grids for plotting
    # Try to use provided grids where possible
    # If x and/or y are not lists, apply to all variables
    if not isinstance(x, (list,tuple)):
        x = [x]*Nvar # Make list of x with length Nvar
    if not isinstance(y, (list,tuple)):
        y = [y]*Nvar # Make list of x with length Nvar
    xnew = []
    ynew = []
    for i in range(0,Nvar):
        xnew.append([])
        try:
            xnew[i].append(x[i])
            if not (x[i].shape==(Nx[i][0],) or x[i].shape==(Nx[i][0],Ny[i][0]) or x[i].shape==(Nt[i][0],Nx[i][0],Ny[i],[0])):
                raise ValueError("For variable number "+str(i)+", "+titles[i]+", the shape of x is not compatible with the shape of the variable. Shape of x should be (Nx), (Nx,Ny) or (Nt,Nx,Ny).")
        except:
            for j in range(0, Nlines[i]):
                xnew[i].append(linspace(0,Nx[i][j]-1, Nx[i][j]))

        #x.append(linspace(0,Nx[i][0]-1, Nx[i][0]))

        if (Ndims[i][0] == 3):
            try:
                ynew.append(y[i])
                if not (y[i].shape==(Ny[i][0],) or y[i].shape==(Nx[i][0],Ny[i][0]) or y[i].shape==(Nt[i][0],Nx[i][0],Ny[i],[0])):
                    raise ValueError("For variable number "+str(i)+", "+titles[i]+", the shape of y is not compatible with the shape of the variable. Shape of y should be (Ny), (Nx,Ny) or (Nt,Nx,Ny).")
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
            absmax =max(abs(array(fmax[i], fmin[i])))
            fmax[i] = absmax
            fmin[i] = -absmax

        for j in range(0,Nlines[i]):
            dummymax[i][j] = max(x[i][j])
        xmax.append(max(dummymax[i]))


        if not (global_colors):
            if isclose(fmin[i], fmax[i]):
                # add/subtract very small constant in case fmin=fmax=0
                thiscontourmin = fmin[i] - 3.e-15*abs(fmin[i]) - 1.e-36
                thiscontourmax = fmax[i] + 3.e-15*abs(fmax[i]) + 1.e-36
                alwayswarn("Contour levels too close, adding padding to colorbar range")
                clevels.append(linspace(thiscontourmin, thiscontourmax, Ncolors))
            else:
                clevels.append(linspace(fmin[i], fmax[i], Ncolors))

    if(global_colors):
        fmaxglobal = max(fmax)
        fminglobal = min(fmin)
        if isclose(fminglobal, fmaxglobal):
            fminglobal = fminglobal - 3.e-15*abs(fminglobal) - 1.e-36
            fmaxglobal = fmaxglobal + 3.e-15*abs(fmaxglobal) + 1.e-36
        for i in range(0,Nvar):
            clevels.append(linspace(fminglobal, fmaxglobal, Ncolors))

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


    fig = plt.figure(window_title, figsize=(w,h))
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
            thisx = x[i][0]
            if len(thisx.shape) == 3:
                thisx = thisx[0]
            thisy = y[i]
            if len(thisy.shape) == 3:
                thisy = thisy[0]
            plots.append(ax[i].contourf(thisx.T,thisy.T,vars[i][0][0,:,:].T, Ncolors, cmap=cmap, lw=0, levels=clevels[i] ))
            plt.axes(ax[i])
            cbars.append(fig.colorbar(plots[i], format='%1.1e'))
            # Pad out unused list variables with zeros
            lines[i].append(0)
            xstride.append(0)
            ystride.append(0)
            r.append(0)
            theta.append(0)

        elif (surf[i] == 1):
            if (len(x[i][0].shape)==1 and len(y[i].shape)==1):
                # plot_wireframe() requires 2d arrays for x and y coordinates
                x[i][0],y[i] = meshgrid(x[i][0],y[i])
            thisx = x[i][0]
            if len(thisx.shape) == 3:
                thisx = thisx[0]
            thisy = y[i]
            if len(thisy.shape) == 3:
                thisy = thisy[0]
            if (Nx[i][0]<= 20):
                xstride.append(1)
            else:
                xstride.append(int(floor(Nx[i][0]/20)))
            if (Ny[i][0]<=20):
                ystride.append(1)
            else:
                ystride.append(int(floor(Ny[i][0]/20)))
            ax.append(fig.add_subplot(row,col,i+1, projection='3d'))
            plots.append(ax[i].plot_wireframe(thisx, thisy, vars[i][0][0,:,:].T, rstride=ystride[i], cstride=xstride[i]))
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
            plots.append(ax[i].contourf(theta[i], r[i], vars[i][0][0,:,:].T, cmap=cmap, levels=clevels[i]))
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
            #Default to clearing axis between frames on all plots except line plots
            if (clear_between_frames is None and lineplot[j] != 1 ) or clear_between_frames is True:
                ax[j].cla() #Clear axis between frames so that masked arrays can be plotted
            if (lineplot[j] == 1):
                for k in range(0,Nlines[j]):
                    lines[j][k].set_data(x[j][k], vars[j][k][index,:])
            elif (contour[j] == 1):
                thisx = x[j][0]
                if len(thisx.shape) == 3:
                    thisx = thisx[index]
                thisy = y[j]
                if len(thisy.shape) == 3:
                    thisy = thisy[index]
                plots[j] = ax[j].contourf(x[j][0].T,y[j].T,vars[j][0][index,:,:].T, Ncolors, cmap=cmap, lw=0, levels=clevels[j])
                ax[j].set_xlabel(r'x')
                ax[j].set_ylabel(r'y')
                ax[j].set_title(titles[j])
            elif (surf[j] == 1):
                thisx = x[j][0]
                if len(thisx.shape) == 3:
                    thisx = thisx[index]
                thisy = y[j][0]
                if len(thisy.shape) == 3:
                    thisy = thisy[index]
                ax[j] = fig.add_subplot(row,col,j+1, projection='3d')
                plots[j] = ax[j].plot_wireframe(thisx, thisy, vars[j][0][index,:,:].T, rstride=ystride[j], cstride=xstride[j])
                ax[j].set_zlim(fmin[j],fmax[j])
                ax[j].set_xlabel(r'x')
                ax[j].set_ylabel(r'y')
                ax[j].set_title(titles[j])
            elif (polar[j] == 1):
                plots[j] = ax[j].contourf(theta[j], r[j], vars[j][0][index,:,:].T,cmap=cmap, levels=clevels[j])
                ax[j].set_rmax(Nx[j][0]-1)
                ax[j].set_title(titles[j])

        if t_array is not None:
            title.set_text('t = %1.2e' % t_array[index])
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

    #If movie is not passed as a string assign the default filename
    if (movie==1):
        movie='animation.mp4'

    # Save movie with given or default name
    if ((isinstance(movie,str)==1)):
        movietype = movie.split('.')[-1]
        if movietype == 'mp4':
            try:
                anim.save(movie,writer = FFwriter, fps=fps, dpi=dpi, extra_args=['-vcodec', 'libx264'])
            except Exception:
            #Try specifying writer by string if ffmpeg not found
                try:
                    anim.save(movie,writer = 'ffmpeg', fps=fps, dpi=dpi, extra_args=['-vcodec', 'libx264'])
                except Exception:
                     print('Save failed: Check ffmpeg path')
                     raise
        elif movietype == 'gif':
            anim.save(movie,writer = 'imagemagick', fps=fps, dpi=dpi)
        else:
            raise ValueError("Unrecognized file type for movie. Supported types are .mp4 and .gif")

    # Show animation if not saved or returned, otherwise close the plot
    if (movie==0 and return_animation == 0):
        plt.show()
    else:
        plt.close()
    # Return animation object
    if(return_animation == 1):
        return(anim)
