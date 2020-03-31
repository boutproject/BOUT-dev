from . import fieldtracer

import numpy as np

import warnings

try:
    import matplotlib.animation as anim
    import matplotlib.pyplot as plt
    plotting_available = True
except ImportError:
    warnings.warn("Couldn't import matplotlib, plotting not available.")
    plotting_available = False


def plot_poincare(magnetic_field, xpos, zpos, yperiod, nplot=3, y_slices=None, revs=40, nover=20,
                  interactive=False):
    """Plot a Poincare graph of the field lines.

    Parameters
    ----------
    magnetic_field : :py:obj:`zoidberg.field.MagneticField`
        Magnetic field object
    xpos, zpos : array_like
        Starting X, Z locations
    yperiod : float
        Length of period in y domain
    nplot : int, optional
        Number of equally spaced y-slices to plot
    y_slices : list of int, optional
        List of y-slices to plot; overrides nplot
    revs : int, optional
        Number of revolutions (times around phi)
    interactive : bool, optional
        If True, plots can be interacted with via the mouse:
        - Left-click on the plot to trace a field-line from that point
        - Right-click to add an additional trace
        - Middle-click to clear added traces

    Returns
    -------
    fig, ax
        The matplotlib figure and axis used
    """

    if not plotting_available:
        warnings.warning("matplotlib not available, unable to plot")
        return

    # Get Poincare plot
    result, y_slices = fieldtracer.trace_poincare(magnetic_field, xpos, zpos, yperiod,
                                                  nplot=nplot, y_slices=y_slices, revs=revs,
                                                  nover=nover)
    nplot = len(y_slices)
    
    #######################################################
    # Plotting

    colours = ["k", "b", "r", "g", "c", "m"]
    if nplot > len(colours):
        colours += colours * np.floor(nplot/len(colours))
        
    fig, ax = plt.subplots(1,1)
    
    for index in range(nplot):
        r = result[:,index,..., 0]
        z = result[:,index,..., 1]
        style = {
            'marker'    : '.',
            'color'     : colours[index],
            'linestyle' : 'None',
            }
        ax.plot(r, z, **style)

    ax.set_xlabel("Radius [m]", fontsize=20)
    ax.set_ylabel("Height [m]", fontsize=20)
    ax.tick_params(axis='both', labelsize=15)
    
    for phi, colour in zip(y_slices, colours):
        ax.plot([], [], color=colour, label=r'$Y = {0:.2f}$'.format(phi))

    ax.legend()

    overplot, = ax.plot([], [], 'ok')

    def onclick(event):
        # Check if user clicked inside the plot
        if event.xdata is None:
            return

        # Middle mouse
        if event.button == 2:
            x_ = []
            z_ = []
        else:
            result = field_tracer.follow_field_lines(event.xdata, event.ydata, y_values)
            # Right mouse
            if event.button == 3:
                x_, z_ = overplot.get_data()
                x_ = np.append(x_, result[:,0])
                z_ = np.append(z_, result[:,1])
            # Left mouse
            else:
                x_ = result[:,0]
                z_ = result[:,1]

        overplot.set_data(x_, z_)
        ax.relim()
        ax.autoscale_view()
        fig.canvas.draw()

    if interactive:
        field_tracer = fieldtracer.FieldTracer(magnetic_field)
        
        revs = int(revs)    
        y_values = y_slices[:]
        for n in np.arange(1, revs):
            y_values = np.append(y_values, n*yperiod + y_values[:nplot])
        
        fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()

    return fig, ax


def plot_3d_field_line(magnetic_field, xpos, zpos, yperiod, cycles=20, y_res=50):
    """Make a 3D plot of field lines

    Parameters
    ----------
    magnetic_field : :py:obj:`zoidberg.field.MagneticField`
        Magnetic field object
    xpos, zpos : array_like
        Starting X, Z locations
    yperiod : float
        Length of period in y domain
    cycles : int, optional
        Number of times to go round in y
    y_res : int, optional
        Number of points in y in each cycle

    Returns
    -------
    fig, ax
        The matplotlib figure and axis used

    """

    if not plotting_available:
        warnings.warning("matplotlib not available, unable to plot")
        return

    yperiod = float(yperiod)

    # Go round toroidally cycles times
    phivals_hires = np.linspace(0, cycles*yperiod, num=y_res*cycles, endpoint=False)
    
    xpos = np.asfarray(xpos)
    zpos = np.asfarray(zpos)

    field_tracer = fieldtracer.FieldTracer(magnetic_field)
    result_hires = field_tracer.follow_field_lines(xpos, zpos, phivals_hires)

    # Get phivals_hires into [0,yperiod]
    phivals_hires_mod = np.remainder(phivals_hires, yperiod)
    # There are cycles sets of field lines y_res points long each
    # and we also need to transpose for reasons
    phivals_hires_mod = phivals_hires_mod.reshape( (cycles, y_res) ).T
    # Same for the result, but only swap first and second indices
    result_hires_mod = result_hires.reshape( (cycles,y_res,2) ).transpose(1,0,2)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for n in range(cycles):
        ax.plot(result_hires_mod[:,n,0], result_hires_mod[:,n,1], phivals_hires_mod[:,n])

    plt.show()

    return fig, ax


def plot_streamlines(grid, magnetic_field, y_slice=0, width=None, **kwargs):
    """Plot streamlines of the magnetic field in the poloidal plane

    Parameters
    ----------
    grid : :py:obj:`zoidberg.grid.Grid`
        Grid generated by Zoidberg
    magnetic_field : :py:obj:`zoidberg.field.MagneticField`
        Zoidberg magnetic field object
    y_slice : int, optional
        y-index to plot streamlines at
    width : float, optional
        If not None, line widths are proportional to the magnitude of
        the `magnetic_field` times `width`

    Returns
    -------
    fig, ax
        The matplotlib figure and axis used

    """

    if not plotting_available:
        warnings.warning("matplotlib not available, unable to plot")
        return

    fig, ax = plt.subplots(1,1)
    full_slice = np.s_[:, y_slice, :]

    if width is not None:
        # Get the B field magnitude in the poloidal plane
        bxz_mag = np.sqrt(magnetic_field.b_mag**2 - magnetic_field.by**2)
        linewidth = width*(bxz_mag[full_slice] / bxz_mag.max()).T
    else:
        linewidth = 1

    ax.streamplot(grid.xarray, grid.zarray,
                  magnetic_field.bx[full_slice].T,
                  magnetic_field.bz[full_slice].T,
                  linewidth=linewidth,
                  **kwargs)

    ax.set_xlabel("Radius [m]", fontsize=20)
    ax.set_ylabel("Height [m]", fontsize=20)
    ax.tick_params(axis='both', labelsize=15)

    plt.show()

    return fig, ax


class AnimateVectorField(object):
    """Very basic/experimental class for animating vector fields

    Transpose U, V to have dimension to animate at index 0,
    e.g. to animate along y, pass:

    >>> AnimateVectorField(X, Y, U.transpose((1,0,2)), V.transpose((1,0,2)))

    Parameters
    ----------
    X, Y : array_like
        The X, Y coordinates
    U, V : ndarray
        Vector components in X, Y respectively

    Examples
    --------

    >>> anim = AnimateVectorField(X, Y, U, V)
    >>> anim.animate()
    """

    def __init__(self, X, Y, U, V):
        self.X = X
        self.Y = Y
        self.U = U
        self.V = V

        self.frames = U.shape[0]

        self.fig, self.ax = plt.subplots(1,1)
        self.quiv_plot = self.ax.quiver(X, Y, U[0,...].T, V[0,...].T, pivot='mid', angles='xy')

    def __update_quiver(self, num):
        self.quiv_plot.set_UVC(self.U[num,...].T, self.V[num,...].T)
        return self.quiv_plot,

    def animate(self):
        self.animation = anim.FuncAnimation(self.fig, self.__update_quiver, frames=self.frames,
                                            interval=100, blit=False)
        plt.show()
        plt.draw()


def plot_forward_map(grid, maps, yslice=0):
    """Plots the forward map from `yslice` to `yslice` + 1

    Parameters
    ----------
    grid : :py:obj:`zoidberg.grid.Grid`
        Grid generated by Zoidberg
    maps : dict
        Dictionary containing the forward FCI maps
    y_slice : int, optional
        Originating y-index to plot map from

    TODO
    ----
    - Take optional axis
    - Take optional show argument

    """

    if not plotting_available:
        warnings.warning("matplotlib not available, unable to plot")
        return
    
    nx, ny, nz = grid.shape
    
    pol, ycoord = grid.getPoloidalGrid(yslice)
    pol_next, ycoord_next = grid.getPoloidalGrid(yslice+1)
    
    # Plot the points on yslice+1 as 'bx'
    # Note: ravel used here only so multiple labels are not created
    plt.plot(np.ravel(pol_next.R), np.ravel(pol_next.Z), 'bx', label="Grid points on slice {0}".format(yslice+1))
    
    # Plot the forward map from yslice to yslice+1 as red 'o'
    forward_R = maps['forward_R'][:,yslice,:]
    forward_Z = maps['forward_Z'][:,yslice,:]
    plt.plot(np.ravel(forward_R), np.ravel(forward_Z), 'ro', label="Forward map from slice {0}".format(yslice))
    
    # Mark the points which hit the inner boundary
    # These are marked with a negative x index
    in_boundary = maps['forward_xt_prime'][:,yslice,:] < 0.0
    plt.plot( np.ravel(forward_R[ in_boundary ]), np.ravel(forward_Z[ in_boundary ]), 'ko', label="Inner boundary points")
    
    # Outer boundary marked with x index nx
    out_boundary = maps['forward_xt_prime'][:,yslice,:] > nx-0.5
    plt.plot( np.ravel(forward_R[ out_boundary ]), np.ravel(forward_Z[ out_boundary ]), 'bo', label="Outer boundary points")
    
    
    plt.legend()
    
    plt.show()


def plot_backward_map(grid, maps, yslice=0):
    """Plots the backward map from yslice to yslice-1

    Parameters
    ----------
    grid : 'zoidberg.grid.Grid`
        Grid generated by Zoidberg
    maps : dict
        Dictionary containing the backward FCI maps
    y_slice : int, optional
        Originating y-index to plot map from

    TODO
    ----
    - Take optional axis
    - Take optional show argument

    """

    if not plotting_available:
        warnings.warning("matplotlib not available, unable to plot")
        return

    nx, ny, nz = grid.shape
    
    pol, ycoord = grid.getPoloidalGrid(yslice)
    pol_last, ycoord_last = grid.getPoloidalGrid(yslice-1)
    
    # Plot the points on yslice-1 as 'bx'
    # Note: ravel used here only so multiple labels are not created
    plt.plot(np.ravel(pol_last.R), np.ravel(pol_last.Z), 'bx', label="Grid points on slice {0}".format(yslice-1))
    
    # Plot the backward map from yslice to yslice-1 as red 'o'
    backward_R = maps['backward_R'][:,yslice,:]
    backward_Z = maps['backward_Z'][:,yslice,:]
    plt.plot(np.ravel(backward_R), np.ravel(backward_Z), 'ro', label="Backward map from slice {0}".format(yslice))
    
    # Mark the points which hit the inner boundary
    # These are marked with a negative x index
    in_boundary = maps['backward_xt_prime'][:,yslice,:] < 0.0
    plt.plot( np.ravel(backward_R[ in_boundary ]), np.ravel(backward_Z[ in_boundary ]), 'ko', label="Inner boundary points")
    
    # Outer boundary marked with x index nx
    out_boundary = maps['backward_xt_prime'][:,yslice,:] > nx-0.5
    plt.plot( np.ravel(backward_R[ out_boundary ]), np.ravel(backward_Z[ out_boundary ]), 'bo', label="Outer boundary points")
    
    
    plt.legend()
    
    plt.show()
