from . import fieldtracer

import numpy as np
from scipy.integrate import odeint

import warnings

try:
    import matplotlib.animation as anim
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    plotting_available = True
except ImportError:
    warnings.warn("Couldn't import matplotlib, plotting not available.")
    plotting_available = False


def plot_poincare(magnetic_field, xpos, zpos, yperiod, nplot=3, y_slices=None, revs=40, nover=20,
                  interactive=False):
    """Plot a Poincare graph of the field lines.

    Inputs
    ------
    magnetic_field   Magnetic field object
    
    xpos             Starting X location. Can be scalar or list/array
    zpos             Starting Z location. Can be scalar or list/array
    
    yperiod          Length of period in y domain

    nplot            Number of equally spaced y-slices to plot
    
    y_slices         List of y-slices to plot; overrides nplot
    
    revs             Number of revolutions (times around phi)
    
    interactive      Left-click on the plot to trace a field-line from that point
                     Right-click to add an additional trace
                     Middle-click to clear added traces
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
        if event.button is 2:
            x_ = []
            z_ = []
        else:
            result = field_tracer.follow_field_lines(event.xdata, event.ydata, y_values)
            # Right mouse
            if event.button is 3:
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

    Inputs
    ------
    magnetic_field - Magnetic field object

    xpos             Starting X location. Can be scalar or list/array
    zpos             Starting Z location. Can be scalar or list/array
    
    yperiod          Length of period in y domain

    cycles         - Number of times to go round in y [20]
    y_res          - Number of points in y in each cycle [50]
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

    Inputs
    ------
    grid           - Grid generated by Zoidberg
    magnetic_field - Zoidberg magnetic field object
    y_slice        - y-index to plot streamlines at
    width          - If not None, line widths are proportional to the
                     magnitude of B times 'width' [ None ]
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
    """Transpose U, V to have dimension to animate at index 0,
    e.g. to animate along y, pass:
        U.transpose( (1,0,2) )
        V.transpose( (1,0,2) )
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
    """
    Plots the forward map from yslice to yslice+1
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
    """
    Plots the backward map from yslice to yslice-1
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
