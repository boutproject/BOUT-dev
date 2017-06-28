from . import fieldtracer

import numpy as np
import matplotlib.animation as anim
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint

def plot_poincare(magnetic_field, xpos, zpos, yperiod, nplot=3, phi_slices=None, revs=100,
                  interactive=False):
    """Plot a Poincare graph of the field lines.

    Inputs
    ------
    magnetic_field   Magnetic field object
    
    xpos             Starting X location. Can be scalar or list/array
    zpos             Starting Z location. Can be scalar or list/array
    
    yperiod          Length of period in y domain

    nplot            Number of equally spaced phi-slices to plot [3]
    
    phi_slices       List of phi-slices to plot; overrides nplot
    
    revs             Number of revolutions (times around phi) [40]
    
    interactive      Left-click on the plot to trace a field-line from that point
                     Right-click to add an additional trace
                     Middle-click to clear added traces
    """

    colours = ["k", "b", "r", "g", "c", "m"]

    # Check arguments are ok
    if nplot is None and phi_slices is None:
        raise ValueError("nplot and phi_slices cannot both be None")
    if phi_slices is not None:
        phi_slices = np.asfarray(phi_slices)
        
        if np.amin(phi_slices) < 0.0 or np.amax(phi_slices) > yperiod:
            raise ValueError("phi_slices must all be between 0.0 and yperiod ({yperiod})"
                             .format(yperiod=yperiod))
        # Make sure phi_slices is monotonically increasing
        phi_slices.sort()
        # If phi_slices is given, then nplot is the number of slices
        nplot = len(phi_slices)
    else:
        # nplot equally spaced phi slices
        nplot = int(nplot)
        phi_slices = np.linspace(0, yperiod, nplot, endpoint=False)
        
    if nplot > len(colours):
        colours += colours * np.floor(nplot/len(colours))

    ########################################################
    # Extend the domain from [0,yperiod] to [0,revs*yperiod]

    
    phi_values = phi_slices[:]
    for n in np.arange(1, revs):
        phi_values = np.append(phi_values, n*yperiod + phi_values[:nplot])
        
    nover = 20 # Over-sample
    phi_values_over = np.zeros( ( nplot * revs * nover - (nover-1)) )
    phi_values_over[::nover] = phi_values
    for i in range(1,nover):
        phi_values_over[i::nover] = (float(i)/float(nover))*phi_values[1:] + (float(nover-i)/float(nover))*phi_values[:-1]
    
    #######################################################
    # Plotting

    # Starting location
    xpos = np.asfarray(xpos)
    zpos = np.asfarray(zpos)

    field_tracer = fieldtracer.FieldTracer(magnetic_field)
    result = field_tracer.follow_field_lines(xpos, zpos, phi_values_over)

    result = result[::nover,...] # Remove unneeded points

    fig, ax = plt.subplots(1,1)
    
    for index in range(nplot):
        r = result[index::nplot,..., 0]
        z = result[index::nplot,..., 1]
        style = {
            'marker'    : '.',
            'color'     : colours[index],
            'linestyle' : 'None',
            }
        ax.plot(r, z, **style)

    ax.set_xlabel("Radius [m]", fontsize=20)
    ax.set_ylabel("Height [m]", fontsize=20)
    ax.tick_params(axis='both', labelsize=15)

    for phi, colour in zip(phi_slices, colours):
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
            result = field_tracer.follow_field_lines(event.xdata, event.ydata, phi_values)
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
