from . import fieldtracer

import numpy as np
import matplotlib.animation as anim
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint

def plot_poincare(grid, magnetic_field, nplot=3):

    sym=[".k", ".b", ".r", ".g"]
    phivals = np.arange(100 * nplot)*2.*np.pi/(nplot)
    # phivals = np.linspace(0, 100*2*np.pi, 1000)

    fig, ax = plt.subplots(1,1)

    xpos = grid.xcentre + np.linspace(0, 0.5*np.max(grid.xarray), 10)
    zpos = grid.zcentre + np.zeros(xpos.shape)

    position = np.column_stack((xpos, zpos)).flatten()

    # for xpos in grid.xcentre + np.linspace(0, 0.5*np.max(grid.xarray),10):
    #     pos = (xpos, grid.zcentre)
    result = odeint(magnetic_field.field_direction, position, phivals, args=(True,))
    result = result.reshape((-1,2))
    for p in range(nplot):
        ax.plot(result[p::nplot,0], result[p::nplot,1],sym[p])

    ax.set_xlabel("Radius [m]", fontsize=20)
    ax.set_ylabel("Height [m]", fontsize=20)
    ax.tick_params(axis='both', labelsize=15)

    for p in range(nplot):
        ax.plot([], [], sym[p], label=r'$Y = \left({0} * L/{1}\right)$'.format(p, nplot))

    ax.legend()

    plt.show()

    return fig, ax

def plot_3d_field_line(grid, magnetic_field, cycles=20):
    # Go round toroidally 20 times
    phivals_hires = np.linspace(0, cycles*2*np.pi, num=50*cycles)

    xpos = grid.xcentre + 0.5*np.max(grid.xarray)
    result_hires = odeint(magnetic_field.field_direction, (xpos, grid.zcentre), phivals_hires)
    # Get phivals_hires into [0,2pi]
    phivals_hires_mod = np.mod(phivals_hires, 2*np.pi)
    # There are 20 sets of field lines 50 points long each
    # and we also need to transpose for reasons
    phivals_hires_mod = phivals_hires_mod.reshape( (cycles, 50) ).T
    # Same for the result, but only swap first and second indices
    result_hires_mod = result_hires.reshape( (cycles,50,2) ).transpose(1,0,2)

    fig = plt.figure()

    ax = fig.gca(projection='3d')
    for n in xrange(cycles):
        ax.plot(result_hires_mod[:,n,0], result_hires_mod[:,n,1], phivals_hires_mod[:,n])

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
