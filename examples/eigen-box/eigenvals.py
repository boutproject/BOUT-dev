
from boutdata.collect import collect

import matplotlib.pyplot as plt

from numpy import argmin, amax, amin, arange


def plot_eigenvals(eigs, data=None):
    """

    """
    fig = plt.figure()

    if data is None:
        # If no data supplied, only plot eigenvalues
        ax = fig.add_subplot(111)
    else:
        # If data, plot two figures
        ax = fig.add_subplot(211)

        # Check that the data has the right size
        if len(data.shape) != 2:
            raise ValueError("Expecting data to be 2D")
        if data.shape[0] != len(eigs):
            raise ValueError(
                "First dimension of data must match length of eigs")

    eigs_r = eigs[:-1:2]
    eigs_i = eigs[1::2]

    range_r = amax(eigs_r) - amin(eigs_r)
    range_i = amax(eigs_i) - amin(eigs_i)

    ax.plot(eigs_r, eigs_i, 'x')
    ax.set_xlabel("Real component")
    ax.set_ylabel("Imaginary component")

    overplot, = ax.plot([], [], 'ok')

    if data is not None:
        # Add a data plot
        ax2 = fig.add_subplot(212)
        vector_r, = ax2.plot([], [], '-k', label="Real")
        vector_i, = ax2.plot([], [], '-r', label="Imag")
        ax2.set_xlabel("X")
        ax2.set_ylabel("Eigenvector")

    def onclick(event):
        # Check if user clicked inside the plot
        if event.xdata is None:
            return

        # Find closest data point, but stretch axes so
        # real and imaginary components are weighted equally
        if(range_r == 0):
            dist = ((eigs_i - event.ydata)/range_i)**2
        elif(range_i == 0):
            dist = ((eigs_r - event.xdata)/range_r)**2
        else:
            dist = ((eigs_r - event.xdata)/range_r)**2 + \
                ((eigs_i - event.ydata)/range_i)**2

        ind = argmin(dist)

        # Update the highlight plot
        overplot.set_data([eigs_r[ind]], [eigs_i[ind]])

        print("Eigenvalue number: %d (%e,%e)" %
              (ind, eigs_r[ind], eigs_i[ind]))

        if data is not None:
            # Update plots
            nx = data.shape[1]
            vector_r.set_data(arange(nx), data[2*ind, :])
            vector_i.set_data(arange(nx), data[2*ind+1, :])
            ax2.relim()
            ax2.autoscale_view()

        fig.canvas.draw()

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()


if __name__ == "__main__":
    path = "data"
    eigs = collect("t_array", path=path)
    data = collect("f", path=path)
    plot_eigenvals(eigs, data=data[:, 2:-2, 0, 0])
