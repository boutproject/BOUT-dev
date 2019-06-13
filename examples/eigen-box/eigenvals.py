#!/usr/bin/env python3

from boutdata.collect import collect

import matplotlib.pyplot as plt

from numpy import argmin, amax, amin, arange, ndarray


def plot_eigenvals(eigenvalues, eigenvectors=None):
    """Plot the eigenvalues and, optionally, the eigenvectors for a 1D simulation

    eigenvalues should be a 1D array where the odd indices are the
    real components, evens are the imaginary

    eigenvectors should be a 2D array whose first dimension is the
    same length as eigenvalues, with the same interpretation

    """

    if eigenvalues.shape[0] % 2 != 0:
        raise ValueError("Odd number of elements in eigenvalues")

    if eigenvectors is not None:
        # Check that the eigenvectors has the right size
        if len(eigenvectors.shape) != 2:
            raise ValueError("Expecting eigenvectors to be 2D")
        if eigenvectors.shape[0] != len(eigenvalues):
            raise ValueError(
                "First dimension of eigenvectors must match length of eigenvalues")

    # If no eigenvectors supplied, only plot eigenvalues, otherwise
    # eigenvalues and eigenvectors
    nrows = 1 if eigenvectors is None else 2
    fig, ax = plt.subplots(nrows=nrows)

    # If we've only made one plot, ax won't be a list
    if not isinstance(ax, (list, ndarray)):
        ax = [ax]

    eigs_r = eigenvalues[:-1:2]
    eigs_i = eigenvalues[1::2]

    range_r = amax(eigs_r) - amin(eigs_r)
    range_i = amax(eigs_i) - amin(eigs_i)

    ax[0].plot(eigs_r, eigs_i, 'x')
    ax[0].set_xlabel("Real component")
    ax[0].set_ylabel("Imaginary component")
    ax[0].set_title("Eigenvalue")

    overplot, = ax[0].plot([], [], 'ok')

    if eigenvectors is not None:
        # Add a eigenvectors plot
        vector_r, = ax[1].plot([], [], '-k', label="Real")
        vector_i, = ax[1].plot([], [], '-r', label="Imag")
        ax[1].legend(loc='upper right')
        ax[1].set_xlabel("X")
        ax[1].set_ylabel("Amplitude")
        ax[1].set_title("Eigenvector")
        plt.subplots_adjust(hspace=0.5)

    def onclick(event):
        # Check if user clicked inside the plot
        if event.xdata is None:
            return

        # Find closest eigenvectors point, but stretch axes so
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

        if eigenvectors is not None:
            # Update plots
            nx = eigenvectors.shape[1]
            vector_r.set_data(arange(nx), eigenvectors[2*ind, :])
            vector_i.set_data(arange(nx), eigenvectors[2*ind+1, :])
            ax[1].relim()
            ax[1].autoscale_view()

        fig.canvas.draw()

    fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()


if __name__ == "__main__":
    path = "data"
    eigenvalues = collect("t_array", path=path, info=False)
    eigenvectors = collect("f", xguards=False, path=path, info=False)
    plot_eigenvals(eigenvalues, eigenvectors[..., 0, 0])
