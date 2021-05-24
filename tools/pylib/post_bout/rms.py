from __future__ import division
from builtins import range
from past.utils import old_div

###
# rms(f) : compute growth rate vs. time based on rms of bout variable f for all grid points
# plot_rms (x,y): plots the graph growth_rate vs. time for grid point x,y
###


import numpy as np
from pylab import plot, show, xlabel, ylabel, tight_layout
from boutdata.collect import collect


def rms(f):

    nt = f.shape[0]

    ns = f.shape[1]
    ne = f.shape[2]
    nz = f.shape[3]

    ar = np.zeros([nz])

    rms = np.zeros([nt, ns, ne])

    for i in range(nt):
        for j in range(ns):
            for k in range(ne):
                ar = f[i, j, k, :]
                valav = np.sum(ar)
                tot = np.sum(old_div(np.power(ar - valav, 2), nz))
                rms[i, j, k] = np.sqrt(tot)
    return rms


def plot_rms(x, y):
    s = plot(np.gradient(np.log(rmsp[:, x, y])))
    ylabel("$\gamma / \omega_A$", fontsize=25)
    xlabel("Time$(\\tau_A)$", fontsize=25)
    tight_layout()
    return s


# test
if __name__ == "__main__":
    path = "../../../examples/elm-pb/data"

    data = collect("P", path=path)

    rmsp = rms(data)

    plot_rms(34, 32)
    tight_layout()
    show()
