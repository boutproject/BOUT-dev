#!/usr/bin/env python
"""Classes and functions for setting style to plots."""

__authors__ = 'Michael Loeiten'
__email__   = 'mmag@fysik.dtu.dk'
__version__ = '0.301beta'
__date__    = '06.11.2014'

from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import numpy as np
import warnings
from _tkinter import TclError

# TODO: 2D solution and error plot
# Consider: Make it possible to change colors on the graphs


#{{{set_style
class set_style():
    """Class for setting the plotting style in plots from BOUT++"""

    # The following is shared by all instances
    title_size = 30
    plt.rc("font", size = 30)
    plt.rc("axes", labelsize = 25, titlesize = title_size)
    plt.rc("xtick", labelsize = 25)
    plt.rc("ytick", labelsize = 25)
    plt.rc("legend", fontsize = 30)
    fig_no = 0

    def __init__(self, plot_type, rows=1):
        plt.close("all")
        self.fig_no += 1
        if plot_type == 'single_plot':
            self.plt_size = (10, 7 + (rows-1)*10)
            plt.rc("lines", linewidth = 3.0)
            plt.rc("lines", markersize = 20.0)
        elif plot_type == 'two_columns':
            self.plt_size = (20, 10 + (rows-1)*10)
            plt.rc("lines", linewidth = 3.0)

        # Try to make a figure with the current backend
        try:
            plt.figure(self.fig_no, figsize = self.plt_size)
        except TclError:
            # Switch if a backend needs the display 
            plt.switch_backend('Agg')  
            plt.figure(self.fig_no, figsize = self.plt_size)



#{{{plot_look_nice
    def plot_look_nice(self, ax):
        """Makes a plot look nice"""
        
        # Set the legend
        # Ignore the warning that there were no legends
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            leg = ax.legend(loc="best", fancybox = True, numpoints=1)
        try:
            leg.get_frame().set_alpha(0.5)
        except:
            pass

        # Plot the grid
        ax.grid()
        # Make sure no collision between the ticks
        ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
        # Includes the xlabel if outside
        plt.tight_layout()
#}}}
#}}}



#{{{plot_look_nice_3d
def plot_look_nice_3d(ax, X, Y):
    """Makes a 3D plot look nice"""
    x = X[0,:]
    x_min = np.min(x)
    x_max = np.max(x)
    y = Y[:,0]
    y_min = np.min(y)
    y_max = np.max(y)
    ax.set_xlim(left=x_min, right=x_max)
    ax.set_ylim(bottom=y_min, top=y_max)        
#}}}



metric_prefixes =\
 {-24:'y'     ,\
  -21:'z'     ,\
  -18:'a'     ,\
  -15:'f'     ,\
  -12:'p'     ,\
  - 9:'n'     ,\
  - 6:'$\mu$' ,\
  - 3:'m'     ,\
  - 2:'c'     ,\
  - 1:'d'     ,\
    3:'k'     ,\
    6:'M'     ,\
    9:'G'     ,\
   12:'T'     ,\
   15:'P'     ,\
   18:'E'     ,\
   21:'Z'     ,\
   24:'Y'
 }

#{{{find_exponents
def find_exponents(max_number, use_metric_prefix=True,\
                   include_centi_deka= False):
    """Takes in a number, returns the exponent together with the metric
    prefix"""

    # Find the exponent of the number
    exponent_in = np.floor(np.log10(max_number))

    if use_metric_prefix:
        # Make a list of the keys in the dictionary
        prefixes = np.sort(metric_prefixes.keys())
        # We would like the numbers to range from 0 to 100
        exponent = -24
        for prefix in prefixes:
            if (include_centi_deka == False) and\
                (prefix ==-2 or prefix == -1):
                # Skip
                continue
            if exponent_in > prefix:
                exponent = prefix
            else:
                break
        metric_prefix = metric_prefixes[exponent]
    else:
        exponent = exponent_in
        metric_prefix = ''

    exponent = 10**(exponent)
    return exponent, metric_prefix
#}}}
