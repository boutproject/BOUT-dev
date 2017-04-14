#!/usr/bin/env python

"""Init file for pre and post processing"""

import os
import matplotlib.pylab as plt

# Set proper backend for the display
try:
    os.environ["DISPLAY"]
except KeyError:
    plt.switch_backend("Agg")
