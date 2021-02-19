from __future__ import absolute_import
from __future__ import division
from builtins import range
from past.utils import old_div

# some function to plot nonlinear stuff
from .pb_corral import LinRes
from .ListDict import ListDictKey, ListDictFilt
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.artist as artist
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages

from reportlab.platypus import *
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.rl_config import defaultPageSize
from reportlab.lib.units import inch
from reportlab.graphics.charts.linecharts import HorizontalLineChart
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.lineplots import LinePlot
from reportlab.graphics.widgets.markers import makeMarker
from reportlab.lib import colors

from replab_x_vs_y import RL_Plot
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, MultipleLocator


class NLinResDraw(LinRes):
    def __init__(self, alldb):
        LinRes.__init__(self, alldb)

    def plotnlrhs(
        self,
        pp,
        field="Ni",
        yscale="linear",
        clip=0,
        xaxis="t",
        xscale="linear",
        xrange=1,
    ):

        colors = ["b", "g", "r", "c", "m", "y", "k", "b", "g", "r", "c", "m", "y", "k"]

        Modes = subset(self.db, "field", [field])  # pick  field
        comp = "ave"

        fig1 = plt.figure()
        adj = fig1.subplots_adjust(hspace=0.4, wspace=0.4)
        fig1.suptitle("Nonlinear contribution for " + field)
        props = dict(alpha=0.8, edgecolors="none")
        Nplots = self.nrun

        k = 0
        for j in list(set(Modes.path).union()):
            s = subset(Modes.db, "path", [j])  # pick a run folder - many modes
            dz = s.dz[0]
            data = s.ave[0]["nl"]
            x = np.array(list(range(data.size)))

            ax = fig1.add_subplot(round(old_div(Nplots, 2.0) + 1.0), 2, k + 1)
            ax.set_ylabel(r"$\frac{ddt_N}{ddt}$", fontsize=12, rotation="horizontal")
            k += 1
            ax.grid(True, linestyle="-", color=".75")
            try:
                ax.set_yscale(yscale, linthreshy=1e-13)
            except:
                ax.set_yscale("linear")
            i = 1
            ax.plot(x, data.flatten(), c=cm.jet(0.2 * i), linestyle="-")

            # data = np.array(ListDictKey(s.db,comp)) #pick component should be ok for a fixed dz key

            # we are not interested in looping over all modes

        fig1.savefig(pp, format="pdf")
        plt.close(fig1)

    # return 0


class subset(NLinResDraw):
    def __init__(self, alldb, key, valuelist, model=False):
        selection = ListDictFilt(alldb, key, valuelist)
        if len(selection) != 0:
            LinRes.__init__(self, selection)
            self.skey = key
            if model == True:
                self.model()
        else:
            LinRes.__init__(self, alldb)
            if model == True:
                self.model()
