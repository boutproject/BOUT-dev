from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from builtins import range
from .pb_draw import LinResDraw, subset
from .pb_corral import LinRes
from .pb_nonlinear import NLinResDraw
from pb_transport import Transport

import numpy as np
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.artist as artist
import matplotlib.ticker as ticker

# from matplotlib.ticker import FuncFormatter
# from matplotlib.ticker import ScalarFormatter

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

# for movie making
from multiprocessing import Queue, Pool
import multiprocessing
import subprocess

# uses  LinResDraw to make a pdf


class LinResPresent(LinResDraw, NLinResDraw, Transport):
    def __init__(self, alldb):
        LinResDraw.__init__(self, alldb)
        NLinResDraw.__init__(self, alldb)
        Transport.__init__(self, alldb)

    def show(
        self,
        filter=True,
        quick=False,
        pdfname="output2.pdf",
        debug=False,
        spectrum_movie=False,
    ):
        colors = ["b", "g", "r", "c", "m", "y", "k", "b", "g", "r", "c", "m", "y", "k"]
        pp = PdfPages("output.pdf")

        # start by removing modes above the maxN threshold
        modelist = []
        [
            modelist.append(list(self.modeid[p]))
            for p in range(self.nmodes)
            if self.mn[p][1] <= self.maxN[p]
        ]

        s = subset(self.db, "modeid", modelist)

        try:
            # fig = Figure(figsize=(6,6))
            # fig = plt.figure()
            dz0 = list(set(s.dz).union())[0]
            ss = subset(s.db, "dz", [dz0])

            # show initial condition and the first step after
            s.plotvsK(
                pp,
                yscale="log",
                xscale="log",
                t=[0, 1, -1],
                overplot=False,
                comp="amp",
                trans=True,
            )

            if spectrum_movie:
                ss.savemovie()

        except:
            print("no scatter")
        # 2D true NM spectrum with color code and boxes around spectral res regions log scale

        plt.figure()
        i = 0
        for j in list(
            set(s.dz).union()
        ):  # looping over runs, over unique 'dz' key values

            ss = subset(s.db, "dz", [j])  # subset where dz = j
            plt.scatter(ss.MN[:, 1], ss.MN[:, 0], c=colors[i])
            plt.annotate(str(j), (ss.MN[0, 1], ss.MN[0, 0]))
            i += 1

        plt.title(" Ni spectrum at t=0, all x")
        plt.ylabel("M -parallel")
        plt.xlabel("N -  axisymmteric")
        plt.xscale("log")
        plt.grid(True, linestyle="-", color=".75")

        try:
            plt.savefig(pp, format="pdf")
        except:
            print("FAILED TO save 1st part")

        plt.close()

        # for elem in self.meta['evolved']['v']:
        #    s.plotnl(pp

        if self.meta["nonlinear"]["v"] == "true":
            self.plotnlrhs(pp)

        if self.meta["transport"] == "true":
            self.plotnlrms(pp)

        for elem in self.meta["evolved"]:
            s.plotmodes(
                pp,
                yscale="symlog",
                comp="phase",
                linestyle=".",
                field=elem,
                summary=False,
            )
            s.plotmodes(pp, yscale="symlog", field=elem, summary=False)
            print(elem)
            try:
                s.plotmodes(
                    pp, yscale="symlog", field=elem, comp="gamma_i", summary=False
                )
            except:
                print("gamma_i plot for " + elem + " failed")

        # s.plotmodes(pp,yscale='symlog',summary=False)

        modelist = []
        # maxZ =
        # [modelist.append([1,p+1]) for p in range(maxZ-1)]
        [
            modelist.append(list(self.modeid[p]))
            for p in range(self.nmodes)
            if self.mn[p][1] <= self.maxN[p]
        ]
        ss = subset(s.db, "mn", modelist)

        if debug:  # just a few problematic slides
            fig1 = plt.figure()
            pp_bug = PdfPages("debug.pdf")
            # ss.plotmodes(pp_bug,yscale='symlog',comp='phase',summary=False)
            s.plotfreq2(pp_bug, xscale="log", yscale="symlog", overplot=True)
            ss.plotgamma(pp_bug, xscale="log", yscale="symlog", overplot=True)
            ss.plottheory(pp_bug)
            ss.plottheory(pp_bug, comp="freq")
            fig1.savefig(pp_bug, format="pdf")
            pp_bug.close()
            pp.close()
            return 0

        dir(ss)
        ss.plotmodes(pp, yscale="log", debug=True, summary=False)
        ss.plotmodes(pp, yscale="symlog", comp="phase", summary=False)
        ss.plotmodes(pp, yscale="symlog", comp="phase", field="rho", summary=False)
        print(dir(ss))

        # ss.plotmodes(pp,yscale='log',comp='phase',clip=True)

        # ss.plotfreq2(pp,xscale='log',yscale='linear',overplot=False)
        for elem in self.meta["evolved"]:
            ss.plotfreq2(
                pp, xscale="log", yscale="symlog", field=elem, overplot=True, trans=True
            )

        # ss.plotfreq2(pp,xscale='log',yscale='symlog',field='rho',overplot=True)

        if quick == True:
            pp.close()
            s.printmeta(pp)

            # plt.savefig(pp, format='pdf')
            return 0

        all_fields = list(set(s.field).union())

        s.plotgamma(pp, xscale="log", yscale="linear", overplot=True, trans=True)

        s.plotgamma(pp, yscale="symlog", xscale="log", overplot=True)
        s.plotgamma(pp, yscale="symlog", xscale="log", field="rho", overplot=True)

        try:
            s.plotfreq2(pp, xscale="log", yscale="linear", overplot=True)
            # s.plotfreq2(pp,xscale='log',yscale='symlog',overplot=False)
            s.plotfreq2(pp, xscale="log", yscale="symlog", field="rho", overplot=True)

            # s.plotfreq2(pp,xscale='log',yscale='linear')
        except:
            print("something terrible")

        s.plotradeigen(pp, yscale="linear")
        # s.plotradeigen(pp,field ='Vi',yscale='linear')
        s.plotradeigen(pp, field="rho", yscale="log")

        pp.close()
        s.printmeta(pp, filename=pdfname)  # append a metadata header
