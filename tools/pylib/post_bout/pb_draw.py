from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div

# some standard analytic stuff to plot, if appending just overplot gam or omeg
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


class LinResDraw(LinRes):
    def __init__(self, alldb):
        LinRes.__init__(self, alldb)

    def plottheory(
        self, pp, m=1, canvas=None, comp="gamma", field="Ni", allroots=False
    ):
        if len(self.models) == 0:
            try:
                self.models = []
                self.models.append(_model(self))  # create a list to contain models
                self.models.append(
                    _model(self, haswak=True, name="haswak")
                )  # another model
            except:
                return 0

        s = subset(self.db, "field", [field])

        modelist = []

        [modelist.append([m, n + 1]) for n in range(min(s.maxN) - 1)]

        s = subset(s.db, "mn", modelist)

        allk = s.k[:, 1, old_div(s.nx, 2)]
        ki = np.argsort(allk)

        ownpage = False
        if canvas is None:
            ownpage = True

        if ownpage:  # if not an overplot
            fig1 = plt.figure()
            canvas = fig1.add_subplot(1, 1, 1)

        label = "gamma analytic"

        # if comp=='gamma':
        #     y = np.array(s.gammamax)[ki]
        # else:
        #     y = np.array(s.omegamax)[ki]

        for m in s.models:
            print(m.name)

        for i, m in enumerate(s.models):
            print(m.name, comp, m.soln[comp].shape)

            if allroots:
                for elem in m.soln[comp]:
                    # y = []
                    # elem has 2 or more elements
                    y = (np.array(elem)[ki]).flatten()  # n values
                    # y = y.astype('float')
                    print(y.shape)
                    canvas.plot(
                        (allk[ki]).flatten(), y, ",", label=label, c=cm.jet(0.2 * i)
                    )
            try:
                ymax = (np.array(m.soln[comp + "max"])[ki]).flatten()
                # ymax =  (np.array(m.gammamax)[ki]).flatten()
                ymax = ymax.astype("float")
                print(comp, " ymax:", ymax)
                canvas.plot(
                    (allk[ki]).flatten(), ymax, "-", label=label, c=cm.jet(0.2 * i)
                )
                # if comp=='gamma':
                #     y = (np.array(m.gammamax)[ki]).flatten()

                # else:
                #     y = (np.array(m.omegamax)[ki]).flatten()

            # print m.name, ':' ,y.astype('float')

            except:
                print("fail to add theory curve")

            canvas.annotate(m.name, (allk[ki[0]], 1.1 * ymax[0]), fontsize=8)
            canvas.annotate(m.name, (1.1 * allk[ki[-1]], 1.1 * ymax[-1]), fontsize=8)

            try:
                for i, m in enumerate(s.ref):
                    if not allroots:
                        y = (np.array(m.soln[comp])[ki]).flatten()
                        y = y.astype("float")
                        canvas.plot(
                            (allk[ki]).flatten(),
                            y,
                            "--",
                            label=label,
                            c=cm.jet(0.2 * i),
                        )
            except:
                print("no reference curve")

        if ownpage:  # set scales if this is its own plot
            # canvas.set_yscale('symlog',linthreshy=1e-13)
            # canvas.set_xscale('log')

            canvas.axis("tight")
            canvas.set_xscale("log")
            canvas.set_yscale("symlog")
            fig1.savefig(pp, format="pdf")
            plt.close(fig1)
        else:  # if not plot its probably plotted iwth sim data, print chi somewhere

            for i, m in enumerate(s.models):
                textstr = r"$\chi^2$" + "$=%.2f$" % (m.chi[comp].sum())
                print(textstr)
                # textstr = '$\L=%.2f$'%(m.chi[comp].sum())
                props = dict(boxstyle="square", facecolor="white", alpha=0.3)
                textbox = canvas.text(
                    0.1,
                    0.1,
                    textstr,
                    transform=canvas.transAxes,
                    fontsize=10,
                    verticalalignment="top",
                    bbox=props,
                )

    def plotomega(
        self,
        pp,
        canvas=None,
        field="Ni",
        yscale="linear",
        clip=0,
        xaxis="t",
        xscale="linear",
        xrange=1,
        comp="gamma",
        pltlegend="both",
        overplot=False,
        gridON=True,
        trans=False,
        infobox=True,
    ):

        colors = [
            "b.",
            "r.",
            "k.",
            "c.",
            "g.",
            "y.",
            "m.",
            "b.",
            "r.",
            "k.",
            "c.",
            "g.",
            "y.",
            "m.",
        ]
        colordash = [
            "b",
            "r",
            "k",
            "c",
            "g",
            "y",
            "m",
            "b",
            "r",
            "k",
            "c",
            "g",
            "y",
            "m",
        ]

        if canvas is None:
            ownpage = True
        else:
            ownpage = False

        if ownpage:
            fig1 = plt.figure()
            fig1.subplots_adjust(bottom=0.12)
            fig1.subplots_adjust(top=0.80)
            fig1.subplots_adjust(right=0.83)
            fig1.subplots_adjust(left=0.17)
            canvas = fig1.add_subplot(1, 1, 1)
            clonex = canvas.twinx()
            # if trans:
            #     cloney = canvas.twiny()

        dzhandles = []
        parhandles = []
        parlabels = []
        dzlabels = []

        m_shift = 1
        for q in np.array(list(range(1))) + m_shift:
            s = subset(self.db, "field", [field])  # pick  field
            maxZ = min(s.maxN)
            modelist = []
            [modelist.append([q, p + 1]) for p in range(maxZ - 1)]
            print(modelist)
            print(q, "in plotgamma")

        s = subset(s.db, "mn", modelist)

        xrange = old_div(s.nx, 2) - 2

        xrange = [old_div(s.nx, 2), old_div(s.nx, 2) + xrange]

        y = np.array(ListDictKey(s.db, comp))

        # y = s.gamma #nmodes x 2 x nx ndarray
        k = s.k  ##nmodes x 2 x nx ndarray k_zeta

        kfactor = np.mean(
            old_div(s.k_r[:, 1, old_div(s.nx, 2)], s.k[:, 1, old_div(s.nx, 2)])
        )  # good enough for now

        print(
            k[:, 1, old_div(s.nx, 2)].shape,
            y[:, 0, old_div(s.nx, 2)].shape,
            len(colors),
            ownpage,
        )  # ,k[:,1,s.nx/2],y[:,0,s.nx/2]

        parhandles.append(
            canvas.errorbar(
                np.squeeze(k[:, 1, old_div(s.nx, 2)]),
                np.squeeze(y[:, 0, old_div(s.nx, 2)]),
                yerr=np.squeeze(y[:, 1, old_div(s.nx, 2)]),
                fmt=colors[q],
            )
        )

        parlabels.append("m " + str(q))

        # loop over dz sets and connect with dotted line  . . .
        jj = 0

        ymin_data = np.max(np.array(ListDictKey(s.db, comp)))
        ymax_data = 0  # for bookeeping

        for p in list(set(s.path).union()):
            print(p, "in plotomega")

            sub_s = subset(s.db, "path", [p])
            j = sub_s.dz[0]
            # print sub_s.amp.shape
            s_i = np.argsort(sub_s.mn[:, 1])  # sort by 'local' m, global m is ok also
            # print s_i, sub_s.mn, sub_s.nx, jj
            y = np.array(ListDictKey(sub_s.db, comp))
            y_alt = 2.0 * np.array(ListDictKey(sub_s.db, comp))

            k = sub_s.k  ##
            if q == m_shift:  # fix the parallel mode
                dzhandles.append(
                    canvas.plot(
                        k[s_i, 1, old_div(sub_s.nx, 2)],
                        y[s_i, 0, old_div(sub_s.nx, 2)],
                        color=colordash[jj],
                        alpha=0.5,
                    )
                )
                # clonex.plot(k[s_i,1,sub_s.nx/2],
                # y_alt[s_i,0,sub_s.nx/2],color=colordash[jj],alpha=.5)
                # if np.any(sub_s.trans) and trans:
                #     comp_r = comp+'_r'
                #     k_r =  sub_s.k_r
                #     y2 = np.array(ListDictKey(sub_s.db,comp_r))
                #     cloney.plot(k[s_i,1,sub_s.nx/2],
                #                 y2[s_i,0,sub_s.nx/2],'k.',ms = 3)

                ymin_data = np.min([np.min(y[s_i, 0, old_div(sub_s.nx, 2)]), ymin_data])
                ymax_data = np.max([np.max(y[s_i, 0, old_div(sub_s.nx, 2)]), ymax_data])

                print("dzhandle color", jj)
                # dzlabels.append("DZ: "+ str(2*j)+r'$\pi$')
                dzlabels.append(j)

                if yscale == "log":
                    factor = 10
                else:
                    factor = 2
                print("annotating")
                canvas.annotate(
                    str(j),
                    (
                        k[s_i[0], 1, old_div(sub_s.nx, 2)],
                        y[s_i[0], 0, old_div(sub_s.nx, 2)],
                    ),
                    fontsize=8,
                )
                p = canvas.axvspan(
                    k[s_i[0], 1, old_div(sub_s.nx, 2)],
                    k[s_i[-1], 1, old_div(sub_s.nx, 2)],
                    facecolor=colordash[jj],
                    alpha=0.01,
                )
                print("done annotating")
            else:
                canvas.plot(
                    k[s_i, 1, old_div(sub_s.nx, 2)],
                    y[s_i, 0, old_div(sub_s.nx, 2)],
                    color=colordash[jj],
                    alpha=0.3,
                )

            jj = jj + 1

        dzhandles = np.array(dzhandles).flatten()
        dzlabels = np.array(dzlabels).flatten()

        dzlabels = list(set(dzlabels).union())

        dz_i = np.argsort(dzlabels)

        dzhandles = dzhandles[dz_i]
        dzlabels_cp = np.array(dzlabels)[dz_i]

        print(type(dzlabels), np.size(dzlabels))
        for i in range(np.size(dzlabels)):
            dzlabels[i] = "DZ: " + str(dzlabels_cp[i])  # +r"$\pi$"

        parlabels = np.array(parlabels).flatten()

        # if pltlegend =='both':    #

        print("legends")

        # l1 = legend(parhandles,parlabels,loc = 3,prop={'size':6})
        # l2 = legend(dzhandles,dzlabels,loc = 1,prop={'size':6})
        # plt.gca().add_artist(l1)

        # else:
        #    legend(dzhandles,dzlabels,loc=3,prop={'size':6})
        if overplot == True:
            try:
                self.plottheory(pp, canvas=canvas, comp=comp, field=field)
                # self.plottheory(pp,comp=comp)
            except:
                print("no theory plot")
        if infobox:
            textstr = "$\L_{\parallel}=%.2f$\n$\L_{\partial_r n}=%.2f$\n$B=%.2f$" % (
                s.meta["lpar"][old_div(s.nx, 2)],
                s.meta["L"][old_div(s.nx, 2), old_div(s.ny, 2)],
                s.meta["Bpxy"]["v"][old_div(s.nx, 2), old_div(s.ny, 2)],
            )
            props = dict(boxstyle="square", facecolor="white", alpha=0.3)
            textbox = canvas.text(
                0.82,
                0.95,
                textstr,
                transform=canvas.transAxes,
                fontsize=10,
                verticalalignment="top",
                bbox=props,
            )
            # leg = canvas.legend(handles,labels,ncol=2,loc='best',prop={'size':4},fancybox=True)
            # textbox.get_frame().set_alpha(0.3)
            # matplotlib.patches.Rectangle
            # p = patches.Rectangle((0, 0), 1, 1, fc="r")
            # p = str('L_par')
            # leg = canvas.legend([p], ["Red Rectangle"],loc='best',prop={'size':4})
            # leg.get_frame().set_alpha(0.3)

        # cloney.set_xlim(xmin,xmax)
        try:
            canvas.set_yscale(yscale)
            canvas.set_xscale(xscale)

            if yscale == "symlog":
                canvas.set_yscale(yscale, linthreshy=1e-13)
            if xscale == "symlog":
                canvas.set_xscale(xscale, linthreshy=1e-13)

            if gridON:
                canvas.grid()
        except:
            try:
                canvas.set_yscale("symlog")
            except:
                print("scaling failed completely")

        # print '[xmin, xmax, ymin, ymax]: ',[xmin, xmax, ymin, ymax]

        clonex.set_yscale(yscale)  # must be called before limits are set

        try:
            if yscale == "linear":
                formatter = ticker.ScalarFormatter()
                formatter.set_powerlimits((-2, 2))  # force scientific notation
                canvas.yaxis.set_major_formatter(formatter)
                clonex.yaxis.set_major_formatter(formatter)
                # canvas.useOffset=False
        except:
            print("fail 1")
        [xmin, xmax, ymin, ymax] = canvas.axis()

        if yscale == "symlog":
            clonex.set_yscale(yscale, linthreshy=1e-9)
        if xscale == "symlog":
            clonex.set_xscale(xscale, linthreshy=1e-9)
            # if np.any(s.trans) and trans:
        [xmin1, xmax1, ymin1, ymax1] = canvas.axis()
        if trans:
            try:
                cloney = canvas.twiny()
                # cloney.set_yscale(yscale)
                cloney.set_xscale(xscale)
                [xmin1, xmax1, ymin2, ymax2] = canvas.axis()

                if xscale == "symlog":
                    cloney.set_xscale(xscale, linthreshy=1e-9)
                if yscale == "symlog":
                    cloney.set_yscale(yscale, linthreshy=1e-9)
                if yscale == "linear":
                    cloney.yaxis.set_major_formatter(formatter)
            except:
                print("fail trans")
                #     cloney.useOffset=False

            # if xscale =='symlog' and trans:
            #     cloney.set_yscale(yscale,linthreshy=1e-9)
            #     cloney.set_xscale(xscale,linthreshy=1e-9)

        Ln_drive_scale = s.meta["w_Ln"][0] ** -1
        # Ln_drive_scale = 2.1e3
        clonex.set_ylim(Ln_drive_scale * ymin, Ln_drive_scale * ymax)

        try:
            if trans:
                # k_factor = #scales from k_zeta to k_perp
                cloney.set_xlim(kfactor * xmin, kfactor * xmax)
                # if np.any(sub_s.trans) and trans:
                #     comp_r = comp+'_r'
                #     y2 = np.array(ListDictKey(sub_s.db,comp_r))
                #    #canvas.plot(k[s_i,1,sub_s.nx/2],
                #    #           y2[s_i,0,sub_s.nx/2],'k.',ms = 3)

                #     cloney.plot(k_r[s_i,1,sub_s.nx/2],
                #                y2[s_i,0,sub_s.nx/2],'k.',ms = 3)
                #     print 'np.sum(np.abs(y-y2)): ',np.sum(np.abs(y-y2)),comp_r
                # kfactor =1.0
                # cloney.set_xlim(xmin,xmax)
                cloney.set_ylim(
                    ymin, ymax
                )  # because cloney shares the yaxis with canvas it may overide them, this fixes that
                cloney.set_xlabel(r"$k_{\perp} \rho_{ci}$", fontsize=18)
        except:
            print("moar fail")
            # clonex.set_xscale(xscale)

        # except:
        #     #canvas.set_xscale('symlog', linthreshx=0.1)
        #     print 'extra axis FAIL'

        # if yscale == 'linear':
        # canvas.yaxis.set_major_locator(ticker.LinearLocator(numticks=8))

        # minorLocator   = MultipleLocator(.005)
        # canvas.yaxis.set_minor_locator(minorLocator)
        # spawn another y label

        # clone = canvas.twinx()
        # s2 = np.sin(2*np.pi*t)
        # ax2.plot(x, s2, 'r.')

        # ion_acoust_str = r"$\frac{c_s}{L_{\partial_r n}}}$"

        if comp == "gamma":
            canvas.set_ylabel(
                r"$\frac{\gamma}{\omega_{ci}}$", fontsize=18, rotation="horizontal"
            )
            clonex.set_ylabel(
                r"$\frac{\gamma}{\frac{c_s}{L_n}}$",
                color="k",
                fontsize=18,
                rotation="horizontal",
            )
        if comp == "freq":
            canvas.set_ylabel(
                r"$\frac{\omega}{\omega_{ci}}$", fontsize=18, rotation="horizontal"
            )
            clonex.set_ylabel(
                r"$\frac{\omega}{\frac{c_s}{L_n}}$",
                color="k",
                fontsize=18,
                rotation="horizontal",
            )

        if comp == "amp":
            canvas.set_ylabel(r"$A_k$", fontsize=18, rotation="horizontal")
            clonex.set_ylabel(
                r"$\frac{A_k}{A_{max}}$", color="k", fontsize=18, rotation="horizontal"
            )

        canvas.set_xlabel(r"$k_{\zeta} \rho_{ci}$", fontsize=18)

        title = comp + " computed from " + field
        # canvas.set_title(title,fontsize=14)
        fig1.suptitle(title, fontsize=14)

        if ownpage:
            try:
                fig1.savefig(pp, format="pdf")
            except:
                print("pyplt doesnt like you")
            plt.close(fig1)

    def plotfreq(
        self, pp, field="Ni", clip=0, xaxis="t", xscale="linear", yscale="linear"
    ):
        # colors = ['b','g','r','c','m','y','k','b','g','r','c','m','y','k']
        colors = [
            "b.",
            "g.",
            "r.",
            "c.",
            "m.",
            "y.",
            "k.",
            "b.",
            "g.",
            "r.",
            "c.",
            "m.",
            "y",
            "k",
        ]
        plt.figure()

        # s = subset(self.db,'field',[field]) #pick  field

        for q in range(4):
            s = subset(self.db, "field", [field])  # pick  field across all dz sets
            modelist = []
            [modelist.append([q + 1, p + 1]) for p in range(5)]
            print(q, "in plotgamma")
            s = subset(s.db, "mn", modelist)

            gamma = s.freq  # nmodes x 2 x nx ndarray
            k = s.k  ##nmodes x 2 x nx ndarray

            plt.errorbar(
                k[:, 1, old_div(s.nx, 2)],
                gamma[:, 0, old_div(s.nx, 2)],
                yerr=gamma[:, 1, old_div(s.nx, 2)],
                fmt=colors[q],
            )
            plt.plot(
                k[:, 1, old_div(s.nx, 2)],
                gamma[:, 0, old_div(s.nx, 2)],
                "k:",
                alpha=0.3,
            )

            # loop over dz sets and connect with dotted line  . . .
            for j in list(set(s.dz).union()):
                # print j,len(s.mn)
                sub_s = subset(s.db, "dz", [j])
                gamma = sub_s.gamma
                k = sub_s.k  ##
                plt.plot(
                    k[:, 1, old_div(sub_s.nx, 2)],
                    gamma[:, 0, old_div(sub_s.nx, 2)],
                    "k:",
                    alpha=0.1,
                )

        try:
            plt.yscale(yscale)
        except:
            print("yscale fail")

        try:
            plt.xscale(yscale)
        except:
            plt.xscale("symlog")
        plt.xlabel(r"$k \rho_{ci}$", fontsize=14)
        plt.ylabel(r"$\frac{\omega}{\omega_{ci}}$", fontsize=14)
        # plt.title(r'$\frac{\omega}\{\omega_{ci}}$ '+ 'computed from'+field+ 'field',fontsize=10)

        plt.savefig(pp, format="pdf")
        plt.close()

    def plotgamma(
        self,
        pp,
        field="Ni",
        yscale="symlog",
        clip=0,
        xaxis="t",
        xscale="linear",
        xrange=1,
        comp="gamma",
        overplot=False,
        trans=True,
    ):
        self.plotomega(
            pp,
            field=field,
            yscale=yscale,
            clip=clip,
            xaxis=xaxis,
            xscale=xscale,
            xrange=xrange,
            comp=comp,
            overplot=overplot,
            trans=trans,
        )

    def plotfreq2(
        self,
        pp,
        field="Ni",
        yscale="symlog",
        clip=0,
        xaxis="t",
        xscale="linear",
        xrange=1,
        comp="freq",
        overplot=False,
        trans=True,
    ):
        self.plotomega(
            pp,
            field=field,
            yscale=yscale,
            clip=clip,
            xaxis=xaxis,
            xscale=xscale,
            xrange=xrange,
            comp=comp,
            overplot=overplot,
            trans=trans,
        )

    def plotvsK(
        self,
        pp,
        rootfig=None,
        field="Ni",
        yscale="log",
        clip=0,
        xaxis="t",
        xscale="linear",
        xrange=1,
        comp="amp",
        pltlegend="both",
        overplot=False,
        gridON=True,
        trans=False,
        infobox=True,
        m=1,
        t=[0],
        file=None,
        save=True,
    ):
        colors = [
            "b.",
            "r.",
            "k.",
            "c.",
            "g.",
            "y.",
            "m.",
            "b.",
            "r.",
            "k.",
            "c.",
            "g.",
            "y.",
            "m.",
        ]
        colordash = [
            "b",
            "r",
            "k",
            "c",
            "g",
            "y",
            "m",
            "b",
            "r",
            "k",
            "c",
            "g",
            "y",
            "m",
        ]

        if rootfig is None:
            ownpage = True
        else:
            ownpage = False

        if ownpage:
            fig1 = plt.figure()
            fig1.subplots_adjust(bottom=0.12)
            fig1.subplots_adjust(top=0.80)
            fig1.subplots_adjust(right=0.83)
            fig1.subplots_adjust(left=0.17)
            canvas = fig1.add_subplot(1, 1, 1)
            clonex = canvas.twinx()
            # if trans:
            #     cloney = canvas.twiny()
        else:
            canvas = rootfig.add_subplot(1, 1, 1)

        dzhandles = []
        parhandles = []
        parlabels = []
        dzlabels = []

        # pick the modes
        m_shift = m
        for q in np.array(list(range(1))) + m_shift:
            s = subset(self.db, "field", [field])  # pick  field
            maxZ = min(s.maxN)
            modelist = []
            [modelist.append([q, p + 1]) for p in range(maxZ - 1)]
        # print q,'in plotgamma'
        s = subset(s.db, "mn", modelist)

        # set x-range
        xrange = old_div(s.nx, 2) - 2
        xrange = [old_div(s.nx, 2), old_div(s.nx, 2) + xrange]

        # pull up the data
        y = np.array(ListDictKey(s.db, comp))
        print("y.shape", y.shape)

        # in case multiple timesteps are indicated
        all_y = []
        all_yerr = []
        if comp == "amp":
            for elem in t:
                all_y.append(np.squeeze(y[:, elem, :]))
                all_yerr.append(np.squeeze(0 * y[:, elem, :]))
            ynorm = np.max(all_y)
            # all_y = np.array(np.squeeze(all_y))
            # all_yerr = np.array(np.squeeze(all_yerr))

        else:
            all_y.append(np.squeeze(y[:, 0, :]))
            all_yerr.append(np.squeeze(y[:, 1, :]))
            ynorm = s.meta["w_Ln"][0]

        k = s.k  ##nmodes x 2 x nx ndarray k_zeta

        kfactor = np.mean(
            old_div(s.k_r[:, 1, old_div(s.nx, 2)], s.k[:, 1, old_div(s.nx, 2)])
        )  # good enough for now

        for elem in range(np.size(t)):
            # print 'printing line' , elem
            errorline = parhandles.append(
                canvas.errorbar(
                    k[:, 1, old_div(s.nx, 2)],
                    all_y[elem][:, old_div(s.nx, 2)],
                    yerr=all_yerr[elem][:, old_div(s.nx, 2)],
                    fmt=colors[q],
                )
            )

        parlabels.append("m " + str(q))

        # loop over dz sets and connect with dotted line  . . .
        jj = 0  # will reference dz color

        ymin_data = np.max(np.array(ListDictKey(s.db, comp)))
        ymax_data = 0  # for bookeeping

        for p in list(set(s.path).union()):
            sub_s = subset(s.db, "path", [p])
            j = sub_s.dz[0]
            # print sub_s.amp.shape
            s_i = np.argsort(sub_s.mn[:, 1])  # sort by 'local' m, global m is ok also
            # print s_i, sub_s.mn, sub_s.nx, jj
            y = np.array(ListDictKey(sub_s.db, comp))
            y_alt = 2.0 * np.array(ListDictKey(sub_s.db, comp))
            all_y = []
            all_yerr = []
            if comp == "amp":
                for elem in t:
                    all_y.append(np.squeeze(y[:, elem, :]))
                    all_yerr.append(np.squeeze(0 * y[:, elem, :]))
            else:
                all_y.append(np.squeeze(y[:, 0, :]))
                all_yerr.append(np.squeeze(y[:, 1, :]))

            k = sub_s.k  ##

            for elem in range(np.size(t)):
                if q == m_shift:  # fix the parallel mode
                    dzhandles.append(
                        canvas.plot(
                            k[s_i, 1, old_div(sub_s.nx, 2)],
                            all_y[elem][s_i, old_div(sub_s.nx, 2)],
                            color=colordash[jj],
                            alpha=0.5,
                        )
                    )

                    ymin_data = np.min(
                        [np.min(y[s_i, old_div(sub_s.nx, 2)]), ymin_data]
                    )
                    ymax_data = np.max(
                        [np.max(y[s_i, old_div(sub_s.nx, 2)]), ymax_data]
                    )

                    dzlabels.append(j)

                    if yscale == "log":
                        factor = 10
                    else:
                        factor = 2
                    # print 'annotating'
                    canvas.annotate(
                        str(j),
                        (
                            k[s_i[0], 1, old_div(sub_s.nx, 2)],
                            y[elem][s_i[0], old_div(sub_s.nx, 2)],
                        ),
                        fontsize=8,
                    )
                    # p = canvas.axvspan(k[s_i[0],1,sub_s.nx/2], k[s_i[-1],1,sub_s.nx/2],
                    #               facecolor=colordash[jj], alpha=0.01)
                    print("done annotating")
                else:
                    canvas.plot(
                        k[s_i, 1, old_div(sub_s.nx, 2)],
                        y[elem][s_i, old_div(sub_s.nx, 2)],
                        color=colordash[jj],
                        alpha=0.3,
                    )

                jj = jj + 1

        dzhandles = np.array(dzhandles).flatten()
        dzlabels = np.array(dzlabels).flatten()

        dzlabels = list(set(dzlabels).union())

        dz_i = np.argsort(dzlabels)

        dzhandles = dzhandles[dz_i]
        dzlabels_cp = np.array(dzlabels)[dz_i]

        # print type(dzlabels), np.size(dzlabels)
        for i in range(np.size(dzlabels)):
            dzlabels[i] = "DZ: " + str(dzlabels_cp[i])  # +r"$\pi$"

        parlabels = np.array(parlabels).flatten()

        if overplot == True:
            try:
                self.plottheory(pp, canvas=canvas, comp=comp, field=field)
                # self.plottheory(pp,comp=comp)
            except:
                print("no theory plot")
        if infobox:
            textstr = "$\L_{\parallel}=%.2f$\n$\L_{\partial_r n}=%.2f$\n$B=%.2f$" % (
                s.meta["lpar"][old_div(s.nx, 2)],
                s.meta["L"][old_div(s.nx, 2), old_div(s.ny, 2)],
                s.meta["Bpxy"]["v"][old_div(s.nx, 2), old_div(s.ny, 2)],
            )
            props = dict(boxstyle="square", facecolor="white", alpha=0.3)
            textbox = canvas.text(
                0.82,
                0.95,
                textstr,
                transform=canvas.transAxes,
                fontsize=10,
                verticalalignment="top",
                bbox=props,
            )
            # leg = canvas.legend(handles,labels,ncol=2,loc='best',prop={'size':4},fancybox=True)

        # cloney.set_xlim(xmin,xmax)
        try:
            canvas.set_yscale(yscale)
            canvas.set_xscale(xscale)

            if yscale == "symlog":
                canvas.set_yscale(yscale, linthreshy=1e-13)
            if xscale == "symlog":
                canvas.set_xscale(xscale, linthreshy=1e-13)

            if gridON:
                canvas.grid()
        except:
            try:
                canvas.set_yscale("symlog")
            except:
                print("scaling failed completely")

        ##################################################################

        if ownpage and rootfig is None:
            clonex.set_yscale(yscale)  # must be called before limits are set

            try:
                if yscale == "linear":
                    formatter = ticker.ScalarFormatter()
                    formatter.set_powerlimits((-2, 2))  # force scientific notation
                    canvas.yaxis.set_major_formatter(formatter)
                    clonex.yaxis.set_major_formatter(formatter)
                # canvas.useOffset=False
            except:
                print("fail 1")
            [xmin, xmax, ymin, ymax] = canvas.axis()

            if yscale == "symlog":
                clonex.set_yscale(yscale, linthreshy=1e-9)
            if xscale == "symlog":
                clonex.set_xscale(xscale, linthreshy=1e-9)
            # if np.any(s.trans) and trans:
            [xmin1, xmax1, ymin1, ymax1] = canvas.axis()
            if trans:
                try:
                    cloney = canvas.twiny()
                    # cloney.set_yscale(yscale)
                    cloney.set_xscale(xscale)
                    [xmin1, xmax1, ymin2, ymax2] = canvas.axis()

                    if xscale == "symlog":
                        cloney.set_xscale(xscale, linthreshy=1e-9)
                    if yscale == "symlog":
                        cloney.set_yscale(yscale, linthreshy=1e-9)
                    if yscale == "linear":
                        cloney.yaxis.set_major_formatter(formatter)
                except:
                    print("fail trans")

            Ln_drive_scale = s.meta["w_Ln"][0] ** -1
            # Ln_drive_scale = 2.1e3
            # clonex.set_ylim(Ln_drive_scale*ymin, Ln_drive_scale*ymax)
            clonex.set_ylim(ynorm ** -1 * ymin, ynorm ** -1 * ymax)

            try:
                if trans:
                    # k_factor = #scales from k_zeta to k_perp
                    cloney.set_xlim(kfactor * xmin, kfactor * xmax)

                    cloney.set_ylim(
                        ymin, ymax
                    )  # because cloney shares the yaxis with canvas it may overide them, this fixes that
                    cloney.set_xlabel(r"$k_{\perp} \rho_{ci}$", fontsize=18)
            except:
                print("moar fail")
            # clonex.set_xscale(xscale)

            # ion_acoust_str = r"$\frac{c_s}{L_{\partial_r n}}}$"

            if comp == "gamma":
                canvas.set_ylabel(
                    r"$\frac{\gamma}{\omega_{ci}}$", fontsize=18, rotation="horizontal"
                )
                clonex.set_ylabel(
                    r"$\frac{\gamma}{\frac{c_s}{L_n}}$",
                    color="k",
                    fontsize=18,
                    rotation="horizontal",
                )
            if comp == "freq":
                canvas.set_ylabel(
                    r"$\frac{\omega}{\omega_{ci}}$", fontsize=18, rotation="horizontal"
                )
                clonex.set_ylabel(
                    r"$\frac{\omega}{\frac{c_s}{L_n}}$",
                    color="k",
                    fontsize=18,
                    rotation="horizontal",
                )

            if comp == "amp":
                canvas.set_ylabel(r"$A_k$", fontsize=18, rotation="horizontal")
                clonex.set_ylabel(
                    r"$\frac{A_k}{A_{max}}$",
                    color="k",
                    fontsize=18,
                    rotation="horizontal",
                )

            canvas.set_xlabel(r"$k_{\zeta} \rho_{ci}$", fontsize=18)

            title = comp + " computed from " + field
            # canvas.set_title(title,fontsize=14)
            fig1.suptitle(title, fontsize=14)

        if not ownpage:
            print("probably for a movie")
            fig1 = rootfig
            # canvasjunk = fig1.add_subplot(1,1,1)
            # canvasjunk = canvas

        if save:
            if file is None:
                try:
                    fig1.savefig(pp, format="pdf")
                except:
                    print("pyplt doesnt like you")
            else:
                try:
                    fig1.savefig(file, dpi=200)
                except:
                    print("no movie for you ;(")

        if ownpage:
            # fig1.close()
            plt.close(fig1)

    def plotmodes(
        self,
        pp,
        field="Ni",
        comp="amp",
        math="1",
        ylim=1,
        yscale="symlog",
        clip=False,
        xaxis="t",
        xscale="linear",
        xrange=1,
        debug=False,
        yaxis=r"$\frac{Ni}{Ni_0}$",
        linestyle="-",
        summary=True,
    ):

        Nplots = self.nrun

        colors = ["b", "g", "r", "c", "m", "y", "k", "b", "g", "r", "c", "m", "y", "k"]
        styles = ["^", "s"]

        fig1 = plt.figure()

        fig2 = plt.figure()

        Modes = subset(self.db, "field", [field])  # pick  field

        adj = fig2.subplots_adjust(hspace=0.4, wspace=0.4)
        fig2.suptitle("Dominant mode " + comp + " for  " + field)
        props = dict(alpha=0.8, edgecolors="none")

        allcurves = fig1.add_subplot(1, 1, 1)
        fig1.suptitle("Dominant mode behavior for  " + field)

        modenames = []
        k = 0

        for j in list(set(Modes.path).union()):  #
            s = subset(Modes.db, "path", [j])  # pick run
            dz = s.dz[0]
            xr = list(
                range(
                    old_div(s.nx, 2) - old_div(xrange, 2),
                    old_div(s.nx, 2) + old_div(xrange, 2) + 1,
                )
            )
            data = np.array(
                ListDictKey(s.db, comp)
            )  # pick component should be ok for a fixed dz key

            data = data  # + 1e-32 #hacky way to deal with buggy scaling
            ax = fig2.add_subplot(round(old_div(Nplots, 3.0) + 1.0), 3, k + 1)

            ax.grid(True, linestyle="-", color=".75")
            handles = []
            # modenames.append(str(j))

            # find the "biggest" mode for this dz
            d = data[:, s.nt[0] - 1, :]  # nmode X nx array
            # d = s.gamma[:,2,:]
            where = d == np.nanmax(d)
            z = where.nonzero()  # mode index and n index
            imax = z[0][0]
            # xi_max = z[1][0]
            xi_max = old_div(s.nx, 2)

            if debug and yscale == "log":
                gamma = np.array(ListDictKey(s.db, "gamma"))  # nmodes x 2 x nx

            for i in range(s.nmodes):
                if math == "gamma":
                    out = old_div(np.gradient(data[i, :, xr])[1], data[i, :, xr])
                else:
                    out = data[i, 2:, xi_max]  # skip the first 2 points

                if xaxis == "t":
                    # print 'out.size', out.size, out.shape
                    x = np.array(list(range(out.size)))
                    # plt.plot(x,out.flatten(),c=colors[k])
                    label = str(s.mn[i])
                    # handles.append(ax.plot(x,out.flatten(),
                    #                 c=cm.jet(1.*k),label = label))
                    ax.plot(
                        x, out.flatten(), c=cm.jet(0.2 * i), label=label, linestyle="-"
                    )

                else:
                    x = np.array(ListDictKey(s.db, xaxis))[i, :, xr]
                    # x #an N? by nx array
                    print(x[:, 1], out[:, 0])
                    plt.scatter(x[:, 1], out[:, 0])  # ,c=colors[k])
                    ax.scatter(
                        x[:, 1], out[:, 0]
                    )  # ,c=colors[k])#,alpha = (1 +i)/s.nmodes)

                    # detect error (bar data
                    print("error bars:", x, out)

            # ax.legend(handles,labels,loc='best',prop={'size':6})

            formatter = ticker.ScalarFormatter()
            formatter.set_powerlimits((0, 0))
            ax.xaxis.set_major_formatter(formatter)
            # ax.axis('tight')
            if yscale == "linear":
                ax.yaxis.set_major_formatter(formatter)
            if yscale == "symlog":
                ax.set_yscale("symlog", linthreshy=1e-13)
            else:
                try:
                    ax.set_yscale(yscale)
                except:
                    print("may get weird axis")
                    ax.set_yscale("symlog")
            # if comp=='phase' or yscale=='linear':
            #         ax.set_xscale('symlog',linthreshx=1.0)

            ax.set_xscale(xscale)

            ax.axis("tight")
            artist.setp(ax.axes.get_xticklabels(), fontsize=6)
            artist.setp(ax.axes.get_yticklabels(), fontsize=8)
            # artist.setp(ax.axes.get_yscale(), fontsize=8)
            ax.set_title(str(dz), fontsize=10)
            ax.set_xlabel(xaxis)
            handles, labels = ax.get_legend_handles_labels()
            leg = ax.legend(
                handles, labels, ncol=2, loc="best", prop={"size": 4}, fancybox=True
            )
            leg.get_frame().set_alpha(0.3)
            # x = s.Rxy[imax,:,s.ny/2]

            t0 = 2
            if clip == True:
                t0 = round(old_div(s.nt[0], 3))
            y = np.squeeze(data[imax, t0:, xi_max])
            x = np.array(list(range(y.size)))

            print(imax, xi_max)

            label = (
                str([round(elem, 3) for elem in s.MN[imax]])
                + str(s.mn[imax])
                + " at x= "
                + str(xi_max)
                + " ,"
                + str(
                    round(
                        old_div(s.gamma[imax, 2, xi_max], s.gamma[imax, 0, xi_max]), 3
                    )
                )
                + "%  "
                + str(round(s.gamma[imax, 0, xi_max], 4))
            )

            short_label = str(dz)
            print(short_label, x.shape, y.shape)
            allcurves.plot(x, y, ".", c=cm.jet(1.0 * k / len(x)), label=label)
            # print len(x), k*len(x)/(Nplots+2),s.nrun
            allcurves.annotate(
                short_label,
                (x[k * len(x) / (Nplots + 1)], y[k * len(x) / (Nplots + 1)]),
                fontsize=8,
            )

            # modenames.append(str([round(elem,3) for elem in s.MN[imax]])
            #        +str(s.mn[imax])+' at x= '+str(xi_max)+' ,'+str(s.gamma[imax,2,xi_max]))

            if debug and yscale == "log":
                gam = gamma[imax, 0, xi_max]
                f0 = gamma[imax, 1, xi_max]
                allcurves.plot(x, f0 * np.exp(gam * s.dt[imax] * x), "k:")

            k += 1

        # if ylim:
        #    allcurves.set_ylim(data[,xi_max].min(),5*data[:,xi_max].max())

        fig2.savefig(pp, format="pdf")

        handles, labels = allcurves.get_legend_handles_labels()
        allcurves.legend(handles, labels, loc="best", prop={"size": 6})
        # allcurves.legend(modenames,loc='best',prop={'size':6})
        allcurves.set_title(
            field + " " + comp + ", all runs, " + yscale + " yscale", fontsize=10
        )
        allcurves.set_ylabel(yaxis)
        allcurves.set_xlabel(xaxis)

        if yscale == "linear":
            allcurves.yaxis.set_major_formatter(formatter)
        else:
            try:
                allcurves.set_yscale(yscale)
            except:
                print("may get weird axis scaling")
            if yscale == "log":
                allcurves.axis("tight")
            # allcurves.set_ylim(data.min(),data.max())
        # allcurves.set_yscale(yscale,nonposy='mask')

        # plt.xscale(xscale)

        # plt.legend(modenames,loc='best')
        if summary:
            fig1.savefig(pp, format="pdf")
            plt.close(fig1)

        plt.close(fig2)

    # except:
    #   print "Sorry you fail"

    def plotradeigen(
        self, pp, field="Ni", comp="amp", yscale="linear", xscale="linear"
    ):

        Nplots = self.nrun
        colors = ["b", "g", "r", "c", "m", "y", "k", "b", "g", "r", "c", "m", "y", "k"]
        fig1 = plt.figure()

        fig2 = plt.figure()
        adj = fig2.subplots_adjust(hspace=0.4, wspace=0.4)

        # canvas = FigureCanvas(fig)

        Modes = subset(self.db, "field", [field])

        k = 0
        fig2.suptitle("Dominant mode behavior for  " + field)
        props = dict(alpha=0.8, edgecolors="none")

        allcurves = fig1.add_subplot(1, 1, 1)
        fig1.suptitle("Dominant mode behavior for  " + field)

        modeleg = []

        for p in list(set(Modes.path).union()):
            print(p)
            s = subset(Modes.db, "path", [p])  # pick run
            # data = np.array(ListDictKey(s.db,comp)) #pick component
            j = s.dz[0]
            ax = fig2.add_subplot(round(old_div(Nplots, 3.0) + 1.0), 3, k + 1)
            ax.grid(True, linestyle="-", color=".75")
            data = np.array(ListDictKey(s.db, comp))  # pick component
            handles = []

            # find the "biggest" mode for this dz
            d = data[:, s.nt[0] - 1, :]  # nmode X nx array
            where = d == d.max()
            z = where.nonzero()  # mode index and n index
            imax = z[0][0]
            modeleg.append(
                str([round(elem, 3) for elem in s.MN[imax]]) + str(s.mn[imax])
            )

            # str(s.k[:,1,:][z])

            for i in range(s.nmodes):
                # out = mode[mode.ny/2,:]
                # print i,s.Rxynorm.shape,s.ny
                x = np.squeeze(s.Rxynorm[i, :, old_div(s.ny, 2)])
                y = data[i, s.nt[0] - 1, :]

                handles.append(ax.plot(x, y, c=cm.jet(1.0 * k / len(x))))

            formatter = ticker.ScalarFormatter()
            formatter.set_powerlimits((0, 0))
            ax.xaxis.set_major_formatter(formatter)

            if yscale == "linear":
                ax.yaxis.set_major_formatter(formatter)
            else:
                ax.set_yscale(yscale)

            artist.setp(ax.axes.get_xticklabels(), fontsize=6)
            artist.setp(ax.axes.get_yticklabels(), fontsize=8)
            # artist.setp(ax.axes.get_yscale(), fontsize=8)
            ax.set_title(str(j), fontsize=10)

            x = np.squeeze(s.Rxynorm[imax, :, old_div(s.ny, 2)])
            y = data[imax, s.nt[0] - 1, :]
            # allcurves.plot(x,y,c= colors[k])
            allcurves.plot(x, y, c=cm.jet(0.1 * k / len(x)))
            print(k)
            k = k + 1

        fig2.savefig(pp, format="pdf")
        if yscale == "linear":
            allcurves.yaxis.set_major_formatter(formatter)
        else:
            allcurves.set_yscale(yscale)
        try:
            allcurves.set_xscale(xscale)
        except:
            allcurves.set_xscale("symlog")

        # allcurves.xaxis.set_major_formatter(ticker.NullFormatter())
        allcurves.legend(modeleg, loc="best", prop={"size": 6})
        allcurves.set_xlabel(r"$\frac{x}{\rho_{ci}}$")
        allcurves.set_ylabel(r"$\frac{Ni}{Ni_0}$")
        fig1.savefig(pp, format="pdf")
        plt.close(fig1)
        plt.close(fig2)

    def plotmodes2(
        self,
        pp,
        field="Ni",
        comp="amp",
        math="1",
        ylim=1,
        yscale="symlog",
        clip=0,
        xaxis="t",
        xscale="linear",
        xrange=1,
        debug=False,
    ):

        Nplots = self.nrun
        Modes = subset(self.db, "field", [field])  # pick  field
        colors = ["b", "g", "r", "c", "m", "y", "k", "b", "g", "r", "c", "m", "y", "k"]

        fig = Figure()
        plt.figure()

        canvas = FigureCanvas(fig)
        k = 0
        nrow = round(old_div(Nplots, 3.0) + 1.0)
        ncol = 3
        # nrow = round(Nplots/3.0 + 1.0)
        # ncol = round(Nplots/3.0 + 1.0)
        f, axarr = plt.subplots(int(nrow), int(ncol))

        for p in list(set(Modes.path).union()):  #
            s = subset(Modes.db, "path", [p])  # pick run
            j = s.dz[0]
            xr = list(
                range(
                    old_div(s.nx, 2) - old_div(xrange, 2),
                    old_div(s.nx, 2) + old_div(xrange, 2) + 1,
                )
            )
            data = np.array(ListDictKey(s.db, comp))  # pick component
            # ax =fig.add_subplot(round(Nplots/3.0 + 1.0),3,k+1)

            for i in range(s.nmodes):
                out = data[i, :, xr]

            print(j, i)
            if xaxis == "t":
                x = list(range(out.size))
                # plt.scatter(x,out.flatten(),c=colors[k])
                plt.scatter(x, out.flatten(), c=cm.jet(1.0 * k / len(x)))
                # axarr[j%(ncol),j/ncol].scatter(x,out.flatten(),c=colors[k])#,alpha = (1 +i)/s.nmodes)
                axarr[old_div(j, ncol), j % (ncol)].scatter(
                    x, out.flatten(), c=cm.jet(1.0 * k / len(x))
                )  #

            else:
                x = np.array(ListDictKey(s.db, xaxis))[i, :, xr]
                # x #an N? by nx array
                print(x[:, 1], out[:, 0])
                plt.scatter(x[:, 1], out[:, 0])  # ,c=colors[k])
                axarr[j % (col), old_div(j, col)].scatter(
                    x[:, 1], out[:, 0]
                )  # ,c=colors[k])#,alpha = (1 +i)/s.nmodes)

                # detect error (bar data
                print("error bars:", x, out)

            axarr[old_div(j, ncol), j % (ncol)].set_yscale(yscale)
            axarr[old_div(j, ncol), j % (ncol)].set_xscale(xscale)
            axarr[old_div(j, ncol), j % (ncol)].set_title(str(j), fontsize=10)
            axarr[old_div(j, ncol), j % (ncol)].set_xlabel(xaxis)

            plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
            plt.setp([a.get_yticklabels() for a in axarr[:, ncol - 1]], visible=False)

            if ylim:
                axarr[j % (ncol), old_div(j, ncol)].set_ylim(data.min(), 5 * data.max())
            k += 1

        plt.title(field + " " + comp + ", all runs, " + yscale + " yscale", fontsize=10)
        plt.xlabel(xaxis)
        if ylim:
            plt.ylim(data.min(), 10 * data.max())

        plt.yscale(yscale, nonposy="mask")
        plt.xscale(xscale)

        fig.savefig(pp, format="pdf")
        plt.savefig(pp, format="pdf")

        plt.close()

        return 0

    def plotMacroDep(
        self,
        pp,
        field="Ni",
        yscale="symlog",
        clip=0,
        xaxis="t",
        xscale="linear",
        xrange=1,
    ):
        colors = ["b", "g", "r", "c", "m", "y", "k", "b", "g", "r", "c", "m", "y", "k"]
        plt.figure()

    def savemovie(
        self, field="Ni", yscale="log", xscale="log", moviename="spectrum.avi"
    ):

        print("Making movie animation.mpg - this make take a while")
        files = []

        for t in range(self.nt[0] - 3):
            print(t)
            filename = str("%03d" % (t + 1) + ".png")
            self.plotvsK(
                "dont need pp",
                yscale="log",
                t=[1, t + 2],
                xscale="log",
                overplot=False,
                comp="amp",
                trans=True,
                file=filename,
            )
            files.append(filename)

        command = (
            "mencoder",
            "mf://*.png",
            "-mf",
            "type=png:w=800:h=600:fps=10",
            "-ovc",
            "lavc",
            "-lavcopts",
            "vcodec=mpeg4",
            "-oac",
            "copy",
            "-o",
            moviename,
        )

        import subprocess, os

        subprocess.check_call(command)
        os.system("rm *png")

    def printmeta(self, pp, filename="output2.pdf", debug=False):

        import os
        from pyPdf import PdfFileWriter, PdfFileReader

        PAGE_HEIGHT = defaultPageSize[1]
        styles = getSampleStyleSheet()
        Title = "BOUT++ Results"
        Author = "Dmitry Meyerson"
        URL = ""
        email = "dmitry.meyerson@gmail.com"
        Abstract = """This document highlights some results from BOUT++ simulation"""
        Elements = []
        HeaderStyle = styles["Heading1"]
        ParaStyle = styles["Normal"]
        PreStyle = styles["Code"]

        def header(
            txt, style=HeaderStyle, klass=Paragraph, sep=0.3
        ):  # return styled text with a space
            s = Spacer(0.2 * inch, sep * inch)
            para = klass(txt, style)
            sect = [s, para]
            result = KeepTogether(sect)
            return result

        def p(txt):  # wrapper for header
            return header(txt, style=ParaStyle, sep=0.1)

        def pre(txt):  # return styled text with a space
            s = Spacer(0.1 * inch, 0.1 * inch)
            p = Preformatted(txt, PreStyle)
            precomps = [s, p]
            result = KeepTogether(precomps)
            return result

        def graphout(name, datain, xaxis=None):
            if xaxis is None:
                xaxis = list(range(datain.size))
            if xlabel is None:
                xlabel = ""
            if ylabel is None:
                ylabel = ""

            drawing = Drawing(400, 200)
            # data = [
            #    ((1,1), (2,2), (2.5,1), (3,3), (4,5)),
            #    ((1,2), (2,3), (2.5,2), (3.5,5), (4,6))
            #    ]
            dataview = [tuple([(xaxis[i], datain[i]) for i in range(datain.size)])]
            lp = LinePlot()
            lp.x = 50
            lp.y = 50
            lp.height = 125
            lp.width = 300
            lp.data = dataview
            lp.xValueAxis.xLabelFormat = "{mmm} {yy}"
            lp.lineLabels.fontSize = 6
            lp.lineLabels.boxStrokeWidth = 0.5
            lp.lineLabels.visible = 1
            lp.lineLabels.boxAnchor = "c"
            # lp.joinedLines = 1
            # lp.lines[0].symbol = makeMarker('FilledCircle')
            # lp.lines[1].symbol = makeMarker('Circle')
            # lp.lineLabelFormat = '%2.0f'
            # lp.strokeColor = colors.black
            # lp.xValueAxis.valueMin = min(xaxis)
            # lp.xValueAxis.valueMax = max(xaxis)
            # lp.xValueAxis.valueSteps = xaxis
            # lp.xValueAxis.labelTextFormat = '%2.1f'
            # lp.yValueAxis.valueMin = min(datain)
            # lp.yValueAxis.valueMax = max(datain)
            # lp.yValueAxis.valueSteps = [1, 2, 3, 5, 6]
            drawing.add(lp)
            return drawing

        def go():
            doc = SimpleDocTemplate("meta.pdf")
            doc.build(Elements)

        mytitle = header(Title)
        myname = header(Author, sep=0.1, style=ParaStyle)
        mysite = header(URL, sep=0.1, style=ParaStyle)
        mymail = header(email, sep=0.1, style=ParaStyle)
        abstract_title = header("ABSTRACT")
        myabstract = p(Abstract)
        head_info = [mytitle, myname, mysite, mymail, abstract_title, myabstract]
        Elements.extend(head_info)

        meta_title = header("metadata", sep=0)
        metasection = []
        metasection.append(meta_title)

        for i, elem in enumerate(self.meta):
            # if type(self.meta[elem]) != type(np.array([])):

            print(elem, type(self.meta[elem]))

            if type(self.meta[elem]) == type({}):
                print("{}")
                data = np.array(self.meta[elem]["v"])
                unit_label = str(self.meta[elem]["u"])
            else:
                data = np.array(self.meta[elem])
                unit_label = ""

            xaxis = np.squeeze(self.meta["Rxy"]["v"][:, old_div(self.ny, 2)])

            if data.shape == (self.nx, self.ny):
                datastr = np.squeeze(data[:, old_div(self.ny, 2)])
                # metasection.append(graphout('stuff',datastr,xaxis=xaxis))
                # metasection.append(RL_Plot(datastr,xaxis))

                metasection.append(RL_Plot(datastr, xaxis, linelabel=str(elem)))
            # metasection.append(RL_Plot(datastr,xaxis,xlabel='xlabel'))
            elif data.shape == self.nx:
                datastr = data
            # metasection.append(graphout('stuff',datastr,xaxis=xaxis))
            # metasection.append(RL_Plot(datastr,xaxis,linelabel=str(elem)))
            elif data.shape == (1,):
                data = data[0]
                metasection.append(
                    header(
                        str(elem) + ": " + str(data) + " " + unit_label,
                        sep=0.1,
                        style=ParaStyle,
                    )
                )
            else:
                print(elem, data, data.shape)
                metasection.append(
                    header(
                        str(elem) + ": " + str(data) + " " + unit_label,
                        sep=0.1,
                        style=ParaStyle,
                    )
                )

        src = KeepTogether(metasection)
        Elements.append(src)

        cxxtitle = header("Equations in CXX")
        cxxsection = []
        # print self.cxx
        cxxsection.append(header(self.cxx[0], sep=0.1, style=ParaStyle))
        cxxsrc = KeepTogether(cxxsection)

        Elements.append(cxxsrc)
        # for i,elem in enumerate(self.cxx):
        #    if type(self.meta[elem])== type({}):
        #       print elem #np.array(self.meta[elem]['v']).shape()
        #       if np.array(self.meta[elem]['v']).shape == (self.nx,self.ny):
        #          datastr = str(self.meta[elem]['v'][:,self.ny/2])
        #          metasection.append(graphout('stuff',
        #                                      self.meta[elem]['v'][:,self.ny/2]))
        #       else:
        #          datastr = str(self.meta[elem]['v'])
        #       metasection.append(header(str(elem)+': '+datastr
        #                                    + ' '+ str(self.meta[elem]['u']),
        #                                    sep=0.1, style=ParaStyle))

        if debug:
            return Elements
        go()

        output = PdfFileWriter()
        metapdf = PdfFileReader(file("meta.pdf", "rb"))
        mainpdf = PdfFileReader(file("output.pdf", "rb"))

        for i in range(0, metapdf.getNumPages()):
            output.addPage(metapdf.getPage(i))

        for i in range(0, mainpdf.getNumPages()):
            output.addPage(mainpdf.getPage(i))

        outputFile = filename
        outputStream = file(outputFile, "wb")
        output.write(outputStream)
        outputStream.close()
        print("Consolidation complete.")


class subset(LinResDraw):
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
