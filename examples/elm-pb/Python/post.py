# -*- coding: utf-8 -*-
# Replicate the graphs from Xu_3fields_hands-on.pdf (BOUT++ workshop 2013)
# needs runs in data & data0 folders
# After running the cases define path0 and path1 accordingly
from __future__ import division
from past.utils import old_div
import numpy as np
from boututils.datafile import DataFile
from boutdata.collect import collect
from pylab import plot, show, figure, xlabel, ylabel, annotate, xlim, ylim
from boututils.moment_xyzt import moment_xyzt
from boututils.mode_structure import mode_structure
from boututils.plotpolslice import plotpolslice
from boututils.calculus import deriv
from mayavi import mlab


path0="./data0/"
path1="./data/"

period=15

gfile='./cbm18_dens8.grid_nx68ny64.nc'


with DataFile(gfile) as f:
    g = {v: f.read(v) for v in f.keys()}


Dphi0 = collect("Dphi0", path=path0)
phi0 = collect("phi0", path=path1) # needs diamagnetic effects
#
psixy=g.get('psixy')
PSI_AXIS=g.get('psi_axis')
PSI_BNDRY=g.get('psi_bndry')
#
psix=old_div((psixy[:,32]-PSI_AXIS),(PSI_BNDRY-PSI_AXIS))
Epsi=-deriv(phi0[:,32],psix)
#
#
fig=figure()
plot(psix,-Dphi0[:,32], 'r', linewidth=5)
plot(psix,Epsi,'k',linewidth=5)
annotate('w/o flow', xy=(.3, .7),  xycoords='axes fraction',horizontalalignment='center', verticalalignment='center', size=30)
annotate('w/ flow', xy=(.7, .4),  xycoords='axes fraction',horizontalalignment='center', verticalalignment='center', color='r', size=30)
xlabel('Radial $\psi$',fontsize=25)
ylabel('$\Omega(\psi)/\omega_A$',fontsize=25)
ylim([-.05,0])
xlim([0.4,1.2])
fig.set_tight_layout(True)
show(block=False)

p_f0 = collect("P", path=path0)
p_f = collect("P", path=path1)
#
rmsp_f0=moment_xyzt(p_f0, 'RMS').rms
rmsp_f=moment_xyzt(p_f, 'RMS').rms
#

fig=figure(figsize=(10, 8))
plot(np.gradient(np.log(rmsp_f0[:,34,32])), color='k',linewidth=3)
plot(np.gradient(np.log(rmsp_f[:,34,32])),color='red',linewidth=3)

ylabel('$\gamma / \omega_A$',fontsize=25)
xlabel('Time$(\\tau_A)$',fontsize=25)
annotate('w/o flow', xy=(.5, .7),  xycoords='axes fraction',horizontalalignment='center', verticalalignment='center', size=30)
annotate('w/ flow', xy=(.5, .4),  xycoords='axes fraction',horizontalalignment='center', verticalalignment='center', color='r', size=30)
ylim([0,0.5])
xlim([0,100])
fig.set_tight_layout(True)
show(block=False)



plotpolslice(p_f0[50,:,:,:],gfile,period=period, fig=1)
mlab.text(.01,.99,"w/o flow")

plotpolslice(p_f[50,:,:,:],gfile,period=period, fig=1)
mlab.text(.01,.99,"w/ flow")

fig=figure()
mode_structure(p_f0[50,:,:,:], g, period=period)
plot([40,40],[0,.014],'k--',linewidth=5)
annotate('w/o flow', xy=(.3, .7),  xycoords='axes fraction',horizontalalignment='center', verticalalignment='center', size=30)
ylim([0,0.014])
xlim([0,80])
fig.set_tight_layout(True)
show(block=False)


figure()
mode_structure(p_f[50,:,:,:], g, period=period)
plot([40,40],[0,.014],'k--',linewidth=5)
annotate('w/ flow', xy=(.3, .7),  xycoords='axes fraction',horizontalalignment='center', verticalalignment='center', color='k', size=30)
ylim([0,0.0001])
xlim([0,80])
show(block=False)
show()
