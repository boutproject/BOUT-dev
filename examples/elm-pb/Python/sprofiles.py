from __future__ import division
from builtins import zip
from builtins import range
from past.utils import old_div
import numpy as np
from boututils.datafile import DataFile
from boututils.surface_average import surface_average
from boutdata.collect import collect
from pylab import plot, show, xlabel, ylabel, figure, legend, gca


path = "./data"

gfile = "../cbm18_dens8.grid_nx68ny64.nc"

with DataFile(gfile) as f:
    g = {v: f.read(v) for v in f.keys()}

var=collect("P", path=path)

sol=surface_average(var, g)
#sol=np.mean(var,axis=3)

p0av=collect("P0", path=path)

q=np.zeros(sol.shape)

for i in range(sol.shape[1]):
    q[:,i]=sol[:,i]+p0av[:,0]


psixy=g.get('psixy')
psi0=g.get('psi_axis')
psix=g.get('psi_bndry')

xarr = psixy[:,0]
xarr = old_div((xarr - psi0), (-psi0 + psix)) #for this grid


fig=figure()

nt=q.shape[1]

plot(xarr, p0av,'k',label='t=0')
plot(xarr,q[:,nt/4],'r',label='t='+np.str(nt/4))
plot(xarr,q[:,nt/2],'b',label='t='+np.str(nt/2))
plot(xarr,q[:,3*nt/4],'g',label='t='+np.str(3*nt/4))
plot(xarr, q[:,-1],'k',label='t='+np.str(nt))

from collections import OrderedDict
handles, labels = gca().get_legend_handles_labels()
by_label = OrderedDict(list(zip(labels, handles)))
legend(list(by_label.values()), list(by_label.keys()))


xlabel(r"$\psi$",fontsize=25)
ylabel(r"$2 \mu_0 <P^2> / B^2$",fontsize=25)

fig.set_tight_layout(True)

show()
