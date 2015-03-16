from __future__ import division
from builtins import range
from past.utils import old_div
import numpy as np
from boututils import file_import, surface_average, showdata
from boutdata import collect
from pylab import plot, show, annotate, xlabel, ylabel, figure, rc, rcParams, arrow, xlim, ylim


path='../data/'
 
gfile='../cbm18_dens8.grid_nx68ny64.nc'
 
g = file_import(gfile)
   
var=collect("P", path=path)   
   
sol=surface_average(var, g)


p0=g.get('pressure')



p0av=np.mean(p0,axis=1)

MU0 = 4.0e-7*np.pi

b0=g.get('bmag')

p0av=p0av*2*MU0/b0**2

q=np.zeros(sol.shape)

for i in range(sol.shape[1]):
    q[:,i]=sol[:,i]+p0av
    
    
psixy=g.get('psixy')
psi0=g.get('psi_axis')
psix=g.get('psi_bndry')

xarr = psixy[:,0]
xarr = old_div((xarr - psi0), (-psi0 + psix)) #for this grid
    
    
fig=figure()    

plot(xarr, p0av,'k-',xarr,q[:,:50:5],'r--')


xlabel(r"$\psi$",fontsize=25)
ylabel(r"$2 \mu_0 <P^2> / B^2$",fontsize=25)

#xlim(.6,.9)
#ylim(0,.008)

fig.set_tight_layout(True)

show(block=False)
showdata(q.T, tslice=1)
