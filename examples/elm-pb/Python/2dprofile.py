from __future__ import division
from builtins import range
from past.utils import old_div
import numpy as np
import matplotlib.pyplot as plt
from boututils.datafile import DataFile
from matplotlib.ticker import FixedFormatter, FormatStrFormatter, AutoLocator, AutoMinorLocator

with DataFile("./cbm18_dens8.grid_nx68ny64.nc") as f:
    g = {v: f.read(v) for v in f.keys()}

majorLocator   = AutoLocator()
majorFormatter = FormatStrFormatter('%3.0e')
minorLocator   = AutoMinorLocator()
Fm = FixedFormatter(['0','$1 \\times 10^4$','$2 \\times  10^4$','$3 \\times 10^4$','$4 \\times 10^4$'])
Fm2 = FixedFormatter(['0','$2 \\times 10^5$','$4 \\times  10^5$','$6 \\times 10^5$'])

bxy=g.get('Bxy')
p=g.get('pressure')
jpar0=g.get('Jpar0')
psixy=g.get('psixy')
btxy=g.get('Btxy')
shiftangle=g.get('ShiftAngle')

nx=g.get('nx')
ny=g.get('ny')

q = np.zeros((nx, ny))
for i in range(ny):
    q[:,i] = old_div(- shiftangle, (2 * np.pi))




xarr = psixy[:,0]
xarr = old_div((xarr + 0.854856), (0.854856 + 0.0760856))


fig=plt.figure()
plt.plot(xarr,q[:,32])
plt.xlabel('normalized $\psi$', fontsize=25)
plt.ylabel('$q$',rotation='horizontal',fontsize=25)

fig.set_tight_layout(True)

fig, ax1 = plt.subplots()

ax1.plot(xarr, p[:,32], 'r-', markevery=1, linewidth=3)
ax1.set_xlabel('normalized $\psi$',fontsize=25)
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('Pressure [Pa]', color='k',fontsize=25)

#set y limit
ax1.set_ylim(0,40000,10000)

#define ticks#
ax1.yaxis.set_ticks(np.arange(0, 40000, 10000))
#ax1.yaxis.set_major_locator(majorLocator)
#ax1.yaxis.set_major_formatter(majorFormatter)
ax1.yaxis.set_major_formatter(Fm)

#for the minor ticks, use no labels; default NullFormatter
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator(10))

#format tick labels
for tl in ax1.get_yticklabels():
    tl.set_color('k')


ax2 = ax1.twinx()

s2 = -jpar0

ax2.plot(xarr, s2[:,32], 'r-',markevery=1,linewidth=3)
ax2.set_ylabel('$J_\parallel  [A/m^2]$', color='k',fontsize=25)
ax2.set_ylim(0,600000)
ax2.yaxis.set_ticks(np.arange(0, 600000, 200000))
ax2.yaxis.set_major_formatter(Fm2)
for tl in ax2.get_yticklabels():
    tl.set_color('k')


fig.set_tight_layout(True)

plt.show()


#plt.savefig('2d.png', transparent=True)
