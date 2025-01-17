#!/usr/bin/env python3
from __future__ import division
from past.utils import old_div

import matplotlib.pyplot as plt
import numpy as np
from boututils.moment_xyzt import moment_xyzt
from boututils.datafile import DataFile
from boutdata.collect import collect
import os
from pathlib import Path

#Dynamic matplotlib settings
from matplotlib import rcParams
rcParams['font.size'] = 20.
rcParams['legend.fontsize'] = 'small'
rcParams['lines.linewidth'] = 2

if not os.path.exists('image'):
   os.makedirs('image')
filename = Path(__file__).with_name("cbm18_dens8.grid_nx68ny64.nc")
with DataFile(str(filename)) as f:
    g = {v: f.read(v) for v in f.keys()}

psi = old_div((g['psixy'][:, 32] - g['psi_axis']), (g['psi_bndry'] - g['psi_axis']))

path = './data'

plt.figure()

p0=collect('P0', path=path)

p=collect('P', path=path)
res = moment_xyzt(p,'RMS','DC')
rmsp = res.rms
dcp = res.dc
nt = dcp.shape[0]

plt.plot(psi, p0[:, 32], 'k--', label='t=0')
plt.plot(psi, p0[:, 32] + dcp[nt//4, :, 32], 'r-', label='t='+np.str(nt//4))
plt.plot(psi, p0[:, 32] + dcp[nt//2, :, 32], 'g-', label='t='+np.str(nt//2))
plt.plot(psi, p0[:, 32] + dcp[3*nt//4, :, 32], 'b-', label='t='+np.str(3*nt//4))
plt.plot(psi, p0[:, 32] + dcp[-1, :, 32], 'c-', label='t='+np.str(nt))

plt.legend()
#plt.xlim(0.6, 1.0)
plt.xlabel(r'Normalized poloidal flux ($\psi$)')
plt.ylabel(r'$\langle p\rangle_\xi$')
plt.title(r'Pressure')
xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()


#plt.savefig('image/plotcollapse.png', bbox_inches='tight')
#plt.savefig('image/plotcollapse.eps', bbox_inches='tight')


plt.tight_layout()

print("Showing plot")
plt.show()
