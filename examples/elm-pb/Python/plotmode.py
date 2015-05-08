from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div

from numpy import *;
#from scipy.io import readsav;
import matplotlib.pyplot as plt;

# Dynamic matplotlib settings
from matplotlib import rcParams;
rcParams['font.size'] = 20;
rcParams['legend.fontsize'] = 'small';
rcParams['legend.labelspacing'] = 0.1;
rcParams['lines.linewidth'] = 2;
rcParams['savefig.bbox'] = 'tight';

# Create image directory if not exists
import os;
if not os.path.exists('image'):
   os.makedirs('image');

#fphi = transpose(readsav('fphi.idl.dat')['fphi'])[:,:,:,];
fphi = load('fp.npy')

plt.figure();
for i in range(1, 9):
   print("Growth rate for mode number", i)
   print(gradient(log(abs(fphi[34, 32, i, :]))))
   plt.semilogy(((abs(fphi[34, 32, i, :]))), label = 'n=' + str(i * 5));

plt.legend(loc=2);
plt.xlabel('Time');
plt.savefig('image/plotmode.png');
plt.savefig('image/plotmode.eps');


plt.show(block=False);
plt.figure();
for i in range(1, 9):
   plt.plot(abs(fphi[:, 32, i, -1]), label = 'n=' + str(i * 5));

plt.legend();
plt.xlabel('X index');

plt.savefig('image/plotmodeamp.png');
plt.savefig('image/plotmodeamp.eps');

plt.show(block=False);

plt.figure();
for i in range(1, 9):
   plt.plot(old_div(abs(fphi[:, 32, i, -1]),abs(fphi[:, 32, i, -1]).max()), label = 'n=' + str(i * 5));

plt.legend();
plt.xlabel('X index');

plt.savefig('image/plotmodenorm.png');
plt.savefig('image/plotmodenorm.eps');

plt.show();

