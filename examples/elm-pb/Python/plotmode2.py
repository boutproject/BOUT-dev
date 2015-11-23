from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div

from numpy import *;
#from scipy.io import readsav;
import matplotlib.pyplot as plt;
from boutdata.collect import collect

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

path='./data/'
data=collect('P',path=path)

#fphi = transpose(readsav('fphi.idl.dat')['fphi'])[:,:,:,];
fphi = fft.fft(data, axis=3)

plt.figure();
for i in range(1, 9):
   print("Growth rate for mode number", i)
   print(gradient(log(abs(fphi[:,34, 32, i]))))
   plt.semilogy(((abs(fphi[:,34, 32, i]))), label = 'n=' + str(i * 5));

plt.legend(loc=2);
plt.xlabel('Time');
plt.savefig('image/plotmode.png');
plt.savefig('image/plotmode.eps');


plt.show(block=False);
plt.figure();
for i in range(1, 9):
   plt.plot(abs(fphi[-1, :, 32, i]), label = 'n=' + str(i * 5));

plt.legend();
plt.xlabel('X index');

plt.savefig('image/plotmodeamp.png');
plt.savefig('image/plotmodeamp.eps');

plt.show(block=False);

plt.figure();
for i in range(1, 9):
   plt.plot(old_div(abs(fphi[-1, :, 32, i]),abs(fphi[-1, :, 32, i]).max()), label = 'n=' + str(i * 5));

plt.legend();
plt.xlabel('X index');

plt.savefig('image/plotmodenorm.png');
plt.savefig('image/plotmodenorm.eps');

plt.show();
