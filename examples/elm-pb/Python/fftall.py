from __future__ import print_function
from numpy import *
from scipy.io import readsav

print('Calculating P..')
a=transpose(readsav('phi.idl.dat')['phi'])
fa=fft.fft(a,axis=2)
save('fp',fa)

