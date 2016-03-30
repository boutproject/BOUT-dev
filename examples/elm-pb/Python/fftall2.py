from __future__ import print_function
from numpy import *
from boutdata.collect import collect

path='./data/'
data=collect('P',path=path)

print('Saving P..')
fa=fft.fft(data,axis=3)
save('fp',rollaxis(fa,0,4))
