#####################################
# Example of analysis of data
#####################################

import numpy as np
import pylab as plt
from boutdata.collect import collect

path='./data/'
var=collect('P', path=path)

dcvar=np.mean(var, axis=3)
rmsvar=np.sqrt(np.mean(var**2,axis=3)-dcvar**2)

plt.figure()
plt.plot(rmsvar[:,34,32])
plt.show(block=False)
fvar=np.fft.rfft(var,axis=3)

plt.figure()
plt.plot(abs(fvar[:,34,32,1:10]))
plt.show(block=False)


plt.figure()
plt.semilogy(abs(fvar[:,34,32,1:7]))
plt.show(block=False)

plt.show()
