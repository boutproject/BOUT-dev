from __future__ import division
from builtins import range
from past.utils import old_div
import numpy as np
from pylab import plot, show
from boutdata import collect

def rms(f):

    nt=f.shape[0]

    ns=f.shape[1]
    ne=f.shape[2]
    nz=f.shape[3]
    
    ar=np.zeros([nz])

    rms=np.zeros([nt,ns,ne])

    for i in range(nt):
        for j in range(ns):
            for k in range(ne):
                ar=f[i,j,k,:]
	        valav=np.sum(ar)
        	tot=np.sum(old_div(np.power(ar-valav,2),nz))
                rms[i,j,k]=np.sqrt(tot)
    return rms
        
def plot_rms(x,y):  
    s=plot(np.gradient(np.log(rmsp[:,x,y])))
    return s
    
#test

path='../data0'

data = collect("P", path=path)

rmsp=rms(data)
        
plot_rms(34,32)
show()