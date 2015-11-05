from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
###
# computes average growth rate for all points at the final timestep
# computes average growth rate for points in the mead plane at the final timestep
###
import numpy as np
from boutdata.collect import collect
from boututils.moment_xyzt import moment_xyzt

path='./data/'

p=collect('P',path=path)

nmpy=old_div(p.shape[2],2)  # define mead plane

ik = 50 # disregard the first ik timesteps

def gr(p):
        rmsp_f=moment_xyzt(p, 'RMS').rms

        ni=np.shape(rmsp_f)[1]
        nj=np.shape(rmsp_f)[2]

        growth=np.zeros((ni,nj))

        for i in range(ni):
                for j in range(nj):
                         growth[i,j]=np.gradient(np.log(rmsp_f[ik::,i,j]))[-1]

        return growth


growth=gr(p)

d=np.ma.masked_array(growth,np.isnan(growth))

# masked arrays
# http://stackoverflow.com/questions/5480694/numpy-calculate-averages-with-nans-removed

print('Total mean value = ', np.mean(np.ma.masked_array(d,np.isinf(d))))
mm=np.ma.masked_array(growth[:,nmpy],np.isnan(growth[:,nmpy]))
if np.isinf(np.mean(mm)) :
    print('There is an Inf value in the mead plane')
    print('Mean value of floating numbers in mead plane is = ', np.mean(np.ma.masked_array(mm,np.isinf(mm))))
else:
    print('Mean value in mead plane= ', np.mean(mm))
