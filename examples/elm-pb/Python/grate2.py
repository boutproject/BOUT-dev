###
# compute average growth rate for all points
###
import numpy as np
from boutdata import collect
from boututils import moment_xyzt

path='./data/'

p=collect('P',path=path)


def gr(p):
	rmsp_f=moment_xyzt(p, 'RMS').rms

	ni=np.shape(rmsp_f)[1]
	nj=np.shape(rmsp_f)[2]

	growth=np.zeros((ni,nj))

	for i in range(ni):
    		for j in range(nj):
       			 growth[i,j]=np.gradient(np.log(rmsp_f[50::,i,j]))[-1]

	return growth


growth=gr(p)

d=np.ma.masked_array(growth,np.isnan(growth))

# masked arrays
# http://stackoverflow.com/questions/5480694/numpy-calculate-averages-with-nans-removed    
                        
print np.mean(np.ma.masked_array(d,np.isinf(d)))
print np.mean(np.ma.masked_array(growth[:,32],np.isnan(growth[:,32])))
