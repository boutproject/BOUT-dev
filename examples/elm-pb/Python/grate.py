from __future__ import print_function
import numpy as np
from boutdata import collect
from boututils import moment_xyzt

path='../data/'

p=collect('P',path=path)
rmsp_f=moment_xyzt(p[:,34:35,32:33,:], 'RMS')
print(np.gradient(np.log(rmsp_f[:,0,0]))[-1])
