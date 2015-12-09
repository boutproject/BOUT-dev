from __future__ import print_function
######
# Computes the rms of a variable and prints the growth rate value at the last timestep
######

import numpy as np
from boutdata.collect import collect
from boututils.moment_xyzt import moment_xyzt



path='./data/'


p=collect('P',path=path)
rmsp_f=moment_xyzt(p[:,34:35,32:33,:], 'RMS').rms


print(np.gradient(np.log(rmsp_f[:,0,0]))[-1])
