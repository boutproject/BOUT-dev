#!/usr/bin/env python

# 
# Run the test, compare results against the benchmark
#

# Variables to compare
from __future__ import print_function
#from builtins import str
vars = [['n',2.e-4],
        ['max_error2',2.e-5],
        ['max_error3',2.e-4],
        ['max_error4',1.e-5],
	['max_error5',2.e-4],
	['max_error6',2.e-5],
	['max_error7',2.e-4],
	['max_error8',2.e-5]]  
#tol = 1e-4                  # Absolute (?) tolerance

from boututils.run_wrapper import shell, launch, getmpirun
from boutdata.collect import collect
#import numpy as np
from sys import stdout, exit
import numpy as np

##temps
import matplotlib.pyplot as plt


MPIRUN=getmpirun()

print("Making time subsampling test")
shell("make > make.log")

print("Running time subsampling test")
success = True

for nproc in [1]: #,2,4]:
  
  cmd = "./blob2d"
  
  shell("rm data/BOUT.dmp.*.nc")

  print("   %d processors...." % nproc)
  s, out = launch(cmd, runcmd=MPIRUN, nproc=nproc, pipe=True,verbose=True)
  f = open("run.log."+str(nproc), "w")
  f.write(out)
  f.close()

  n0 = collect("n", path="data", info=False)
  t  = collect("t_array", path="data", info=False)

  fig1 = plt.figure(num = 1, dpi=None, facecolor='w', edgecolor='k')
  plt.plot(t,n0[:,0,0,0],'.-') # plot an example

  for i in [2,4,8,16,32]:
      cmd = "./blob2d n:subsample="+str(i)
      s, out = launch(cmd, runcmd=MPIRUN, nproc=nproc, pipe=True,verbose=True)
      n = collect("n", path="data", info=False)
      t = collect("t_array", path="data", info=False)
      plt.plot(t[::i],n[:len(n[:,0,0,0])/i+1,0,0,0],'.-',label="subsample="+str(i))

      error = np.sum( abs( n[:len(n[:,0,0,0])/i+1,:,:,:] - n0[::i,:,:,:] ) )  
      #error = np.sum( abs( n[:len(n[:,0,0,0])/i+1,0,0,0] - n0[::i,0,0,0] ) )  # use 0,0,0 otherwise won't pass on multiple procs
      if error != 0.0:
        print("Fail, error is = "+str(error))
        success = False
      else:
        print("Pass")

  plt.legend(loc=0)
  fig1.savefig("test-subsample.pdf", dpi=1200, facecolor='w', edgecolor='w', orientation='portrait')

if success:
  print(" => All time subsampling tests passed")
  exit(0)
else:
  print(" => Some failed tests")
  exit(1)
