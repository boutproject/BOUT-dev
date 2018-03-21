#!/usr/bin/python3

import boutcore as bc
import numpy as np
bc.init("-d input -f ../data/BOUT.inp".split(" "))

f=bc.Field3D.fromMesh(None)
zeros=np.zeros((6,6,2))
f.set(zeros+1.)
f2=bc.Field3D.fromMesh(None)
f2.set(zeros+2.)
dump=bc.Datafile() # get global :(
dump.add(f3d=f,f2d=f2,save_repeat=False)
dump.write()

