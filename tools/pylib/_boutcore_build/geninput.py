#!/usr/bin/env python3

import boutcore as bc
import numpy as np
bc.init("-d input -f ../data/BOUT.inp".split(" "))

f=bc.Field3D.fromMesh(None)
f.setAll(np.array([[[1.]]]))
f2=bc.Field3D.fromMesh(None)
f2.setAll(np.array([[[2.]]]))
print(f.getAll())
dump=bc.Datafile() # get global :(
dump.add(f3d=f,f2d=f2,save_repeat=True)
dump.write()

