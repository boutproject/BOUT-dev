import boutcore as bc
import sys
from boutdata.collect import collect
import numpy as np

bc.init(sys.argv[1:])


opt = bc.Options()
i = 0
for mesh in [bc.Mesh.getGlobal(), bc.Mesh(section="mesh2")]:
    infuncs = [bc.create3D(f"n{k}:function", mesh=mesh) for k in range(4)]
    dump = bc.Datafile.new(opt, mesh, asGlobal=True)
    dump.openw(f"run{i}.dmp.nc")
    i += 1
    out = bc.create3D("0", mesh=mesh)
    dump.add(save_repeat=True, n=out)
    for useFFT in [True, False]:
        for f in infuncs:
            out.set(bc.Delp2(f, "CELL_DEFAULT", useFFT).get())
            dump.write()
