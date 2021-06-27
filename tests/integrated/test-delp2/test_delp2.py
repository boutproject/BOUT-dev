import boutcore as bc
import sys
from boutdata.collect import collect
import numpy as np

bc.init(sys.argv[2:])


class MyModel(bc.PhysicsModel):
    def __init__(self, useFFT, mesh):
        self.useFFT = useFFT
        self.mesh = mesh
        bc.PhysicsModel.__init__(self)

    def init(self, restart):
        opt = bc.Options("diffusion")
        self.D = opt.get("D", 0.1)
        self.n = bc.create3D("n:function", msh=self.mesh)
        self.solve_for(n=self.n)

    def rhs(self, time):
        self.mesh.communicate(self.n)

        self.n.ddt(self.D * bc.Delp2(self.n, "CELL_DEFAULT", self.useFFT))


opt = bc.Options()
i = 0
for mesh in [bc.Mesh.getGlobal(), bc.Mesh(section="mesh2")]:
    for useFFT in [True, False]:
        dump = bc.Datafile.new(opt, mesh, asGlobal=True)
        dump.openw(f"run{i}.dmp.nc")
        MyModel(useFFT, mesh).solve()
        i += 1
