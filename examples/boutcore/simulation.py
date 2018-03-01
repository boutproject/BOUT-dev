#!/bin/python3
import boutcore

boutcore.init("mesh:n=2")
mesh=boutcore.Mesh.getGlobal();
dens=boutcore.Field3D.fromMesh(mesh)
#dens.set(0)

def rhs(time):
    n_ddt=dens.ddt()
    n_ddt.set(1)

model=boutcore.PhysicsModel()
model.setRhs(rhs)
model.solve_for(n=dens)
model.solve()

