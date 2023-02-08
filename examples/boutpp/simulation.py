#!/bin/python3
import boutpp as bc

bc.init("mesh:n=48")

class Model(bc.PhysicsModel):
    def init(self,restart):
        self.dens = bc.create3D("sin(x)")
        self.solve_for(n=self.dens)

    def rhs(self,time):
        self.dens.ddt(bc.DDX(self.dens))


model = Model()
model.solve()
