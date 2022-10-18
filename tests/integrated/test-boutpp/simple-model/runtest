#!/usr/bin/env python3

# requires boutpp
# requires not make

from boutpp import *

init("-d mini")


class MyModel(PhysicsModel):
    def init(self, restart):
        self.n = create3D("dens:function")
        self.solve_for(dens=self.n)

    def rhs(self, time):
        self.n.ddt(DDX(self.n))


model = MyModel()
model.solve()
