#!/usr/bin/env python3

import boutpp


def test_boutpp_simple_model(run_isolated):

    if run_isolated():
        return

    boutpp.init("-d mini")

    class MyModel(boutpp.PhysicsModel):
        def init(self, restart):
            self.n = boutpp.create3D("dens:function")
            self.solve_for(dens=self.n)

        def rhs(self, time):
            self.n.ddt(boutpp.DDX(self.n))

    model = MyModel()
    model.solve()
