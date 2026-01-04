#!/usr/bin/env python3
import boutpp
import sys

# requires boutpp
# requires not make


def test_boutpp_legacy_model(request):

    boutpp.init(sys.argv[1:])
    request.node.boutpp_initialized = True

    dens = boutpp.create3D("0")


    def rhs(time):
        n_ddt = dens.ddt()
        n_ddt[:, :, :] = dens * 0
        n_ddt += 1


    model = boutpp.PhysicsModelBase()
    model.setRhs(rhs)
    model.solve_for(n=dens)
    model.solve()
