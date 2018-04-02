# BOUT++ Library - Write fluid simulations in curviilinear geometry
# Copyright (C) 2016, 2017, 2018 David Schwörer
#

from stencils import get_interp_vals
import numpy as np
x = np.linspace(-5, 5, 9)
xx = x * x

for f in [xx, x * x * x]:
    fn = xx[1::2]
    fs = xx[::2]
    for i in range(5):
        vals = get_interp_vals(4, i - 2)
        interp = np.dot(vals, fn)
        err = interp - fs[i]
        rel_err = err / np.max(fn)
        assert(rel_err < 1e-10)
        #print(vals)
