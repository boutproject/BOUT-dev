from __future__ import print_function
from builtins import range

###
# compute average growth rate bout variable f and plane y
# prints value in plane y and total average
# optional tind excludes initial 'tind' time steps
# Note it masks the values != Inf
###
import numpy as np
from boutdata.collect import collect
from boututils.moment_xyzt import moment_xyzt


def avgrate(p, y=None, tind=None):

    if tind is None:
        tind = 0

    rmsp_f = moment_xyzt(p, "RMS").rms

    ni = np.shape(rmsp_f)[1]
    nj = np.shape(rmsp_f)[2]

    growth = np.zeros((ni, nj))

    with np.errstate(divide="ignore"):

        for i in range(ni):
            for j in range(nj):
                growth[i, j] = np.gradient(np.log(rmsp_f[tind::, i, j]))[-1]

    d = np.ma.masked_array(growth, np.isnan(growth))

    # masked arrays
    # http://stackoverflow.com/questions/5480694/numpy-calculate-averages-with-nans-removed

    print("Total average growth rate= ", np.mean(np.ma.masked_array(d, np.isinf(d))))
    if y is not None:
        print(
            "Growth rate in plane",
            y,
            "= ",
            np.mean(np.ma.masked_array(growth[:, y], np.isnan(growth[:, y]))),
        )


# test
if __name__ == "__main__":
    path = "/Users/brey/BOUT/bout/examples/elm-pb/data"

    data = collect("P", path=path)

    avgrate(data, 32)
