#!/usr/bin/env python3

#
# Run the test, compare results against the benchmark
#

from sys import stdout

import boutconfig
import pytest
from boutdata import collect
from boututils.run_wrapper import launch_safe, shell
from numpy import abs, array, log, max, mean, polyfit, sqrt

# Display the plots as well as saving to file
show_plot = False

# List of NX values to use
nxlist = [16, 32, 64, 128]

# Variables to compare
varlist = ["a", "b", "c"]
markers = ["bo", "r^", "kx"]
labels = [r"$" + var + r"$" for var in varlist]

methods = {
    "hermitespline": 3,
    "hermitesplinelegacy": 3,
    "hermitesplineserial": 3,
    "lagrange4pt": 3,
    "bilinear": 2,
}


@pytest.mark.parametrize("method", methods)
def test_interpolate(method):

    print("Running Interpolation test")
    success = True

    print("------------------------------")
    print(f"Using {method} interpolation")

    error_2 = {}
    error_inf = {}
    for var in varlist:
        error_2[var] = []  # L2 error (RMS)
        error_inf[var] = []  # Maximum error

    for nx in nxlist:
        dx = 1.0 / (nx)

        args = f" mesh:nx={nx + 4} mesh:dx={dx} MZ={nx} xzinterpolation:type={method}"
        nproc = 2 if method == "hermitespline" and boutconfig.has["petsc"] else 1
        args += f" NXPE={nproc}"

        cmd = "./test_interpolate" + args

        shell(["rm -rf data/BOUT.dmp.*"])

        _s, out = launch_safe(cmd, nproc=nproc, pipe=True)
        with open(f"run.log.{method}.{nx}", "w") as f:
            f.write(out)

        # Collect output data
        for var in varlist:
            interp = collect(var + "_interp", path="data", xguards=False, info=False)
            solution = collect(
                var + "_solution", path="data", xguards=False, info=False
            )

            E = interp - solution

            if False:
                import matplotlib.pyplot as plt

                def myplot(f, lbl=None):
                    plt.plot(f[:, 0, 6], label=lbl)

                myplot(interp, "interp")
                myplot(solution, "sol")
                plt.legend()
                plt.show()

            l2 = float(sqrt(mean(E**2)))
            linf = float(max(abs(E)))

            error_2[var].append(l2)
            error_inf[var].append(linf)

            print(f"{var:s} : l-2 {l2:.8f} l-inf {linf:.8f}")

    dx = 1.0 / array(nxlist)

    for var in varlist:
        fit = polyfit(log(dx), log(error_2[var]), 1)
        order = fit[0]
        stdout.write(f"{var:s} Convergence order = {order:.2f}")

        # Make sure scaling is at least 90% of expected order
        if order < 0.9 * methods[method]:
            print("............ FAIL")
            success = False
        else:
            print("............ PASS")

    if False:
        try:
            import matplotlib.pyplot as plt

            # Plot errors
            plt.figure()

            for var, mark, label in zip(varlist, markers, labels):
                plt.plot(dx, error_2[var], "-" + mark, label=label)
                plt.plot(dx, error_inf[var], "--" + mark)

            plt.legend(loc="upper left")
            plt.grid()

            plt.yscale("log")
            plt.xscale("log")

            plt.xlabel(r"Mesh spacing $\delta x$")
            plt.ylabel("Error norm")
            plt.title(f"Error scaling for {method}")

            name = f"error_scaling_{method}.pdf"
            plt.savefig(name)
            print(f"Plot saved to {name}")

            if show_plot:
                plt.show()
            plt.close()
        except ImportError:
            print("No matplotlib available")
    else:
        print("Plotting disabled")

    print("------------------------------")
    assert success, " => Some failed tests"
