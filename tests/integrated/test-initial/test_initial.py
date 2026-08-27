#!/usr/bin/env python3

# Test initial conditions

import pytest
import configparser
import itertools
import platform
import warnings
import numpy as np
from scipy.special import erf
from pathlib import Path

from boututils.run_wrapper import shell, launch_safe
from boutdata.collect import collect

########################################
# Implementations of BOUT++ functions


def bout_round(x):
    """
    BOUT++ rounding
    """
    return x + 0.5 if x > 0.0 else x - 0.5


def genRand(seed):
    """
    BOUT++ psuedo random number generator

    This PRNG has no memory, i.e. you need to call it with a different
    seed each time
    """
    # Make sure seed is
    if seed < 0.0:
        seed *= -1

    # Round the seed to get the number of iterations
    niter = int(11 + (23 + bout_round(seed)) % 79)

    # Start x between 0 and 1
    A = 0.01
    B = 1.23456789
    x = (A + np.mod(seed, B)) / (B + 2.0 * A)

    # Iterate logistic map
    for i in range(niter):
        x = 3.99 * x * (1.0 - x)

    return x


def ballooning(x, ball_n=3):
    """
    Ballooning function. Currently too tricky to implement
    """
    raise NotImplementedError("ballooning")


def gauss(x, width=1.0):
    """
    Normalised gaussian
    """
    return np.exp(-(x**2) / (2 * width**2)) / np.sqrt(2 * np.pi)


def mixmode(x, seed=0.5):
    """
    14 modes with random phases
    """
    result = 0.0
    for i in range(14):
        phase = np.pi * (2.0 * genRand(seed + i) - 1.0)
        result += (1.0 / (1.0 + np.abs(i - 4.0)) ** 2) * np.cos(i * x + phase)
    return result


def heaviside(x):
    """
    Heaviside step function
    """
    return 1 * (x > 0)


def tanhhat(x, width, centre, steepness):
    """
    BOUT++ TanhHat function
    """
    return 0.5 * (
        np.tanh(steepness * (x - (centre - width / 2.0)))
        + np.tanh(-steepness * (x - (centre + width / 2.0)))
    )


def atan(x, y=None):
    """
    Resolves to either np.arctan or np.arctan depending on the number of arguments
    """
    if y is not None:
        return np.arctan2(x, y)
    else:
        return np.arctan(x)


def max(*args):
    """
    Maximum of *args at each point
    """
    current = args[0]
    for arg in args:
        current = np.maximum(arg, current)
    return current


def min(*args):
    """
    Minimum of *args at each point
    """
    current = args[0]
    for arg in args:
        current = np.minimum(arg, current)
    return current


def fmod(x, denominator=1.0):
    """
    Modulo operator using fmod convention, (rem in Matlab)
    """
    return np.fmod(x, denominator)


# Rename functions to match BOUT++ naming
# Mostly just alternative names to numpy functions
abs = np.abs
asin = np.arcsin
acos = np.arccos
ballooning = ballooning
cos = np.cos
cosh = np.cosh
exp = np.exp
tanh = np.tanh
H = heaviside
log = np.log
power = np.power
sin = np.sin
sinh = np.sinh
sqrt = np.sqrt
tan = np.tan
TanhHat = tanhhat
pi = np.pi
erf = erf


@pytest.mark.parametrize("nproc", [1, 2, 3, 4])
def test_initial(nproc):
    tolerance = 1e-13
    cmd = "./test_initial"
    datadir = Path("data")
    inputfile = datadir / "BOUT.inp"

    shell(["rm -f data/BOUT.dmp.*.nc"])

    # Read the input file
    config = configparser.ConfigParser()
    with open(inputfile, "r") as f:
        config.read_file(itertools.chain(["[global]"], f), source=str(inputfile))

    # Find the variables that have a "function" option
    varlist = [key for key, values in config.items() if "function" in values]

    # Remove the coordinate arrays
    for coord in ["var_x", "var_y", "var_z"]:
        varlist.remove(coord)

    _, out = launch_safe(cmd, nproc=nproc, pipe=True, verbose=True)
    with open(f"run.log.{nproc}", "w") as f:
        f.write(out)

    # Collect the coordinate arrays separately
    x = collect("var_x", xguards=True, yguards=True, path=str(datadir), info=False)
    y = collect("var_y", xguards=True, yguards=True, path=str(datadir), info=False)
    z = collect("var_z", xguards=True, yguards=True, path=str(datadir), info=False)

    failures = []

    # Evaluate the functions
    for var in varlist:
        function = config[var]["function"]
        function = function.replace("^", "**")

        if ":" in function:
            print(
                f"{var} contains reference to variable - not possible to resolve at this time"
            )
            continue

        context = {"x": x, "y": y, "z": z}

        try:
            analytic = eval(function, globals(), context)
        except NotImplementedError as err:
            print(f"{err.args[0]} not implemented, skipping")
            continue

        data = collect(var, xguards=True, yguards=True, path=str(datadir), info=False)
        E2 = np.sqrt(np.mean((analytic - data) ** 2))
        if E2 >= tolerance:
            if var in ("mixmode", "mixmode_seed") and E2 < 1e-3:
                arch = platform.machine()
                if arch == "i686":
                    # This can happen due tue excess precision e.g. on X87 architecture
                    warnings.warn(
                        f"WARNING: Excess precision detected on i686 architecture for {var} (E2={E2:.4e})"
                    )
                    continue
                else:
                    failures.append(
                        f"{var} (arch: {arch}): l-2 error {E2:.4e} >= {tolerance}"
                    )
            else:
                failures.append(f"{var}: l-2 error {E2:.4e} >= {tolerance}")

    # Assert that the failures list is empty, otherwise print all collected failures
    assert not failures, (
        "The following variables failed tolerance checks:\n" + "\n".join(failures)
    )
