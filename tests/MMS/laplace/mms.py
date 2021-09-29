from __future__ import print_function
from builtins import str
from sympy import symbols, sin, diff

###


def Delp2(f):
    """Laplacian in X-Z"""
    d2fdx2 = diff(f, x, 2)
    d2fdz2 = diff(f, z, 2)

    return d2fdx2 + d2fdz2


###

# Define symbols

x = symbols("x")
z = symbols("z")
pi = symbols("pi")

# Define manufactured solution

solution = sin(x * pi)

# Calculate input

input = Delp2(solution)

# Print

print("solution = " + str(solution))
print("input = " + str(input))
