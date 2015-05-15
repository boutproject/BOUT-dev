from __future__ import division
from past.utils import old_div
#
# Perform a linear regression fit
#

from numpy import mean

def linear_regression(x, y):
    """ Simple linear regression of two variables
    
    y = a + bx

    a, b = linear_regression(x, y)
    
    """
    
    if x.size != y.size:
        raise ValueError("x and y inputs must be the same size")
    
    mx = mean(x)
    my = mean(y)

    b = old_div((mean(x*y) - mx*my), (mean(x**2) - mx**2))
    a = my - b*mx
    
    return a, b

