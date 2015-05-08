from builtins import range
import numpy
# Find the closest contour line to a given point
def closest_line(n, x, y, ri, zi, mind=None):
    
    mind = numpy.min( (x[0] - ri)**2 + (y[0] - zi)**2 )
    ind = 0
    
    for i in range (1, n) :
        d = numpy.min( (x[i] - ri)**2 + (y[i] - zi)**2 )
        if d < mind :
            mind = d
            ind = i
    return ind 
