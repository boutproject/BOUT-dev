from __future__ import print_function
import numpy
from gen_surface import gen_surface
import copy
from boututils.int_func import int_func

 # Integrate a function over y
def int_y( var, mesh, loop=None, nosmooth=None, simple=None):

    nr=0
    if loop!=None : nr=1

    f = copy.deepcopy(var)

    s = numpy.shape(var)
    nx = s[0]
    loop = numpy.zeros(nx)

    status = gen_surface(mesh=mesh) # Start generator

    while True:

        period, yi, xi, last = gen_surface(period=None, last=None, xi=None )

        if period:
            # Add the first point onto the end to complete the loop
            yi = numpy.concatenate((yi,[yi[0]]))
        
        f1d = int_func(var[xi,yi], simple=simple)

        if nosmooth==None :
            print('no smooth yet')
            #f1d = SMOOTH(SMOOTH(f1d, 5, /edge_truncate), 5, /edge_truncate)

        loop[xi] = f1d[-1] - f1d[0]

        if period:
            # Remove the last point
            f1d = f1d[:-1]
            yi = yi[:-1]

        # Put data into result
        f[xi,yi] = f1d
        
        if last == 1 : break

    if nr == 1 :
            return f, loop
    else:
        return f
