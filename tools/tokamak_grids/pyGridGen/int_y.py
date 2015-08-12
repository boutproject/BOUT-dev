from __future__ import print_function
import numpy
from gen_surface import gen_surface
import copy
from boututils import int_func

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

        f[xi,yi] = int_func(var[xi,yi], simple=simple)
                
        if nosmooth==None :
            print('no smooth yet')
            #f[xi,yi] = SMOOTH(SMOOTH(f[xi,yi], 5, /edge_truncate), 5, /edge_truncate)
          
        loop[xi] = f[xi,yi[numpy.size(yi)-1]] - f[xi,yi[0]]
        
        if last == 1 : break
                   
    if nr == 1 :   
            return f, loop  
    else:
        return f
    


