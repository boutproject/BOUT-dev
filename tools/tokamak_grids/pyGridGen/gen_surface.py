from builtins import range
#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# Generator of continuous surfaces
#
# First call: 
#   status = gen_surface(mesh=mesh)     - Initialisation
# Subsequent calls
#   yi = gen_surface(period=period, last=last, xi=xi)
#
# period - Set to 1 if the surface is periodic, 0 otherwise
# last   - Set to 1 if this is the last surface
# xi     - The X index of this surface
#

import numpy
import sys


global m,ys,xind,nd,domain,visited


def range ( first, last ):
    if first < last :
        return  first + numpy.arange(last - first + 1)
    else:
        return last + numpy.arange(first - last + 1)[::-1]




def gen_surface ( mesh=None, period=None, last=None, xi=None):
    
    global m, ys, xind, nd, domain, visited
    
    if mesh != None :
        # initializing gen_surface 
        # the global patemeters are initialized if mesh is not None
        m = mesh
        xind = 0 # Radial surface
        nd = numpy.size(mesh.npol) # Number of domains
        domain = 0 # The domain to start in
    
        # Running total of npol to get starting y index
        ys = numpy.zeros(nd).astype(int)
        
        if numpy.size(ys) > 1 :
            for i in range (1, nd+1) : ys[i] = ys[i-1] + mesh.npol[i-1]
        
        # visited marks which domains have been used
        visited = numpy.zeros(nd).astype(int)
    
        return 0
    
    if xind >= numpy.sum(m.nrad) :
        last = 1
        return 1 # Error
   
     
    # Get the next surface
    ny = 0
    period = 0 # Mark as non-periodic
    last = 0   # mark as not the last
    xi = xind
    while True:
        if visited[domain] == 1 :
            # Already visited this domain
            period = 1 # Means this domain is periodic
            break
     
    
        # Get the range of indices for this domain        
        
        #yi = list(range(ys[domain], ys[domain]+m.npol-1)) # H.SETO (QST)
        yi = list(range(ys[domain], ys[domain]+m.npol[domain]-1))
        
        if ny == 0 :
            yinds = yi 
        else: yinds = [yinds, yi]
        
        #ny = ny + m.npol # H.SETO (QST)
        ny = ny + m.npol[domain]

        visited[domain] = 1 # Mark domain as visited
    
        # Find next domain
        if xind < m.yup_xsplit[domain] :
            domain = m.yup_xin[domain]
        else:
            domain = m.yup_xout[domain]

                                
        if domain < 0 :# Keep going until hit a boundary
            break
            
            
    # Find a domain which hasn't been visited
    w = numpy.size(numpy.where(visited == 0))
    
    if w != 0 :
        # See if there are any regions with boundaries on lower side
        domain = -1
        for i in range (w) :
            if xind < m.ydown_xsplit[w[i]] :
                d = m.ydown_xin[w[i]]
            else:
                d = m.ydown_xout[w[i]]
         
            if d < 0 :
                domain = w[i]
                break
       
     
        if domain < 0 : domain = w[0] # Set the domain to the first one
    
    else:
        # No domains left - increase x index (if possible)

        xind = xind + 1
        visited = numpy.zeros(nd).astype(int) # Set all to zeros again
        domain = 0 # Start again with the first domain
        if xind == numpy.sum(m.nrad) : last = 1 # No more left
     
  
    if ny == 0 : return 2 # This shouldn't happen
  
 
    return period, yinds, xi, last
 
