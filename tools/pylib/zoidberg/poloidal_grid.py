"""
Routines for generating structured meshes on poloidal domains
"""

from numpy import pi, linspace, sqrt, cos, sin, append, zeros
from scipy.interpolate import splrep, splev, interp1d
from scipy.integrate import cumtrapz

import matplotlib.pyplot as plt

class RZline:
    """
    Represents (R,Z) coordinates of a periodic line
    
    Members
    -------

    R      Major radius [m]
    Z      Height [m]
    theta  Angle variable [radians]

    R,Z and theta all have the same length
 
    """
    def __init__(self, r, z):
        """
        r,z   1D arrays of the major radius (r) and height (z)
              which are of the same length. A periodic domain
              is assumed, so the last point connects to the first.

        Note that the last point in (r,z) arrays should not be the same
        as the first point. The (r,z) points are in [0,2pi)
        
        """
        self.R = r
        self.Z = z
        
        # Check the sizes of the variables
        n = len(r)
        assert len(z) == n
        assert r.shape == (n,)
        assert z.shape == (n,)
        
        # Define an angle variable
        self.theta = linspace(0,2*pi,n, endpoint=False)
        
        # Create a spline representation
        # Note that the last point needs to be passed but is not used
        self._rspl = splrep(append(self.theta,2*pi), append(r,r[0]), per=True)
        self._zspl = splrep(append(self.theta,2*pi), append(z,z[0]), per=True)
        
    def Rvalue(self, theta):
        """
        Calculate the value of R at given theta locations
        """
        return splev(theta, self._rspl, der=0)

    def Zvalue(self, theta):
        """
        Calculate the value of Z at given theta locations
        """
        return splev(theta, self._zspl, der=0)
        
    def dRdtheta(self):
        """
        Calculate the derivative of R w.r.t. theta
        """
        return splev(self.theta, self._rspl, der=1)

    def dZdtheta(self):
        """
        Calculate the derivative of Z w.r.t. theta
        """
        return splev(self.theta, self._zspl, der=1)
        
    def distance(self):
        """
        Integrates the distance along the line.
        
        Returns
        -------
        An array one longer than theta. The first element is zero,
        and the last element is the total distance around the loop
        """
        
        # variation of length with angle dl/dtheta
        dldtheta = sqrt(self.dRdtheta()**2 + self.dZdtheta()**2)
        
        # Make periodic by repeating the first point
        dldtheta = append(dldtheta, dldtheta[0])

        # Integrate cumulatively
        return cumtrapz(dldtheta, append(self.theta, 2.*pi), initial=0.0)

    def orderByDistance(self, n=None):
        """
        Returns a new RZline which has a theta uniform
        in distance along the line
        
        n   Number of points. Default is the same as the current line
        """
        if n is None:
            n = len(self.theta)
        
        # Distance along the line
        dist = self.distance()
        
        # Positions where points are desired
        positions = linspace(dist[0], dist[-1], n, endpoint=False)
        
        # Find which theta value these correspond to
        thetavals = interp1d(dist, append(self.theta, 2.*pi), copy=False, assume_sorted=True)
        new_theta = thetavals(positions)
        
        return RZline( self.Rvalue(new_theta), self.Zvalue(new_theta) )
        
def circular_boundaries(R0=1.0, rin=1.0, rout=2.0, n=20):
    """
    Creates a pair of RZline objects, for inner and outer boundaries

    >>> inner, outer = circular_boundaries()
    """
    # Define an angle coordinate
    theta = linspace(0,2*pi,n, endpoint=False)
    
    return RZline( R0 + rin*cos(theta), rin*sin(theta) ), RZline( R0 + rout*cos(theta), rout*sin(theta) )
    
    

def grid_annulus(inner, outer, nx, ny, show=True):
    """
    Grid an annular region, given inner and outer boundaries
    both of which are RZline objects
    """
    
    assert nx >= 2
    assert ny > 1
    
    R = zeros((nx, ny))
    Z = zeros((nx, ny))

    # Order both boundaries by distance
    inner = inner.orderByDistance()
    outer = outer.orderByDistance()
    
    # Generate angle values, which should now be equally spaced
    # in distance along inner and outer boundaries
    thetavals = linspace(0, 2*pi, ny, endpoint=False)
    
    # Radial coordinate
    xvals = linspace(0, 1.0, nx, endpoint=True)
    
    for i, x in enumerate(xvals):
        # Get the R and Z coordinates of this line
        R[i,:] = x*outer.Rvalue(thetavals) + (1.-x)*inner.Rvalue(thetavals)
        Z[i,:] = x*outer.Zvalue(thetavals) + (1.-x)*inner.Zvalue(thetavals)
        
    if show:
        plt.plot(inner.R, inner.Z, '-o')
        plt.plot(outer.R, outer.Z, '-o')
        
        plt.plot(R, Z, 'x')
        
        plt.show()
        

if __name__ == "__main__":
    inner,outer = circular_boundaries(R0=1.0, rin=1.0, rout=2.0, n=50)
    
    grid_annulus(inner, outer, 10, 20, show=True)
    
