"""
Routines for generating structured meshes on poloidal domains
"""

from numpy import pi, linspace, sqrt, cos, sin, append, zeros, argmin
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
        
    def Rvalue(self, theta=None, deriv=0):
        """
        Calculate the value of R at given theta locations
        """
        if theta is None:
            theta = self.theta
        return splev(theta, self._rspl, der=deriv)

    def Zvalue(self, theta=None, deriv=0):
        """
        Calculate the value of Z at given theta locations
        """
        if theta is None:
            theta = self.theta
        return splev(theta, self._zspl, der=deriv)
        

    def distance(self):
        """
        Integrates the distance along the line.
        
        Returns
        -------
        An array one longer than theta. The first element is zero,
        and the last element is the total distance around the loop
        """
        
        # variation of length with angle dl/dtheta
        dldtheta = sqrt(self.Rvalue(deriv=1)**2 + self.Zvalue(deriv=1)**2)
        
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
 
    def closestPoint(self, R,Z):
        """
        Find the closest point on the curve to the given (R,Z) point
        
        Returns the value of theta (angle)
        """
        
        # First find the closest control point
        ind = argmin((self.R - R)**2 + (self.Z - Z)**2)
        theta = self.theta[ind]

        # Refine using Newton's method to find where derivative goes to zero
        # distance = dR**2 + dZ**2
        # d(distance)/dtheta = 2.*dR*(dR/dtheta) + 2.*dZ*(dZ/dtheta)
        # d^2(distance)/dtheta^22 = 2*( (dR/dtheta)^2 + dR*(d^2R/dtheta^2) + (dZ/dtheta)^2 + dZ*(d^2Z/dtheta^2) )
        for i in range(3):
            dR = self.Rvalue(theta) - R
            dZ = self.Zvalue(theta) - Z
            
            dRdtheta = self.Rvalue(theta, deriv=1)
            dZdtheta = self.Zvalue(theta, deriv=1)
            
            d2Rdtheta2 = self.Rvalue(theta, deriv=2)
            d2Zdtheta2 = self.Zvalue(theta, deriv=2)

            d_dtheta = 2.*(dR*dRdtheta + dZ*dZdtheta)
            d2_dtheta2 = 2.*( dRdtheta**2 + dR*d2Rdtheta2 + dZdtheta**2 + dZ*d2Zdtheta2 )
            
            theta = theta - d_dtheta/d2_dtheta2
            
        return theta
        
       
def circle(R0=1.0, r= 0.5, n=20):
    """
    Creates a pair of RZline objects, for inner and outer boundaries

    >>> inner, outer = circular_boundaries()
    """
    # Define an angle coordinate
    theta = linspace(0,2*pi,n, endpoint=False)
    
    return RZline( R0 + r*cos(theta), r*sin(theta) )
    
def shaped_line(R0=3.0, a=1.0, elong=0.0, triang=0.0, indent=0.0, n=20):
    """
    Parametrisation of plasma shape from J. Manickam, Nucl. Fusion 24 595 (1984)
    
    R0  Major radius
    a   Minor radius
    elong   Elongation, 0 for a circle
    triang  Triangularity, 0 for a circle
    indent  Indentation, 0 for a circle
    """
    theta = linspace(0,2*pi,n, endpoint=False)
    return RZline( R0 - indent + (a + indent*cos(theta))*cos(theta + triang*sin(theta)),
                   (1.+elong)*a*sin(theta) )


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
        

def grid_weighted_distance(inner, outer, nx, ny, show=True, equidistant=False):
    """
    Grid region using weighted distance from inner and outer
    boundaries as a streamline function
    
    nx   Number of radial surfaces including edges. Must be >=2
    ny   Number of poloidal surfaces. Must be > 1

    equidistant  Poloidal points on each surface are equally spaced
    
    """
    
    assert nx >= 2
    assert ny > 1
    
    R = zeros((nx, ny))
    Z = zeros((nx, ny))
    
    # Order both boundaries by distance
    inner = inner.orderByDistance()
    outer = outer.orderByDistance()
    
    # Radial coordinate
    xvals = linspace(0, 1.0, nx, endpoint=True)
    
    # Generate angle values, which should now be equally spaced
    # in distance along inner and outer boundaries
    thetavals = linspace(0, 2*pi, ny, endpoint=False)
    
    plt.plot(inner.R, inner.Z)
    plt.plot(outer.R, outer.Z)
    
    # Inner boundary
    R[0,:] = inner.Rvalue(thetavals)
    Z[0,:] = inner.Zvalue(thetavals)

    for i in range(1,nx):
        for j, theta in enumerate(thetavals):
            # Start on inner boundary and find shortest path to outer boundary
            Rpos = R[i-1,j]
            Zpos = Z[i-1,j]
        
            plt.plot(Rpos, Zpos, 'x')
            
            # Find the closest location on the outer boundary
            outtheta = outer.closestPoint(Rpos, Zpos)
            outR = outer.Rvalue(outtheta)
            outZ = outer.Zvalue(outtheta)
            
            # Take a small step towards the outer boundary
            frac = 1./(nx-i) # Fraction of distance to move
            Rpos += frac*(outR - Rpos)
            Zpos += frac*(outZ - Zpos)
            
            R[i,j] = Rpos
            Z[i,j] = Zpos
        
        if equidistant:
            # Put all points on a surface into an RZline object
            surf = RZline(R[i,:], Z[i,:])
            # Remap so equidistant 
            surf = surf.orderByDistance()
            R[i,:] = surf.R
            Z[i,:] = surf.Z
        
        
    plt.plot(Rpos, Zpos, 'o')
    plt.show()
    
    

if __name__ == "__main__":
    inner = circle(R0=1.5, r=1.0, n=50)
    outer = circle(R0=1.0, r=2.0, n=50)

    #inner = shaped_line(R0=3.0, a=0.5, elong=1.0, triang=0.0, indent=1.0, n=50)
    #outer = shaped_line(R0=2.8, a=1.5, elong=1.0, triang=0.0, indent=0.2, n=50)
    #outer = shaped_line(R0=3.0, a=1.0, elong=1.0, triang=0.0, indent=1.0, n=50)
    
    grid_weighted_distance(inner, outer, 10, 20, show=True, equidistant=False)
        
    
