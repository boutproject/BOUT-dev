"""
Routines and classes for representing periodic lines in R-Z poloidal planes

"""
from numpy import pi, linspace, sqrt, cos, sin, append, zeros, argmin
from scipy.interpolate import splrep, splev, interp1d
from scipy.integrate import cumtrapz

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

    >>> inner = circle()
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

