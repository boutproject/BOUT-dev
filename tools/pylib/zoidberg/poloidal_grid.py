"""
Routines for generating structured meshes on poloidal domains

Classes
-------

RectangularPoloidalGrid   Simple rectangles in R-Z
StructuredPoloidalGrid    Curvilinear structured grids in R-Z

Functions
---------

grid_elliptic   Create a StructuredPoloidalGrid 
                from inner and outer RZLine objects
                using elliptic meshing method.

"""

import numpy as np
from numpy import pi, linspace, zeros
from scipy.interpolate import RectBivariateSpline
from scipy.spatial import cKDTree as KDTree

import matplotlib.pyplot as plt


class RectangularPoloidalGrid(object):
    """
    Represents a poloidal grid consisting of a rectangular domain
    
    Note: Here the 2D plane (R,Z) is labelled by (x,y) indices

    Members
    -------
    
    nx, ny  Number of points in x and y

    R  2D Numpy array of R coordinates
    Z  2D Numpy array of Z coordinates
    
    """
    
    def __init__(self, nx, ny, Lx, Ly, Rcentre=0.0, Zcentre=0.0):
        """
        Inputs
        ------
        
        nx  Number of points in major radius (including boundaries)
        ny  Number of points in height (including boundaries)
        Lx  Radial domain size  [m]
        Ly  Vertical domain size [m]
        
        Rmid 
        """
        
        self.nx = nx
        self.ny = ny

        self.Lx = Lx
        self.Ly = Ly
        
        self.Rcentre = Rcentre
        self.Zcentre = Zcentre

        # Some useful derived quantities
        self.dR = self.Lx/(self.nx-1)
        self.dZ = self.Ly/(self.ny-1)
        self.Rmin = self.Rcentre - 0.5*self.Lx
        self.Zmin = self.Zcentre - 0.5*self.Ly

        # Generate 2D arrays
        # Using getCoordinate to ensure consistency
        xind, yind = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
        self.R, self.Z = self.getCoordinate(xind, yind)
        
    def __repr__(self):
        return "RectangularPoloidalGrid({0},{1},{2},{3},Rcentre={4},Zcentre={5})".format(self.nx, self.ny, self.Lx, self.Ly, self.Rcentre, self.Zcentre)
        
    def getCoordinate(self, xind, yind, dx=0, dy=0):
        """
        Get coordinates (R,Z) at given (xind,yind) index
    
        Inputs
        ------

        xind, yind   Indices in X and Y. These should be the same shape
        dx  Order of x derivative
        dy  Order of y derivative
        
        Returns
        -------
        
        R, Z  Locations of point
        or derivatives of R,Z with respect to indices if dx,dy != 0
        
        """
        
        # Convert to NumPy arrays if not already
        xind = np.asfarray(xind)
        yind = np.asfarray(yind)
        # Make sure dx and dy are integers
        dx = int(dx)
        dy = int(dy)
        
        assert xind.shape == yind.shape
        assert dx >= 0
        assert dy >= 0
        
        shape = xind.shape
        
        if dx + dy > 2:
            # Second derivatives and mixed derivatives all zero
            return np.zeros(shape), np.zeros(shape)
        
        if dx == 1:
            # dR/dx, dZ/dx
            return np.full(shape, self.dR), np.zeros(shape)
        elif dy == 1:
            # dR/dy, dZ/dy
            return np.zeros(shape), np.full(shape, self.dZ)
        # Return (R,Z) location
        return self.Rmin + xind*self.dR,  self.Zmin + yind*self.dZ
        

    def findIndex(self, R, Z):
        """
        Finds the (x,y) index corresponding to the given (R,Z) coordinate
        
        Inputs
        ------

        R,Z    Locations. Can be scalar or array, must be the same shape 
        
        Returns
        -------
        
        x,y index as a float, same shape as R,Z

        """

        # Make sure inputs are NumPy arrays
        R = np.asfarray(R)
        Z = np.asfarray(Z)
        
        # Check that they have the same shape
        assert R.shape == Z.shape
        
        return (R - self.Rmin)/self.dR, (Z - self.Zmin)/self.dZ

class StructuredPoloidalGrid(object):
    """
    Represents a structured poloidal grid in R-Z

    Members
    -------

    nx, ny  Number of points in x and y
    R, Z    2D NumPy arrays (nx,ny) of coordinates
    
    """
    def __init__(self, R, Z):
        """
        
        Inputs
        ------
        
        R, Z   2D Numpy arrays of R,Z points. Must be the same shape

        Note: R,Z are not copied, so these arrays should not be modified afterwards
        """
        
        assert R.shape == Z.shape
        
        self.R = R
        self.Z = Z
        
        # Create a KDTree for quick lookup of nearest points
        n = R.size
        data = np.concatenate( (R.reshape((n,1)), Z.reshape((n,1)) ), axis=1)
        self.tree = KDTree(data)
        
        # Create splines for quick interpolation of coordinates
        nx,ny = R.shape

        self.nx = nx
        self.ny = ny
        
        xinds = np.arange(nx)
        yinds = np.arange(ny+1)
        # Repeat the final point in y since periodic in y
        R_ext = np.concatenate((R, np.reshape(R[:,0], (nx, 1))), axis=1)
        Z_ext = np.concatenate((Z, np.reshape(Z[:,0], (nx, 1))), axis=1)
        
        self._spl_r = RectBivariateSpline(xinds, yinds, R_ext)
        self._spl_z = RectBivariateSpline(xinds, yinds, Z_ext)
        
    def __repr__(self):
        return "StructuredPoloidalGrid()"

    def getCoordinate(self, xind, yind, dx=0, dy=0):
        """
        Get coordinates (R,Z) at given (xind,yind) index
    
        Inputs
        ------

        xind, yind   Indices in X and Y. These should be the same shape
        dx  Order of x derivative
        dy  Order of y derivative
        
        Returns
        -------
        
        R, Z  Locations of point
        or derivatives of R,Z with respect to indices if dx,dy != 0
        
        """
        nx,ny = self.R.shape
        if (np.amin(xind) < 0) or (np.amax(xind) > nx-1):
            raise ValueError("x index out of range")
        if (np.amin(yind) < 0) or (np.amax(yind) > ny-1):
            raise ValueError("y index out of range")
            
        R = self._spl_r(xind, yind, dx=dx, dy=dy, grid=False)
        Z = self._spl_z(xind, yind, dx=dx, dy=dy, grid=False)
        
        return R,Z
        
    def findIndex(self, R, Z, tol=1e-10, show=False):
        """
        Finds the (x,y) index corresponding to the given (R,Z) coordinate
        
        Inputs
        ------

        R,Z    Locations. Can be scalar or array, must be the same shape 
        tol    Maximum tolerance on the square distance
        
        Returns
        -------
        
        x,y index as a float, same shape as R,Z

        """

        # Make sure inputs are NumPy arrays
        R = np.asfarray(R)
        Z = np.asfarray(Z)
        
        # Check that they have the same shape
        assert R.shape == Z.shape
        
        input_shape = R.shape  # So output has same shape as input
        
        # Get distance and index into flattened data
        # Note ind can be an integer, or an array of ints
        # with the same number of elements as the input (R,Z) arrays
        n = R.size
        position = np.concatenate( (R.reshape((n,1)), Z.reshape((n,1)) ), axis=1)
        dists, ind = self.tree.query(position)
        
        # Calculate (x,y) index
        nx,ny = self.R.shape
        xind = np.floor_divide(ind, ny)
        yind = ind - xind*ny
        
        # Convert indices to float
        xind = np.asfarray(xind)
        yind = np.asfarray(yind)
        
        if show:
            plt.plot(self.R, self.Z, '.')
            plt.plot(R, Z, 'x')
        
        while True:
            # Use Newton iteration to find the index
            # dR, dZ are the distance away from the desired point
            Rpos,Zpos = self.getCoordinate(xind, yind)
            if show:
                plt.plot(Rpos, Zpos, 'o')
            dR = Rpos - R
            dZ = Zpos - Z
            
            # Check if close enough
            if np.amax(dR**2 + dZ**2) < tol:
                break
            
            # Calculate derivatives
            dRdx, dZdx = self.getCoordinate(xind, yind, dx=1)
            dRdy, dZdy = self.getCoordinate(xind, yind, dy=1)
        
            # Invert 2x2 matrix to get change in coordinates
            #
            # (x) -=  ( dR/dx   dR/dy )^-1  (dR) 
            # (y)     ( dZ/dx   dZ/dy )     (dz) 
            #
            #
            # (x) -=  ( dZ/dy  -dR/dy ) (dR) 
            # (y)     (-dZ/dx   dR/dx ) (dZ) / (dR/dx*dZ/dy - dR/dy*dZ/dx)
            determinant = dRdx*dZdy - dRdy*dZdx
            
            xind -= (dZdy*dR - dRdy*dZ) / determinant
            yind -= (dRdx*dZ - dZdx*dR) / determinant
            
        if show:
            plt.show()
            
        return xind.reshape(input_shape), yind.reshape(input_shape)

def grid_annulus(inner, outer, nx, ny, show=True, return_coords=False):
    """
    Grid an annular region, given inner and outer boundaries
    both of which are RZline objects
    
    This is a very simple algorithm which just draws straight lines
    between inner and outer boundaries.
    """
    
    assert nx >= 2
    assert ny > 1
    
    R = zeros((nx, ny))
    Z = zeros((nx, ny))
    
    # Generate angle values, which should now be equally spaced
    # in distance along inner and outer boundaries
    thetavals = linspace(0, 2*pi, ny, endpoint=False)
    
    # Radial coordinate
    xvals = linspace(0, 1.0, nx, endpoint=True)
    
    innerR = inner.Rvalue(thetavals)
    innerZ = inner.Zvalue(thetavals)

    outerR = outer.Rvalue(thetavals)
    outerZ = outer.Zvalue(thetavals)
    for i, x in enumerate(xvals):
        # Get the R and Z coordinates of this line
        R[i,:] = x*outerR + (1.-x)*innerR
        Z[i,:] = x*outerZ + (1.-x)*innerZ
        
    if show:
        plt.plot(inner.R, inner.Z, '-o')
        plt.plot(outer.R, outer.Z, '-o')
        
        plt.plot(R, Z, 'x')
        
        plt.show()

    if return_coords:
        return R, Z
    return StructuredPoloidalGrid(R,Z)
        
def grid_elliptic(inner, outer, nx, ny, show=False, tol=1e-10, restrict_size=20, restrict_factor=2, return_coords=False):
    """
    Create a structured grid between inner and outer boundaries
    using elliptic method
    
    Input
    -----

    inner, outer   RZline objects describing inner and outer domain boundaries
    nx             The required radial resolution, including boundaries
    ny             The required poloidal resolution
    show           Display plots of intermediate results
    tol            Controls when iteration stops
    restrict_size    The size (nx or ny) above which the grid is coarsened
    restrict_factor  The factor by which the grid is divided if coarsened
    

    Returns
    -------
    If return_coords is true, returns R,Z as arrays. 
    If return_coords is false, returns a StructuredPoloidalGrid object

    Details
    -------

    Coordinates x = x(R, Z) and y = y(R,Z)
    obey an elliptic equation

    d^2x/dR^2 + d^2x/dZ^2 = 0
    
    d^2y/dR^2 + d^2y/dZ^2 = 0
    
    where here x is in in the domain (0,1) and y in (0,2pi)
    
    The above equations are inverted, giving:

    a*R_xx - 2*b*R_xy + c*R+yy = 0
    
    a*Z_xx - 2*b*Z_xy + c*Z_yy = 0

    where 

    a = R_y^2 + Z_y^2
    b = R_y*R_x + Z_x*Z_y
    c = R_x^2 + Z_x^2
    
    This is a nonlinear system of equations which is solved
    iteratively.
    
    See:
    https://www.nada.kth.se/kurser/kth/2D1263/l2.pdf
    https://en.wikipedia.org/wiki/Principles_of_grid_generation
    """
    
    assert nx >= 2
    assert ny > 1
    
    # Generate angle values (y coordinate), 
    # which should now be equally spaced
    # in distance along inner and outer boundaries
    thetavals = linspace(0, 2*pi, ny, endpoint=False)

    # Radial coordinate
    xvals = linspace(0, 1.0, nx, endpoint=True)

    if (nx > restrict_size) or (ny > restrict_size):
        # Create a coarse grid first to get a starting guess
        # Only restrict the dimensions which exceed restrict_size
        # Note that this might result in multiple levels of resolution
        
        nx_r = nx
        if nx > restrict_size:
            nx_r = int(nx / restrict_factor)
            
        ny_r = ny
        if ny > restrict_size:
            ny_r = int(ny / restrict_factor)

        # Create the coarse mesh
        R_r, Z_r = grid_elliptic(inner, outer, nx_r, ny_r, 
                                 tol=tol, restrict_size=restrict_size, restrict_factor=restrict_factor, return_coords=True)

        y_r = linspace(0, 2*pi, ny_r+1, endpoint=True) # Add on the final point duplicating the first
        x_r = linspace(0, 1.0, nx_r, endpoint=True)
        
        R_r = np.concatenate((R_r, np.reshape(R_r[:,0], (nx_r, 1))), axis=1)
        Z_r = np.concatenate((Z_r, np.reshape(Z_r[:,0], (nx_r, 1))), axis=1)
        
        # Now interpolate
        spl = RectBivariateSpline(x_r, y_r, R_r)
        R = spl(xvals, thetavals, grid=True)
        spl = RectBivariateSpline(x_r, y_r, Z_r)
        Z = spl(xvals, thetavals, grid=True)
        
    else:
        # Interpolate coordinates of inner and outer boundary
        Rinner = inner.Rvalue(thetavals)
        Zinner = inner.Zvalue(thetavals)
        
        Router = outer.Rvalue(thetavals)
        Zouter = outer.Zvalue(thetavals)
        
        # Interpolate in x between inner and outer
        # to get starting guess for a grid
        R = zeros((nx, ny))
        Z = zeros((nx, ny))
        for i in range(nx):
            R[i,:] = xvals[i]*Router + (1.-xvals[i])*Rinner
            Z[i,:] = xvals[i]*Zouter + (1.-xvals[i])*Zinner
    
    dx = xvals[1] - xvals[0]
    dy = thetavals[1] - thetavals[0]

    if show:
        plt.plot(inner.R, inner.Z)
        plt.plot(outer.R, outer.Z)
        #plt.plot(R, Z, 'o') # Starting locations
        
    # Start solver loop
    while True:
        # Calculate coefficients, which exclude boundary points
        # Note that the domain is periodic in y so roll arrays
        
        R_xm = R[:-2,:]   # R(x-1,y)
        R_xp = R[2:, :]   # R(x+1,y)
        R_ym = np.roll(R, 1) # R(x, y-1)
        R_yp = np.roll(R,-1) # R(x, y+1)
        R_xmym = R_ym[:-2,:] # R(x-1, y-1)
        R_xpym = R_ym[2:,:]  # R(x+1, y-1)
        R_xmyp = R_yp[:-2,:] # R(x-1, y+1)
        R_xpyp = R_yp[2:,:]  # R(x+1, y+1)
        R_ym = R_ym[1:-1,:]  # Now chop off boundaries
        R_yp = R_yp[1:-1,:]  # This is done to minimise number of rolls
        
        Z_xm = Z[:-2,:]
        Z_xp = Z[2:, :]
        Z_ym = np.roll(Z, 1)
        Z_yp = np.roll(Z,-1)
        Z_xmym = Z_ym[:-2,:] 
        Z_xpym = Z_ym[2:,:]
        Z_xmyp = Z_yp[:-2,:]
        Z_xpyp = Z_yp[2:,:]
        Z_ym = Z_ym[1:-1,:]
        Z_yp = Z_yp[1:-1,:]
        

        dRdy = (R_yp - R_ym)/(2.*dy)
        dRdx = (R_xp - R_xm)/(2.*dx)
        
        dZdy = (Z_yp - Z_ym)/(2.*dy)
        dZdx = (Z_xp - Z_xm)/(2.*dx)
        
        a = dRdy**2 + dZdy**2
        b = dRdy*dRdx + dZdx*dZdy
        c = dRdx**2 + dZdx**2
        
        # Now solve a*R_xx - 2*b*R_xy + c*R_yy = 0
        # For now using Jacobi update
    
        a_dx2 = a/dx**2
        b_dxdy = b/(2.*dx*dy)
        c_dy2 = c/dy**2
        inv_diag = 1./(2*a/dx**2 + 2*c/dy**2)
        
        Rold = R.copy()
        Zold = Z.copy()
        
        R[1:-1,:] = ( a_dx2*(R_xm + R_xp) - b_dxdy*(R_xpyp - R_xmyp - R_xpym + R_xmym) + c_dy2*(R_ym + R_yp) ) * inv_diag
        
        Z[1:-1,:] = ( a_dx2*(Z_xm + Z_xp) - b_dxdy*(Z_xpyp - Z_xmyp - Z_xpym + Z_xmym) + c_dy2*(Z_ym + Z_yp) ) * inv_diag
        

        maxchange_sq = np.amax((R-Rold)**2 + (Z-Zold)**2)
        #print(maxchange_sq)
        
        if maxchange_sq < tol:
            break

    if show:
        plt.plot(R,Z)
        plt.plot(np.transpose(R), np.transpose(Z))
        plt.show()
        
    if return_coords:
        return R, Z
    return StructuredPoloidalGrid(R,Z)

if __name__ == "__main__":
    
    import rzline
    
    #inner = circle(R0=1.5, r=1.0, n=100)
    #outer = circle(R0=1.0, r=2.0, n=100)

    inner = rzline.shaped_line(R0=3.0, a=0.5, elong=1.0, triang=0.0, indent=1.0, n=50)
    outer = rzline.shaped_line(R0=2.8, a=1.5, elong=1.0, triang=0.0, indent=0.2, n=50)
    #outer = shaped_line(R0=3.0, a=1.0, elong=1.0, triang=0.0, indent=1.0, n=50)
    
    grid = grid_elliptic(inner, outer, 100, 100, show=True)
    
    #grid.findIndex(2.0, 1.5)
    x,y = grid.findIndex([2.0,1.9], [1.5,2.0])
    print(x,y)
    
    
    
