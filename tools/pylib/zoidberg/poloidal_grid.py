"""
Routines for generating structured meshes on poloidal domains
"""

import numpy as np
from numpy import pi, linspace, zeros
from scipy.interpolate import RectBivariateSpline
from scipy.spatial import cKDTree as KDTree

import matplotlib.pyplot as plt

class PoloidalGrid(object):
    """
    Represents a poloidal grid in R-Z
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
        
        n = R.size
        print R.shape, n
        data = np.concatenate( (R.reshape((n,1)), Z.reshape((n,1)) ), axis=1)
        self.tree = KDTree(data)
        
    def nearestCoordinate(self, R, Z):
        """
        Finds the nearest coordinate to the given point
        
        """
        try:
            shape = R.shape
            position = np.concatenate( (R.reshape((n,1)), Z.reshape((n,1)) ), axis=1)
        except:
            # Probably just a single number
            shape = None
            position = [R,Z]
        
        # Returns distances and indices into data
        dists, inds = self.tree.query(position)
        
        # inds now contains nearest index
        print inds
        plt.plot(self.R, self.Z, '.')
        plt.plot(R, Z, 'x')
        plt.plot(self.R.ravel()[inds], self.Z.ravel()[inds], 'o')
        plt.show()

def grid_annulus(inner, outer, nx, ny, show=True):
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
    If return_coords is false, returns a PoloidalGrid object

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
    return PoloidalGrid(R,Z)

if __name__ == "__main__":
    
    import rzline
    
    #inner = circle(R0=1.5, r=1.0, n=100)
    #outer = circle(R0=1.0, r=2.0, n=100)

    inner = rzline.shaped_line(R0=3.0, a=0.5, elong=1.0, triang=0.0, indent=1.0, n=50)
    outer = rzline.shaped_line(R0=2.8, a=1.5, elong=1.0, triang=0.0, indent=0.2, n=50)
    #outer = shaped_line(R0=3.0, a=1.0, elong=1.0, triang=0.0, indent=1.0, n=50)
    
    grid = grid_elliptic(inner, outer, 100, 100, show=True)
    
    grid.nearestCoordinate(2.0, 1.5)
    
    
