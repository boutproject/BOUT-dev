from __future__ import division
#from builtins import object

import numpy as np

# Import classes representing poloidal grids
from .poloidal_grid import RectangularPoloidalGrid, StructuredPoloidalGrid

class OldGrid(object):
    def __init__(self, nx, ny, nz,
                 Lx=0.1, Ly=10., Lz = 1.,
                 name="fci_grid", MXG=2):
        """Initialise a grid object

        Inputs
        ------
        nx  - Number of radial points
        ny  - Number of toroidal points (NOTE: Different to BOUT++ standard)
        nz  - Number of poloidal points

        Lx  - Radial domain size  [m]
        Ly  - Toroidal domain size [m]
        Lz  - Poloidal domain size [m]

        name - Grid name (default "fci_grid")

        MXG - Number of x guard cells
        """

        self.MXG = MXG

        self.name = name

        self.nx = int(nx)
        self.ny = int(ny)
        self.nz = int(nz)

        self.Lx = float(Lx)
        self.Ly = float(Ly)
        self.Lz = float(Lz)

        self.delta_x = self.Lx/float(nx-2.*MXG)
        self.delta_y = self.Ly/float(ny)
        self.delta_z = self.Lz/float(nz)

        # Coord arrays
        self.xarray = Lx * (np.arange(nx) - MXG + 0.5)/(nx - 2.*MXG)  # 0 and 1 half-way between cells
        self.yarray = np.linspace(0,Ly,ny,endpoint=False)
        self.zarray = np.linspace(0,Lz,nz,endpoint=False)

        self.xcentre = Lx/2.#0.5*max(self.xarray)
        self.zcentre = 0.5*max(self.zarray)

        self.x_3d, self.y_3d, self.z_3d = np.meshgrid(self.xarray, self.yarray, self.zarray,
                                                      indexing='ij')

        # How to do this properly?
        def Rmaj(x,z,phi):
            return 1.0

        self.Rmaj = Rmaj

    def __repr__(self):
        return "Grid(nx={nx}, ny={ny}, nz={nz}, Lx={Lx}, Ly={Ly}, Lz={Lz}, name='{name}')".format(
            nx=self.nx, ny=self.ny, nz=self.nz, Lx=self.Lx, Ly=self.Ly, Lz=self.Lz, name=self.name)


class Grid(object):
    """
    Represents a 3D grid, consisting of a collection of poloidal grids

    Members
    -------
    
    shape     tuple of grid sizes (nx,ny,nz)
    
    """
    def __init__(self, poloidal_grids, ycoords, Ly, yperiodic=False, name="fci_grid"):
        """
        Inputs
        ------
        
        poloidal_grids
        ycoords

        Examples
        --------

        grid = Grid()
        
        To iterate over the poloidal grids, and get the grids to either side:
        
        for i in range(grid.numberOfPoloidalGrids()):
            pol, y = grid.getPoloidalGrid(i)
            pol_next, y_next = grid.getPoloidalGrid(i+1)
            pol_last, y_last = grid.getPoloidalGrid(i-1)
        
        The getPoloidalGrid() method ensures that y_last <= y <= y_next
        """
        
        try:
            ngrids = len(poloidal_grids)
            
            # Check this is the same length as ycoords
            assert len(ycoords) == ngrids

            nx = poloidal_grids[0].nx
            nz = poloidal_grids[0].nz
        except TypeError:
            # No len(), assume single poloidal grid
            ngrids = 1
            nx = poloidal_grids.nx
            nz = poloidal_grids.nz

        self.poloidal_grids = poloidal_grids
        self._ngrids = ngrids  # This is an implementation detail, whether we have one or multiple separate grids
        self.ycoords = np.asfarray(ycoords)
        self.Ly = Ly
        self.yperiodic = yperiodic
        
        # Define the shape of the grid
        self.shape = (nx, len(ycoords), nz)
        
    def __repr__(self):
        return "Grid({0}, {1}:{2}:{3}, yperiodic={4})".format(self.poloidal_grids, 
                                                              np.amin(self.ycoords), np.amax(self.ycoords), self.ycoords.size, 
                                                              self.yperiodic)
    def numberOfPoloidalGrids(self):
        """
        Returns the number of poloidal grids i.e. number of points in Y
        """
        return self.ycoords.size

    def getPoloidalGrid(self, yindex):
        """
        Returns the poloidal grid and y value at the given y index
        
        This handles negative values and values out of range, if the domain is periodic
        
        """
        yindex = int(yindex)
        
        ny = self.ycoords.size

        if (yindex >= 0) and (yindex < ny):
            # Within index range, so just return
            if self._ngrids == 1:
                # Only one grid
                return self.poloidal_grids, self.ycoords[yindex]
            return self.poloidal_grids[yindex], self.ycoords[yindex]
                
        # Out of range
        
        if self.yperiodic:
            # Periodic domain

            # Map index into domain
            y_remap = np.remainder(yindex, ny) # 0 <= yremap < ny

            # Get number of periods around the domain. Note this can be negative
            nperiods = np.floor( float(yindex) / float(ny) )

            ycoord = self.ycoords[y_remap] + nperiods*self.Ly
            
            if self._ngrids == 1:
                return self.poloidal_grids, ycoord
            return self.poloidal_grids[y_remap], ycoord
                        
        # Not periodic
        
        if yindex < 0:
            return None, 0.0   # Hit the lower end in Y
        return None, self.Ly  # Hit the upper end in Y
        
    def metric(self):
        """
        Return the metric tensor, dx and dz
        """
        
        # Gather dx,dz and x-z metrics from poloidal slices
        dx = np.zeros(self.shape)
        dz = np.zeros(self.shape)
        
        gxx = np.zeros(self.shape)
        gxz = np.zeros(self.shape)
        gzz = np.zeros(self.shape)
        
        g_xx = np.zeros(self.shape)
        g_xz = np.zeros(self.shape)
        g_zz = np.zeros(self.shape)
        
        # Separate grids for each slice
        for y in range(self.shape[1]):
            pol_metric = self.getPoloidalGrid(y)[0].metric()
            dx[:,y,:] = pol_metric["dx"]
            dz[:,y,:] = pol_metric["dz"]
                
            gxx[:,y,:] = pol_metric["gxx"]
            gxz[:,y,:] = pol_metric["gxz"]
            gzz[:,y,:] = pol_metric["gzz"]
            
            g_xx[:,y,:] = pol_metric["g_xx"]
            g_xz[:,y,:] = pol_metric["g_xz"]
            g_zz[:,y,:] = pol_metric["g_zz"]
        
        # Calculate the gradient of the y coordinate w.r.t index
        # To avoid edge effects, repeat array three times then take the middle
        ycoords = np.concatenate( (self.ycoords - self.Ly, self.ycoords, self.ycoords + self.Ly) )
        ny = self.ycoords.size
        dy = np.gradient(ycoords[ny:(2*ny)])
        
        dy3d = np.zeros(self.shape)
        for i in range(self.shape[1]):
            dy3d[:,i,:] = dy[i]

        # Note: These y metrics are for Cartesian coordinates
        # If in cylindrical coordinates then these should be different
        g_yy = 1.0   # Rmaj**2
        gyy = 1.0    # 1/Rmaj**2

        return {"dx":dx, "dy":dy3d, "dz": dz,
                "gyy": gyy,  "g_yy":g_yy,
                "gxx": gxx,  "g_xx":g_xx,
                "gxz": gxz,  "g_xz":g_xz,
                "gzz": gzz,  "g_zz":g_zz}

def rectangular_grid(nx, ny, nz,
                     Lx=1.0, Ly=10., Lz=1.0,
                     xcentre=0.0, zcentre=0.0,
                     yperiodic=False):
    """
    Create a rectangular grid in (x,y,z). Here y is along the
    magnetic field (typically toroidal angle), and (x,z) are in
    the poloidal plane
    
    nx,ny,nz   Number of points in x,y,z
    Lx,Ly,Lz   Size of the domain in x,y,z
    xcentre, zcentre   The middle of the domain
    yperiodic   Determines if the y direction is periodic
    
    
    """
    
    # In this simple case we only need one RectangularPoloidalGrid
    # to represent all poloidal planes
    
    poloidal_grid = RectangularPoloidalGrid(nx, nz, Lx, Lz, 
                                            xcentre,
                                            zcentre)

    if yperiodic:
        ycoords = np.linspace(0.0, Ly, ny, endpoint=False)
    else:
        # Doesn't include the end points
        ycoords = (np.arange(ny) + 0.5)*Ly/float(ny)
    
    return Grid( poloidal_grid, 
                 ycoords, 
                 Ly,
                 yperiodic=yperiodic )

    
if __name__ == "__main__":
    
    grid = rectangular_grid(10,10,10)

    p = grid.getPoloidalGrid(-2)
    
    print(grid)
