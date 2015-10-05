from builtins import object

import numpy as np
from scipy.integrate import odeint

# from . import grid
# from . import field

class FieldTracer(object):
    def __init__(self, grid, field):
        """Create a FieldTracer object

        Inputs
        ------
        grid - An FCI grid
        field - A function specifying the magnetic field function
        """

        self.grid = grid
        self.field_direction = field.field_direction

    def follow_field_line(self, phi, dphi):
        """
        Uses field_direction to follow the magnetic field
        from every grid (x,z) point at toroidal angle phi
        through a change in toroidal angle dphi

        Inputs
        ------
        phi  = Starting toroidal angle [radians]
        dphi = Change in toroidal angle [radians]

        Returns
        -------

        Field line ending coordinates

        result[x,z,0] = x end location [m] from start index (x,z)
        result[x,z,1] = z end location [m] from start index (x,z)

        """
        result = np.zeros( (self.grid.nx, self.grid.nz, 2) )

        # TODO: flatten loop
        for i in np.arange(0,self.grid.nx):
            for k in np.arange(0,self.grid.nz):
                result[i,k,:] = odeint(self.field_direction,         # Function to integrate
                                       [self.grid.xarray[i], self.grid.zarray[k]],  # x [m], z[m]
                                       [phi, phi+dphi])[1,:]

        return result
