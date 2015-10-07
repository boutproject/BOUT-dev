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

    def follow_all_field_lines(self, phi, dphi):
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

        x2d, z2d = np.meshgrid(self.grid.xarray, self.grid.zarray, indexing='ij')
        result = self.follow_single_field_line(x2d.flatten(), z2d.flatten(), [phi, phi+dphi])

        return result.reshape( (2, self.grid.nx, self.grid.nz, 2) )[1,...]

    def follow_single_field_line(self, x_values, z_values, phi_values):

        position = np.column_stack((x_values, z_values)).flatten()
        result = odeint(self.field_direction, position, phi_values, args=(True,))

        return result.reshape((-1, 2))
