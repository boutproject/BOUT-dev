from builtins import object

import numpy as np
from scipy.integrate import odeint

# from . import grid
# from . import field

class FieldTracer(object):
    def __init__(self, field):
        """Create a FieldTracer object

        Inputs
        ------
        field - A function specifying the magnetic field function
        """

        self.field_direction = field.field_direction

    def follow_field_lines(self, x_values, z_values, phi_values):
        """Uses field_direction to follow the magnetic field
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
        # Check len(x) == len(z)
        position = np.column_stack((x_values, z_values)).flatten()
        result = odeint(self.field_direction, position, phi_values, args=(True,))

        return result.reshape((len(phi_values), len(x_values), 2))
