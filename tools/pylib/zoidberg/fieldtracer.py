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

    def follow_field_lines(self, x_values, z_values, y_values):
        """Uses field_direction to follow the magnetic field
        from every grid (x,z) point at toroidal angle y
        through a change in toroidal angle dy

        Inputs
        ------
        x_values - Array-like or scalar of starting x coordinates
        z_values - Array-like or scalar of starting z coordinates
        y_values - Array-like of y coordinates to follow the field line to.
                   y_values[0] is the starting position

        Returns
        -------

        Field line ending coordinates

        result - The first dimension is y, the last is (x,z). The
                 middle dimensions are the same shape as [x|z]:
                 [0,...] is the initial position
                 [...,0] are the x-values
                 [...,1] are the z-values
                 If x_values is a scalar and z_values a 1D array, then result
                 has the shape [len(y), len(z), 2], and vice-versa.
                 If x_values and z_values are 1D arrays, then result has the shape
                 [len(y), len(x), 2].
                 If x_values and z_values are 2D arrays, then result has the shape
                 [len(y), x.shape[0], x.shape[1], 2].

        """

        if not isinstance(x_values, np.ndarray):
            x_values = np.array(x_values)
        if not isinstance(y_values, np.ndarray):
            y_values = np.array(y_values)
        if not isinstance(z_values, np.ndarray):
            z_values = np.array(z_values)
        if len(y_values) < 2:
            raise ValueError("There must be at least two elements in y_values")
        if len(y_values.shape) > 1:
            raise ValueError("y_values must be 1D")

        if x_values.shape != z_values.shape:
            # Make the scalar the same shape as the array
            if x_values.size is 1:
                x_values = np.zeros(z_values.shape) + x_values
            elif z_values.size is 1:
                z_values = np.zeros(x_values.shape) + z_values
            else:
                raise ValueError("x_values and z_values must be the same size, or one must be a scalar")

        array_shape = x_values.shape

        # Position vector must be 1D - so flatten before passing to
        # integrator, then reshape after
        if len(x_values.shape) > 1:
            x_values = x_values.flatten()
            z_values = z_values.flatten()
        position = np.column_stack((x_values, z_values)).flatten()
        result = odeint(self.field_direction, position, y_values, args=(True,))

        return result.reshape(y_values.shape + array_shape + (2,))
