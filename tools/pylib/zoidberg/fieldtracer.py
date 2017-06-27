try:
    from builtins import object
except:
    pass

import numpy as np
from scipy.integrate import odeint

class FieldTracer(object):
    def __init__(self, field):
        """Create a FieldTracer object

        Inputs
        ------
        field - A function specifying the magnetic field function
        """

        self.field_direction = field.field_direction

    def follow_field_lines(self, x_values, z_values, y_values, rtol=None):
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
        result = odeint(self.field_direction, position, y_values, args=(True,),
                        rtol=rtol)

        return result.reshape(y_values.shape + array_shape + (2,))

class FieldTracerReversible(object):
    """
    Traces magnetic field lines in a reversible way by using trapezoidal integration:
    
    pos_{n+1} = pos_n + 0.5*( f(pos_n) + f(pos_{n+1}) )*dy
    
    This requires a Newton iteration to solve the nonlinear set of equations
    for the unknown pos_{n+1}.

    """
    def __init__(self, field, rtol=1e-8, eps=1e-5, nsteps=20):
        """Create a FieldTracer object

        Inputs
        ------
        field - A function specifying the magnetic field function
        
        rtol    Tolerance applied to changes in dx**2 + dz**2
        eps     Change in x,z used to calculate finite differences of magnetic field direction
        nsteps  Number of sub-steps between outputs
        """

        self.field_direction = field.field_direction
        self.rtol = float(rtol)
        self.eps = float(eps)
        self.nsteps = int(nsteps)

    def follow_field_lines(self, x_values, z_values, y_values, rtol=None, eps=None, nsteps=None):
        """Uses field_direction to follow the magnetic field
        from every grid (x,z) point at toroidal angle y
        through a change in toroidal angle dy

        Inputs
        ------
        x_values - Array-like or scalar of starting x coordinates
        z_values - Array-like or scalar of starting z coordinates
        y_values - Array-like of y coordinates to follow the field line to.
                   y_values[0] is the starting position

        rtol    Tolerance applied to changes in dx**2 + dz**2
        eps   Change in x,z used to calculate finite differences of magnetic field direction
        nsteps  Number of sub-steps between outputs
        
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

        # Check settings, use defaults if not given
        if rtol is None:
            rtol = self.rtol # Value set in __init__
        rtol = float(rtol)
        if eps is None:
            eps = self.eps
        eps = float(eps)
        if nsteps is None:
            nsteps = self.nsteps
        nsteps = int(nsteps)

        # Ensure all inputs are NumPy arrays
        x_values = np.asfarray(x_values)
        y_values = np.asfarray(y_values)
        z_values = np.asfarray(z_values)
        
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
        
        result = np.zeros((len(y_values),)+ array_shape + (2,))

        # Starting position
        x_pos = x_values
        z_pos = z_values
        y_pos = y_values[0]
        
        result[0,...,0] = x_pos
        result[0,...,1] = z_pos
        
        for yindex, y_next in enumerate(y_values[1:]):
            yindex += 1 # Since we chopped off the first y_value
            
            # Split into sub-steps
            dy = (y_next - y_pos) / float(nsteps)
            for step in range(nsteps):
                # Evaluate gradient at current position
                dxdy, dzdy = self.field_direction( (x_pos, z_pos), y_pos )
                
                # Half-step values
                x_half = x_pos + 0.5*dxdy*dy
                z_half = z_pos + 0.5*dzdy*dy
        
                # Now need to find the roots of a nonlinear equation
                #
                # f = 0.5*(dpos/dy)(pos_next)*dy - (pos_next - pos_half) = 0
                # 
        
                # Use Euler method to get starting guess
                x_pos += dxdy*dy
                z_pos += dzdy*dy
                y_pos += dy
                
                while True:
                    dxdy, dzdy = self.field_direction( (x_pos, z_pos), y_pos )
                    
                    # Calculate derivatives (Jacobian matrix) by finite difference
                    dxdy_xe, dzdy_xe = self.field_direction( (x_pos+eps, z_pos), y_pos )
                    dxdy_x = (dxdy_xe - dxdy)/eps
                    dzdy_x = (dzdy_xe - dzdy)/eps
                    
                    dxdy_ze, dzdy_ze = self.field_direction( (x_pos, z_pos+eps), y_pos )
                    dxdy_z = (dxdy_ze - dxdy)/eps
                    dzdy_z = (dzdy_ze - dzdy)/eps
        
                    # The function we are trying to find the roots of:
                    fx = 0.5*dxdy*dy - x_pos + x_half 
                    fz = 0.5*dzdy*dy - z_pos + z_half
                    
                    # Now have a linear system to solve
                    #
                    # (x_pos)  -= ( dfx/dx   dfx/dz )^-1 (fx)
                    # (z_pos)     ( dfz/dx   dfz/dz )    (fz)
                    
                    dfxdx = 0.5*dxdy_x*dy - 1.0
                    dfxdz = 0.5*dxdy_z*dy
                    dfzdx = 0.5*dzdy_x*dy
                    dfzdz = 0.5*dzdy_z*dy - 1.0
                    
                    determinant = dfxdx*dfzdz - dfxdz*dfzdx
                    # Note: If determinant is too small then dt should be reduced
                    
                    dx = (dfzdz * fx - dfxdz * fz) / determinant
                    dz = (dfxdx * fz - dfzdx * fx) / determinant
                
                    x_pos -= dx
                    z_pos -= dz
                    # Check for convergence within tolerance
                    if np.amax(dx**2 + dz**2) < rtol:
                        break
                # Finished Newton iteration, taken step to (x_pos,y_pos,z_pos)
            # Finished sub-steps, reached y_pos = y_next
            result[yindex,...,0] = x_pos
            result[yindex,...,1] = z_pos
            
        return result
