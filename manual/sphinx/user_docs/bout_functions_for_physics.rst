BOUT++ functions (alphabetical)
===============================

This is a list of functions which can be called by users writing a
physics module. For a full list of functions, see the Reference manual,
DOxygen documentation, and source code.

-  ``Field = abs(Field | Vector)``

-  | ``(Communicator).add(Field | Vector)``
   | Add a variable to a communicator object.

-  ``apply_boundary(Field. name)``

-  ``Field = b0xGrad_dot_Grad(Field, Field, CELL_LOC)``

-  ``bout_solve(Field, Field, name)``

-  ``bout_solve(Vector, Vector, name)``

-  | ``(Communicator).clear()``
   | Remove all variables from a Communicator object

-  ``Field = cos(Field)``

-  ``Field = cosh(Field)``

-  ``Vector = Curl(Vector)``

-  | ``Field = Delp2(Field)``
   | :math:`\nabla_\perp^2` operator

-  | ``Field = Div(Vector)``
   | Divergence of a vector

-  | ``Field = Div_par(Field f)``
   | Parallel divergence :math:`B_0\mathbf{b}\cdot\nabla(f / B_0)`

-  ``dump.add(Field, name, 1/0)``

-  ``Field = filter(Field, modenr)``

-  | ``geometry_derivs()``
   | Calculates useful quantities from the metric tensor. Call this
     every time the metric tensor is changed.

-  ``Vector = Grad(Field)``

-  ``Field = Grad_par(Field)``

-  ``Field = Grad2_par2(Field)``

-  | ``grid_load(BoutReal, name)``
   | Load a scalar real from the grid file

-  | ``grid_load2d(Field2D, name)``
   | Load a 2D scalar field from the grid file

-  | ``grid_load3d(Field3D, name)``
   | Load a 3D scalar field from the grid file

-  ``invert_laplace(Field input, Field output, flags, Field2D *A)``

-  | ``Field = invert_parderiv(Field2D|BoutReal A, Field2D|BoutReal B, Field3D r)``
   | Inverts an equation ``A*x + B*Grad2_par2(x) = r``

-  ``Field = Laplacian(Field)``

-  ``Field3D = low_pass(Field3D, max_modenr)``

-  ``BoutReal = max(Field)``

-  ``BoutReal = min(Field)``

-  | ``msg_stack.pop( |int)``
   | Remove a message from the top of the stack. If a message ID is
     passed, removes all messages back to that point.

-  | ``int = msg_stack.push(format, ...)``
   | Put a message onto the stack. Works like ``printf`` (and
     ``output.write``).

-  | ``options.get(name, variable, default)``
   | Get an integer, real or boolean value from the options file. If not
     in the file, the default value is used. The value used is printed
     to log file.

-  ``options.setSection(name)`` Set the section name in the input file

-  | ``output < < values``
   | Behaves like cout for stream output

-  | ``output.write(format, ...)``
   | Behaves like printf for formatted output

-  | ``(Communicator).receive()``
   | Receive data from other processors. Must be preceded by a ``send``
     call.

-  | ``(Communicator).run()``
   | Sends and receives data.

-  | ``(Communicator).send()``
   | Sends data to other processors (and posts receives). This must be
     followed by a call to ``receive()`` before calling send again, or
     adding new variables.

-  ``(Field3D).setLocation(CELL_LOC)``

-  ``(Field3D).ShiftZ(bool)``

-  ``Field = sin(Field)``

-  ``Field = sinh(Field)``

-  | ``solver.setPrecon(PhysicsPrecon)``
   | Set a preconditioner function

-  ``Field = sqrt(Field)``

-  ``Field = tan(Field)``

-  ``Field = tanh(Field)``

-  | ``Field = V_dot_Grad(Vector v, Field f)``
   | Calculates an advection term :math:`\mathbf{v}\cdot\nabla f`

-  | ``Vector = V_dot_Grad(Vector v, Vector u)``
   | Advection term :math:`\mathbf{v}\cdot\nabla\mathbf{u}`

-  ``Field = Vpar_Grad_par(Field v, Field f)``

-  | ``Field3D = where(Field2D test, Field|BoutReal gt0, Field|BoutReal lt0)``
   | Chooses between two values, depending on sign of ``test``.

