Fluid equations 2: reduced MHD
==============================

The MHD example presented previously covered some of the functions
available in BOUT++, which can be used for a wide variety of models.
There are however several other significant functions and classes which
are commonly used, which will be illustrated using the
``reconnect-2field`` example. This is solving equations for
:math:`A_{||}` and vorticity :math:`U`

.. math::

   {{\frac{\partial U}{\partial t}}} =& -\frac{1}{B}\mathbf{b}_0\times\nabla\phi\cdot\nabla U + B^2
       \nabla_{||}(j_{||} / B) \\ {{\frac{\partial A_{||}}{\partial t}}} =&
       -\frac{1}{\hat{\beta}}\nabla_{||}\phi - \eta\frac{1}{\hat{\beta}} j_{||}

with :math:`\phi` and :math:`j_{||}` given by

.. math::

   U =& \frac{1}{B}\nabla_\perp^2\phi \\ j_{||} =& -\nabla_\perp^2 A_{||}

First create the variables which are going to be evolved, ensure
they’re communicated

::

    Field3D U, Apar; // Evolving variables

    int physics_init(bool restarting) {

      SOLVE_FOR2(U, Apar);
    }

    int physics_run(BoutReal t) {
      mesh->communicate(U, Apar);

    }

In order to calculate the time derivatives, we need the auxiliary
variables :math:`\phi` and :math:`j_{||}`. Calculating :math:`j_{||}`
from :math:`A_{||}` is a straightforward differential operation, but
getting :math:`\phi` from :math:`U` means inverting a Laplacian.

::

    Field3D U, Apar;
    Field3D phi, jpar; // Auxilliary variables

    int physics_init(bool restarting) {
      SOLVE_FOR2(U, Apar);
      SAVE_REPEAT2(phi, jpar); // Save variables in output file
      return 0;
    }

    int physics_run(BoutReal t) {
      phi = invert_laplace(mesh->Bxy*U, phi_flags); // Solve for phi
      mesh->communicate(U, Apar, phi);  // Communicate phi
      jpar = -Delp2(Apar);     // Calculate jpar
      mesh->communicate(jpar); // Communicate jpar
      return 0;
    }

Note that the Laplacian inversion code takes care of boundary regions,
so ``U`` doesn’t need to be communicated first. The differential
operator ``Delp2`` , like all differential operators, needs the values
in the guard cells and so ``Apar`` needs to be communicated before
calculating ``jpar`` . Since we will need to take derivatives of
``jpar`` later, this needs to be communicated as well.

::

    int physics_run(BoutReal t) {
      ...
      mesh->communicate(jpar);

      ddt(U) = -b0xGrad_dot_Grad(phi, U) + SQ(mesh->Bxy)*Grad_par(Jpar / mesh->Bxy)
      ddt(Apar) = -Grad_par(phi) / beta_hat - eta*jpar / beta_hat; }

.. _sec-printing:

Printing messages/warnings
--------------------------

In order to print to screen and/or a log file, the object ``output`` is
provided. This provides two different ways to write output: the C
(``printf``) way, and the C++ stream way. This is because each method
can be clearer in different circumstances, and people have different
tastes in these matters.

The C-like way (which is the dominant way in BOUT++) is to use the
``write`` function, which works just like ``printf``, and takes all the
same codes (it uses ``sprintf`` internally).

::

    output.write(const char *format, ...)

For example:

::

    output.write("This is an integer: %d, and this a real: %e\n", 5, 2.0)

For those who prefer the C++ way of doing things, a completely
equivalent way is to treat ``output`` as you would ``cout``:

::

    output << "This is an integer: " << 5 << ", and this a real: " << 2.0 << endl;

which will produce exactly the same result as the ``output.write`` call
above.

On all processors, anything sent to ``output`` will be written to a log
file called ``BOUT.log.#`` with # replaced by the processor number. On
processor 0, anything written to the output will be written to screen
(stdout), in addition to the log file. Unless there is a really good
reason not to, please use this ``output`` object when writing text
output.

More details are given in section :ref:`sec-logging`.

Error handling
--------------

Finding where bugs have occurred in a (fairly large) parallel code is a
difficult problem. This is more of a concern for developers of BOUT++
(see the developers manual), but it is still useful for the user to be
able to hunt down bug in their own code, or help narrow down where a bug
could be occurring.

If you have a bug which is easily reproduceable i.e. it occurs almost
immediately every time you run the code, then the easiest way to hunt
down the bug is to insert lots of ``output.write`` statements (see
:ref:`sec-printing`). Things get harder when a bug only occurs after
a long time of running, and/or only occasionally. For this type of
problem, a useful tool can be the message stack. An easy way to use this message
stack is to use the ``TRACE`` macro:

::

	{
      	  TRACE("Some message here"); // message pushed
	
	} // Scope ends, message popped

This will push the message, then pop the message when the current scope ends
(except when an exception occurs).
The error message will also have the file name and line number appended, to help find
where an error occurred. The run-time overhead of this should be small,
but can be removed entirely if the compile-time flag ``-DCHECK`` is not defined or set to ``0``. This turns off checking,
and ``TRACE`` becomes an empty macro.
It is possible to use standard ``printf`` like formatting with the trace macro, for example.
 
::

	{
      	  TRACE("The value of i is %d and this is an arbitrary %s", i, "string"); // message pushed
	} // Scope ends, message popped


