.. _sec-newapi:

Object-orientated interface
===========================

| If you prefer to create classes rather than global variables and C
  functions for your physics model, this can be done using a (somewhat
  experimental) interface. To see the difference, compare
  ``examples/advect1d/gas_compress.cxx`` with
| ``examples/advect1d-newapi/gas_compress.cxx``. The disadvantage of
  this interface is that it’s marginally more complicated to set up, but
  it has several advantages: It makes splitting the model into multiple
  files easier (sharing global variables is a pain), models can be
  combined together to enable coupling of models, and BOUT++ can be more
  easily used alongside other libraries. For large models, it’s
  recommended to use this method. Converting C-style interface to a
  class is also quite straightforward, and discussed below.

In a header file (e.g. ``examples/advect1d-newapi/gas_compress.hxx``),
first put

::

    #include <bout/physicsmodel.hxx>

(do NOT include ``boutmain.hxx``, as that defines the C-like interface
and a ``main()`` function).

Next define a class which inherits from ``PhysicsModel``

::

    class GasCompress : public PhysicsModel {
    protected:
      int init(bool restarting);
      int rhs(BoutReal t);
    private:
      // Evolving variables, parameters etc. here
    };

As a minimum, you need to define the initialisation function ``init``
(it’s a pure virtual member of PhysicsModel, so if you don’t you’ll get
a compile-time error). Any variables being evolved should now be members
of this class. If you are converting a C-style model, just move all the
global variables into the ``private`` section.

Next create a source file (e.g.
``examples/advect1d-newapi/gas_compress.cxx``, which includes your
header file

::

    #include "gas_compress.hxx"

Then implement the init and rhs functions:

::

    int GasCompress::init(bool restarting) {
      ...
    }

    int GasCompress::rhs(BoutReal t) {
      ...
    }

To convert simple physics models, just rename ``physics_init`` to
``YourModel::init`` , and ``physics_run`` to ``YourModel::run`` .

Finally, you need to create a ``main()`` function for your code. The
easiest way to do this is to use the macro ``BOUTMAIN`` :

::

    BOUTMAIN(GasCompress);

This is defined in ``include/bout/physicsmodel.hxx``, and expands to

::

      int main(int argc, char **argv) {
        BoutInitialise(argc, argv); // Initialise BOUT++

        GasCompress *model = new GasCompress(); // Create a model

        Solver *solver = Solver::create(); // Create a solver
        solver->setModel(model); // Specify the model to solve
        solver->addMonitor(bout_monitor); // Monitor the solver

        solver->solve(); // Run the solver

        delete model;
        delete solver;
        BoutFinalise(); // Finished with BOUT++
        return 0;
      }

If you like, you can define your own ``main()`` function, making it
easier to combine BOUT++ with other libraries.

