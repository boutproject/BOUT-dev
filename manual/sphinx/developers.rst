Introduction
============

.. toctree::
   :maxdepth: 1
   :caption: Contents:

This is a manual describing the core BOUT++ [1]_, and is intended for
anyone who wants to work on improving BOUT++. It does its best to
describe the details of how BOUT++ works, and assumes that the user is
very comfortable with C++. For a general introduction, and
instructions for using BOUT++ see the users’ guide. The user’s guide
assumes only minimal knowledge of C++, and provides only those details
needed to use BOUT++.

Since BOUT++ is a scientific code, it is constantly changing and
(hopefully) being improved. This provides a moving target for
documentation, and means that describing all the details of how BOUT++
works in one document is probably impossible. This is particularly
true since often our first priority is to write papers and code - not
documentation - and so whatever is documented is likely to be slightly
out of date. A source of up-to-date documentation of the BOUT++ code
is the comments and Doxygen tags: running ``doxygen`` on the source
should produce a set of HTML documentation. See www.doxygen.org for
more details.

.. [1] http://www.sciencedirect.com/science/article/B6TJ5-4VTCM95-3/2/ed200cd23916d02f86fda4ce6887d798

Using the BOUT++ repository
---------------------------

The BOUT++ distribution is hosted on Github:

https://github.com/boutproject/BOUT-dev

For a full guide to using Git, see the `git website`_ or online
tutorials. This manual just explains some basic ways to use Git, and the
recommended work flow when working with BOUT++.

.. _git website: http://git-scm.com/

If you’re just starting with BOUT++, current developers will want to
check your changes before submitting them to the repository. In this
case you should fork the git repository, make any changes and then
submit a pull request on Github. Fortunately Git makes this process
quite easy: First get a copy of BOUT++

.. code-block:: bash

    $ git clone https://github.com/boutproject/BOUT-dev.git

The BOUT++ repository will now be in a directory called “BOUT-dev”
(sorry - github doesn’t like ’+’ in project names). To get the latest
changes, use

.. code-block:: bash

    $ git pull

To see the status of the repository, commits etc. in a GUI use

.. code-block:: bash

    $ gitk

This is also useful for showing what changes you’ve made which need to
be committed, or which haven’t yet been sent to the main repository.

You can make edits as normal, and commit them using

.. code-block:: bash

    $ git commit -a

which is pretty much the equivalent of ``svn commit`` in that it commits
all changes, though importantly it doesn’t send them to a central
server. To see which changes will be committed, use

.. code-block:: bash

    $ git status

To choose which files you want to commit, use

.. code-block:: bash

    $ git add file1, file2, ...
    $ git commit

(Git can actually only commit selected parts of files if you want). To
make using Git easier, you can create a config file ``$HOME/.gitconfig``
containing:

.. code-block:: cfg

    [user]
        name = A. Developer
        email = a.developer@example.com

    [alias]
        st = status
        ci = commit
        br = branch
        co = checkout
        df = diff
        lg = log -p
        who = shortlog -s --

(though obviously you should change the name and email).

Once you’re done making changes, you should first pull the latest
changes from the server:

.. code-block:: bash

    $ git pull

**Read carefully** what git prints out. If there are conflicts then git
will try to resolve them, but in some cases you will have to resolve
them yourself. To see a list of conflicting changes run ``git status``
(or ``git st`` if you’re using the above ``.gitconfig`` file). Once
you’ve finished resolving conflicts, run ``git commit -a`` to commit the
merge.


Developing BOUT++
~~~~~~~~~~~~~~~~~

If you are doing a lot of development of BOUT++, it will probably make
sense for you to push changes directly to the online repository. In this
case you’ll need to sign up for an account on ``github.com``, then
upload an ssh key and ask to be added. The process is then almost
identical except that you clone using SSH:

.. code-block:: bash

    $ git clone https://github.com/boutproject/BOUT-dev.git

and rather than creating a patch, you push changes to the repository:

.. code-block:: bash

    $ git push

Accessing github from behind a firewall
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you’re working on a machine which can’t access github directly (such
as grendel, smaug etc. at LLNL), you can still seamlessly access github
by using another machine as a proxy over SSH. To do this, edit your SSH
config file `` /.ssh/config`` and add the following lines:

.. code-block:: aconf

    Host            gh
    HostName        github.com
    User            git
    ProxyCommand    ssh -q -x user@euclid.nersc.gov nc %h %p

where ``euclid.nersc.gov`` can be replaced by any machine you can access
which has netcat (``nc``) installed, and which can access github.com. If
you have set up a github account with SSH keys, you should now be able
to get a copy of BOUT++ by running

.. code-block:: bash

    $ git clone gh:boutproject/BOUT-dev.git

Creating a private repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Whilst we would prefer it if improvements to BOUT++ were shared,
sometimes you might want to keep changes private for a while before
publishing them. Creating a private repository with Git is very simple,
because every clone of a repository is itself a repository. Git doesn’t
have the concept of a central repository, which can seem strange coming
from the world of SVN and CVS. What it means is that you can create your
own private repository anywhere you have access to. Sharing it with only
some people means as giving them read or write access to the repository
directory.

The following assumes you have a NERSC account and want to create a
private repository on Franklin. To apply this to a different machine
just replace ``franklin.nersc.gov`` with the machine you want to put the
repository on.

#. SSH to ``franklin.nersc.gov``, or wherever you want your repository

   .. code-block:: bash

           $ ssh username@franklin.nersc.gov
         

#. Create a “bare” Git repository by cloning a repository with the
   ``–bare`` option:

   .. code-block:: bash

           $ cd ~
           $ git clone --bare git@github.com:boutproject/BOUT-dev.git  bout_private
         

   where you can replace ``git@github.com:boutproject/BOUT-dev.git`` with any
   other repository you can access. ``bout_private`` will be the name of
   the directory which will be created. This will make a repository
   without a working version. This means you can’t modify the code in it
   directly, but can pull and push changes to it. If you want to work on
   the code on Franklin, make a clone of your private repository:

   .. code-block:: bash

           $ git clone bout_private bout
         

   which creates a repository ``bout`` from your private repository.
   Running ``git pull`` and ``git push`` from within this new repository
   will exchange patches with your ``bout_private`` repository.

#. You can now clone, pull and push changes to your private repository
   over SSH e.g.

   .. code-block:: bash

           $ git clone username@franklin.nersc.gov:bout_private
         

#. To keep your private repository up to date you may want to pull
   changes from github into your private repository. To do this, you
   need to use a third repository. Log into Franklin again:

   .. code-block:: bash

           $ cd ~
           $ git clone bout_private bout_tmp
         

   This creates a repository ``bout_tmp`` from your private repository.
   Now cd to the new directory and pull the latest changes from github:

   .. code-block:: bash

           $ cd bout_tmp
           $ git pull git://github.com/boutproject/BOUT-dev.git
         

   Note: You should be able to access this repository from Franklin, but
   if not then see the previous subsection for how to access github from
   behind a firewall.

#. This pull might result in some conflicts which need to be resolved.
   If so, git will tell you, and running

   .. code-block:: bash

           $ git status
         

   will give a list of files which need to be resolved. Edit each of the
   files listed, and when you’re happy commit the changes

   .. code-block:: bash

           $ git commit -a
         

#. Your ``bout_tmp`` directory now contains a merge of your private
   repository and the repository on github. To update your private
   repository, just push the changes back:

   .. code-block:: bash

           $ git push
         

   You can now delete the ``bout_tmp`` repository if you want.

House rules
-----------

BOUT++ consists of about 60,000 lines of C/C++  [2]_ along with 18,500
lines of IDL and 12,000 of Python. Of this, about 40,000 lines is the
core BOUT++ code, and the remainder a mix of pre- and post-processors,
and physics modules. As production codes go, this is not particularly
huge, but it is definitely large enough that keeping the code ‘clean’
and understandable is necessary. This is vital if many people are going
to work on the code, and also greatly helps code debugging and
verification. There are therefore a few house rules to keep in mind when
modifying the BOUT++ code.

When modifying the core BOUT++ code, please keep in mind that this
portion of the code is intended to be general (i.e. independent of any
particular physical system of equations), and to be used by a wide range
of users. Making code clear is also more important in this section than
the physics model since the number of developers is potentially much
greater.

Here are some rules for editing the core BOUT++ code:

- **NO FORTRAN**. EVER. Though it may be tempting for scientific
   programmers to use a little Fortran now and then, please please
   don’t put any into BOUT++.  Use of Fortran, particularly when mixed
   with C/C++, is the cause of many problems in porting and modifying
   codes.

-  If a feature is needed to study a particular system, only include it
   in the core code if it is more generally applicable, or cannot be put
   into the physics module.

Coding conventions
------------------

See CONTRIBUTING.md for guidelines on naming and other coding
conventions used within BOUT++.

Code layout
===========

BOUT++ is organised into classes and groups of functions which operate
on them: It’s not purely object-oriented, but takes advantage of many of
C++’s object-oriented features.

Figure [fig:layout1] shows the most important parts of BOUT++ and how
they fit together.

.. figure:: figs/layout1.pdf
   :alt: Overview of BOUT++ control flow

   Overview of BOUT++ control flow during initialisation (red), and
   running (blue)

The initialisation process is shown in red: basic information is first
read from the grid file (e.g. size of the grid, topology etc.), then the
user-supplied initialisation code is called. This code can read other
variables from the grid, and makes at least one call to ``bout_solve``
to specify a variable to be evolved. The main thing ``bout_solve`` does
is to add these variables to the solver.

The process of running a timestep is shown in blue in
figure [fig:layout1]: The main loop calls the solver, which in turn
calls PVODE. To evolve the system PVODE makes calls to the RHS function
inside solver. This moves data between PVODE and BOUT++, and calls the
user-supplied ``physics_run`` code to calculate time-derivatives. Much
of the work calculating time-derivatives involves differential
operators.

Calculation of the RHS function ``physics_run``, and handling of data in
BOUT++ involves many different components. Figure [fig:layout2] shows
(most) of the classes and functions involved, and the relationships
between them. Some thought was put into how this should be organised,
but it has also changed over time, so some parts could be cleaner.

.. figure:: figs/layout2.pdf
   :alt: Relationships used in calculating the RHS function

   Relationship between important classes and functions used in
   calculating the RHS function

Directories
-----------

The source code for the core of BOUT++ is divided into include files
(which can be used in physics models) in ``bout++/include``, and source
code and low-level includes in ``bout++/src``. Many parts of the code
are defined by their interface, and can have multiple different
implementations. An example is the time-integration solvers: many
different implementations are available, some of which use external
libraries, but all have the same interface and can be used
interchangeably. This is reflected in the directory structure inside
``bout++/src``. A common pattern is to store individual implementations
of an interface in a subdirectory called ``impls``.

::

    include/foo.hxx
    src/.../foo.cxx
    src/.../foo_factory.hxx
    src/.../foo_factory.cxx
    src/.../impls/one/one.hxx
    src/.../impls/one/one.cxx

where ``foo.hxx`` defines the interface, ``foo.cxx`` implements common
functions used in several implementations. ``foo_factory`` creates new
implementations, and is the only file which includes all the
implementations. Individual implementations are stored in their own
subdirectories of ``impls``. Components which follow this pattern
include ``fileio`` formats, ``invert/laplace`` and ``invert/parderiv``
inversion codes, ``mesh``, and ``solver``.

The current source code files are:

- :doc:`bout++.cxx<_breathe_autogen/file/bout_09_09_8cxx>`: Main file
  which initialises, runs and finalises BOUT++. Currently contains a
  :cpp:func:`main()` function, though this is being removed shortly.

- field

   - :doc:`field2d.cxx<_breathe_autogen/file/field2d_8cxx>` implements
     the :cpp:class:`Field2D` class. This is a scalar field which
     varies only in :math:`x` and :math:`y` and is used for things
     like metric tensor components and initial profiles. It supplies
     lots of overloaded operators and functions on these objects.

   - :doc:`field3d.cxx<_breathe_autogen/file/field3d_8cxx>` implements
     the :cpp:class:`Field3D` class, which varies in :math:`x`,
     :math:`y` and :math:`z`. Since these handle a lot more memory
     than Field2D objects, the memory management is more complicated
     and includes reference counting. See section [sec:memorymanage]
     for more details.

   - :doc:`field_data.cxx<_breathe_autogen/file/field__data_8cxx>`
     Implements some functions in the :cpp:class:`FieldData`
     class. This is a mainly pure virtual interface class which is
     inherited by :cpp:class:`Field2D` and :cpp:class:`Field3D`.

   - :doc:`fieldperp.cxx<_breathe_autogen/file/fieldperp_8cxx>`
     implements a :cpp:class:`FieldPerp` class to store slices
     perpendicular to the magnetic field i.e. they are a function of
     :math:`x` and :math:`z` only. This is mainly used for Laplacian
     inversion routines, and needs to be integrated with the other
     fields better.

   - :doc:`initialprofiles.cxx<_breathe_autogen/file/initialprofiles_8cxx>`
     routines to set the initial values of fields when a simulation
     first starts. Reads settings from the option file based on the name
     of the variable.

   - :doc:`vecops.cxx<_breathe_autogen/file/vecops_8cxx>` a collection
     of function to operate on vectors.  Contains things like
     ``Grad``, ``Div`` and ``Curl``, and uses a combination of field
     differential operators (in
     :doc:`difops.cxx<_breathe_autogen/file/difops_8cxx>`) and metric
     tensor components (in :cpp:class:`Mesh`).

   - :doc:`vector2d.cxx<_breathe_autogen/file/vector2d_8cxx>`
     implements the :cpp:class:`Vector2D` class, which uses a
     :cpp:class:`Field2D` object for each of its 3
     components. Overloads operators to supply things like dot and
     cross products.

   - :doc:`vector3d.cxx<_breathe_autogen/file/vector3d_8cxx>`
     implements :cpp:class:`Vector3D` by using a :cpp:class:`Field3D`
     object for each component.

   - :doc:`where.cxx<_breathe_autogen/file/where_8cxx>` supplies
     functions for choosing between values based on selection
     criteria.

- fileio

   - :doc:`datafile.cxx<_breathe_autogen/file/datafile_8cxx>`
     supplies an abstract :cpp:class:`DataFile` interface for data
     input and output. Handles the conversion of data in fields and
     vectors into blocks of data which are then sent to a specific
     file format.

   - :doc:`formatfactory.cxx<_breathe_autogen/file/formatfactory_8cxx>`

   - :doc:`formatfactory.hxx<_breathe_autogen/file/formatfactory_8hxx>`

   - impls

      - :doc:`emptyformat.hxx<_breathe_autogen/file/emptyformat_8hxx>`

      - hdf5

         - :doc:`h5_format.cxx<_breathe_autogen/file/h5__format_8cxx>` implements an
           interface to the HDF5 library

         - :doc:`h5_format.hxx<_breathe_autogen/file/h5__format_8hxx>`

      - netcdf

         - :doc:`nc_format.cxx<_breathe_autogen/file/nc__format_8cxx>` implements an
           interface to the NetCDF-4 library

         - :doc:`nc_format.hxx<_breathe_autogen/file/nc__format_8hxx>`

      - netcdf4

         - :doc:`ncxx<_breathe_autogen/file/ncxx4_8cxx>`
           implements an interface to the NetCDF-4 library using the
           C++ API

         - :doc:`ncxx<_breathe_autogen/file/ncxx4_8hxx>`

      - pnetcdf

         - :doc:`pnetcdf.cxx<_breathe_autogen/file/pnetcdf_8cxx>`
           Parallel NetCDF interface

         - :doc:`pnetcdf.hxx<_breathe_autogen/file/pnetcdf_8hxx>`

- invert

   - :doc:`fft_fftw.cxx<_breathe_autogen/file/fft__fftw_8cxx>`
     implements the :doc:`fft.hxx<_breathe_autogen/file/fft_8hxx>`
     interface by calling the Fastest Fourier Transform in the West
     (FFTW) library.

   - :doc:`full_gmres.cxx<_breathe_autogen/file/full__gmres_8cxx>`

   - :doc:`inverter.cxx<_breathe_autogen/file/inverter_8cxx>` is a
     :cpp:class:`FieldPerp` inversion class currently under
     development. It is intended to provide a way to solve nonlinear
     problems using a GMRES iterative method.

   - :doc:`invert_gmres.cxx<_breathe_autogen/file/invert__gmres_8cxx>`

   - :doc:`invert_laplace_gmres.cxx<_breathe_autogen/file/invert__laplace__gmres_8cxx>` inherits
     the :cpp:class:`Inverter` class and will solve more general
     Laplacian problems, using the :cpp:func:`invert_laplace`
     routines as preconditioners.

- invert / laplace

   - :doc:`invert_laplace.cxx<_breathe_autogen/file/invert__laplace_8cxx>` uses Fourier
      decomposition in :math:`z` combined with tri- and band-diagonal
      solvers in :math:`x` to solve Laplacian problems.

   - :doc:`laplacefactory.hxx<_breathe_autogen/file/laplacefactory_8hxx>`

   - :doc:`laplacefactory.cxx<_breathe_autogen/file/laplacefactory_8cxx>`

   - impls

      - serial\_tri

         - :doc:`serial_tri.hxx<_breathe_autogen/file/serial__tri_8hxx>`

         - :doc:`serial_tri.cxx<_breathe_autogen/file/serial__tri_8cxx>`

      - serial\_band

         - :doc:`serial_band.hxx<_breathe_autogen/file/serial__band_8hxx>`

         - :doc:`serial_band.cxx<_breathe_autogen/file/serial__band_8cxx>`

      - spt

         - :doc:`spt.hxx<_breathe_autogen/file/spt_8hxx>`

         - :doc:`spt.cxx<_breathe_autogen/file/spt_8cxx>`

      - pdd

         - :doc:`pdd.hxx<_breathe_autogen/file/pdd_8hxx>`

         - :doc:`pdd.cxx<_breathe_autogen/file/pdd_8cxx>`

- invert / parderiv

   -
     :doc:`invert_parderiv.cxx<_breathe_autogen/file/invert__parderiv_8cxx>`
     inverts a problem involving only parallel :math:`y`
     derivatives. Intended for use in some preconditioners.

   - :doc:`parderiv_factory.hxx<_breathe_autogen/file/parderiv__factory_8hxx>`

   - :doc:`parderiv_factory.cxx<_breathe_autogen/file/parderiv__factory_8cxx>`

   - impls

      - serial

         - :doc:`serial.cxx<_breathe_autogen/file/serial_8cxx>`

         - :doc:`serial.hxx<_breathe_autogen/file/serial_8hxx>`

      - cyclic

         - :doc:`cyclic.cxx<_breathe_autogen/file/cyclic_8cxx>`

         - :doc:`cyclic.hxx<_breathe_autogen/file/cyclic_8hxx>`

- :doc:`lapack_routines.cxx<_breathe_autogen/file/lapack__routines_8cxx>` supplies an
   interface to the LAPACK linear solvers, which are used by the
   ``invert_laplace`` routines.

- mesh

   - :doc:`boundary_factory.cxx<_breathe_autogen/file/boundary__factory_8cxx>` creates boundary
     condition operators which can then be applied to
     fields. Described in section [sec:BoundaryFactory].

   - :doc:`boundary_region.cxx<_breathe_autogen/file/boundary__region_8cxx>` implements a way
     to describe and iterate over boundary regions. Created by the
     mesh, and then used by boundary conditions. See
     section [sec:BoundaryRegion] for more details.

   - :doc:`boundary_standard.cxx<_breathe_autogen/file/boundary__standard_8cxx>` implements some
     standard boundary operations and modifiers such as ``Neumann``
     and ``Dirichlet``.

   - :doc:`difops.cxx<_breathe_autogen/file/difops_8cxx>` is a
     collection of differential operators on scalar fields. It uses
     the differential methods in :doc:`derivs.cxx<_breathe_autogen/file/derivs_8cxx>` and the metric tensor
     components in :cpp:class:`Mesh` to compute operators.

   - :doc:`interpolation.cxx<_breathe_autogen/file/interpolation_8cxx>` contains functions
     for interpolating fields

   - :doc:`mesh.cxx<_breathe_autogen/file/mesh_8cxx>` is the base
     class for the :cpp:class:`Mesh` object. Contains routines useful
     for all :cpp:class:`Mesh` implementations.

   - impls

      - :doc:`domain.cxx<_breathe_autogen/file/domain_8cxx>`

      - :doc:`domain.hxx<_breathe_autogen/file/domain_8hxx>`

      - :doc:`partition.cxx<_breathe_autogen/file/partition_8cxx>`

      - :doc:`partition.hxx<_breathe_autogen/file/partition_8hxx>`

      - bout

         - :doc:`boutmesh.cxx<_breathe_autogen/file/boutmesh_8cxx>`
           implements a mesh interface which is compatible with BOUT
           grid files.

         - :doc:`boutmesh.hxx<_breathe_autogen/file/boutmesh_8hxx>`

- physics

   - :doc:`gyro_average.cxx<_breathe_autogen/file/gyro__average_8cxx>`
      gyro-averaging operators

   - :doc:`smoothing.cxx<_breathe_autogen/file/smoothing_8cxx>`
     provides smoothing routines on scalar fields

   - :doc:`sourcex.cxx<_breathe_autogen/file/sourcex_8cxx>` contains
     some useful routines for creating sources and sinks in physics
     equations.

- precon

   - :doc:`jstruc.cxx<_breathe_autogen/file/jstruc_8cxx>` is an
     experimental code for preconditioning using PETSc

- solver

   - :doc:`solver.cxx<_breathe_autogen/file/solver_8cxx>` is the
     interface for all solvers

   - :doc:`solverfactory.cxx<_breathe_autogen/file/solverfactory_8cxx>` creates solver
     objects

   - :doc:`solverfactory.hxx<_breathe_autogen/file/solverfactory_8hxx>`

   - impls

      - cvode

         - :doc:`cvode.cxx<_breathe_autogen/file/cvode_8cxx>` is the
           implementation of :cpp:class:`Solver` which interfaces with
           the SUNDIALS CVODE library.

         - :doc:`cvode.hxx<_breathe_autogen/file/cvode_8hxx>`

      - ida

         - :doc:`ida.cxx<_breathe_autogen/file/ida_8cxx>` is the
           implementation which interfaces with the SUNDIALS IDA
           library

         - :doc:`ida.hxx<_breathe_autogen/file/ida_8hxx>`

      - petsc

         - :doc:`petsc.cxx<_breathe_autogen/file/petsc_8cxx>` is the
           interface to the PETSc time integration routines

         - :doc:`petsc.hxx<_breathe_autogen/file/petsc_8hxx>`

      - pvode

         - :doc:`pvode.cxx<_breathe_autogen/file/pvode_8cxx>`
           interfaces with the 1998 (pre-SUNDIALS) version of PVODE
           (which became CVODE).

         - :doc:`pvode.hxx<_breathe_autogen/file/pvode_8hxx>`

- sys

   - :doc:`boutcomm.cxx<_breathe_autogen/file/boutcomm_8cxx>`

   - :doc:`boutexception.cxx<_breathe_autogen/file/boutexception_8cxx>`
     is an exception class which are used for error handling

   - :doc:`comm_group.cxx<_breathe_autogen/file/comm__group_8cxx>`
     provides routines for non-blocking collective MPI
     operations. These are not available in MPI-2, though are planned
     for MPI-3.

   - :doc:`derivs.cxx<_breathe_autogen/file/derivs_8cxx>` contains
     basic derivative methods such as upwinding, central difference
     and WENO methods. These are then used by
     :doc:`difops.cxx<_breathe_autogen/file/difops_8cxx>`. Details are
     given in section [sec:derivatives].

   - :doc:`msg_stack.cxx<_breathe_autogen/file/msg__stack_8cxx>` is
     part of the error handling system. It maintains a stack of
     messages which can be pushed onto the stack at the start of a
     function, then removed (popped) at the end. If an error occurs or
     a segmentation fault is caught then this stack is printed out and
     can help to find errors.

   - :doc:`options.cxx<_breathe_autogen/file/options_8cxx>` provides
     an interface to the BOUT.inp option file and the command-line
     options.

   - :doc:`optionsreader.cxx<_breathe_autogen/file/optionsreader_8cxx>`

   - :doc:`output.cxx<_breathe_autogen/file/output_8cxx>`

   - :doc:`range.cxx<_breathe_autogen/file/range_8cxx>` Provides the
     RangeIterator class, used to iterate over a set of
     ranges. Described in section [sec:rangeiterator]

   - :doc:`stencils.cxx<_breathe_autogen/file/stencils_8cxx>` contains
     methods to operate on stencils which are used by differential
     methods.

   - :doc:`timer.cxx<_breathe_autogen/file/timer_8cxx>` a class for
     timing parts of the code like communications and file
     I/O. Described in section [sec:timerclass]

   - :doc:`utils.cxx<_breathe_autogen/file/utils_8cxx>` contains
     miscellaneous small useful routines such as allocating and
     freeing arrays.

   - options

      - :doc:`optionparser.hxx<_breathe_autogen/file/optionparser_8hxx>`

      - :doc:`options_ini.cxx<_breathe_autogen/file/options__ini_8cxx>`

      - :doc:`options_ini.hxx<_breathe_autogen/file/options__ini_8hxx>`

Data types
==========

The classes outlines in red in figure [fig:layout2] are data types
currently implemented in BOUT++.

``FieldData``
-------------

All BOUT++ data types implement a standard interface for accessing their
data, which is then used in communication and file I/O code. This
interface is in ``src/field/field_data.hxx``. The mandatory (pure
virtual) functions are:

::

    bool isReal(); // Returns true if field consists of real values
    bool is3D() const;   // True if variable is 3D
      
    int byteSize() const; // Number of bytes for a single point
    int realSize() const; // Number of reals (not implemented if not real)

To support file I/O there are also some additional functions which may
be implemented. A code can check if they are implemented by calling
``ioSupport``. If one of them is implemented then they all should be.

::

    bool  ioSupport();  // Return true if these functions are implemented
    const string getSuffix(int component) const; // For vectors e.g. "_x"
    void* getMark() const; // Store current settings (e.g. co/contra-variant)
    void  setMark(void *setting); // Return to the stored settings
    BoutReal* getData(int component); 
    void  zeroComponent(int component); // Set a component to zero

For twist-shift conditions, the optional function ``shiftZ`` is called
in the communication routines.

::

    void shiftZ(int jx, int jy, double zangle);

``Field``
---------

The two main types are ``Field2D``, and ``Field3D``. Their main
functions are to provide an easy way to manipulate data; they take care
of all memory management, and most looping over grid-points in algebraic
expressions. The 2D field implementation is relatively simple, but more
optimisations are used in the 3D field implementation because they are
much larger (factor of :math:`\sim 100`).

To handle time-derivatives, and enable expressions to be written in the
following form:

::

    ddt(Ni) = -b0xGrad_dot_Grad(phi, Ni);

fields (and vectors, see below) have a function:

::

    Field3D* timeDeriv();

which returns a pointer to the field holding the time-derivative of this
variable. This function ensures that this field is unique using a
singleton pattern.

``Vector``
----------

Vector classes build on the field classes, just using a field to
represent each component.

To handle time-derivatives of vectors, some care is needed to ensure
that the time-derivative of each vector component points to the same
field as the corresponding component of the time-derivative of the
vector:

::

    ddt(v.x) = ddt(v).x

``dcomplex``
------------

Several parts of the BOUT++ code involve FFTs and are therefore much
easier to write using complex numbers. Unfortunately, the C++ complex
library also tries to define a ``real`` type, which is already defined
by PVODE. Several work-arounds were tried, some of which worked on some
systems, but it was easier in the end to just implement a new class
``dcomplex`` to handle complex numbers.

Memory management
-----------------

This code has been thoroughly tested/debugged, and should only be
altered with great care, since just about every other part of BOUT++
depends on this code working correctly. Two optimisations used in the
data objects to speed up code execution are memory recycling, which
eliminates allocation and freeing of memory; and copy-on-change, which
minimises unnecessary copying of data.

Both of these optimisations are done “behind the scenes”, hidden from
the remainder of the code, and are illustrated in figure [fig:memory]:

.. figure:: figs/memory.pdf
   :alt: Memory handling in BOUT++

   Memory handling in BOUT++. Memory allocation and freeing is
   eliminated by recycling memory blocks, and assignments without
   changes (``A = B``) do not result in copying data, only pointers to
   the data. Both these optimisations are handled internally, and are
   invisible to the programmer.

The objects (A,B,C) accessed by the user in operations discussed in the
previous section act as an interface to underlying data (a,b). Memory
recycling can be used because all the scalar fields are the same size
(and vector fields are implemented as a set of 3 scalar fields). Each
class implements a global stack of available memory blocks. When an
object is assigned a value, it attempts to grab one of these memory
blocks, and if none are available then a new block is allocated. When an
object is destroyed, its memory block is not freed, but is put onto the
stack. Since the evaluation of the time-derivatives involves the same
set of operations each time, this system means that memory is only
allocated the first time the time-derivatives are calculated, after
which the same memory blocks are re-used. This eliminates the often slow
system calls needed to allocate and free memory, replacing them with
fast pointer manipulation.

Copy-on-change (reference counting) further reduces memory useage and
unnecessary copying of data. When one field is set equal to another
(e.g. ``Field3D A = B`` in figure [fig:memory]), no data is copied, only
the reference to the underlying data (in this case both A and B point to
data block a). Only when one of these objects is modified is a second
memory block used to store the different value. This is particularly
useful when returning objects from a routine. Usually this would involve
copying data from one object to another, and then destroying the
original copy. Using reference counting this copying is eliminated.

Derivatives
===========

This is probably the part of the code most people will want to alter,
and is in ``bout++/src/sys/derivs.cxx``. The main task of this module is
to map functions on fields like ``DDX`` to direction-independent
differential methods on stencils such as :math:`4^{th}`-order central
differencing. This mapping depends on global settings in ``BOUT.inp``
and is illustrated in figure [fig:diffOverview].

.. figure:: figs/diffOverview.pdf
   :alt: Overview of ``derivs`` module

   Overview of ``derivs`` module, mapping derivative functions on fields
   to direction-independent differential methods

Four kinds of differencing methods are supported

#. | First derivative ``DDX``, ``DDY``, ``DDZ``
   | Central differencing type schemes for first-order derivatives

#. | Second derivatives ``D2DX2``, ``D2DZ2``, ``D2DZ2``
   | Central differencing second derivatives e.g. for :math:`\nabla^2`

#. | Upwinding ``VDDX``, ``VDDY``, ``VDDZ``
   | Terms like :math:`\mathbf{v}\cdot\nabla`

#. | Flux methods ``FDDX``, ``FDDY``, ``FDDZ``
   | Flux conserving, limiting methods for terms like
     :math:`\nabla\cdot\left(\mathbf{v}f\right)`

The differencing methods themselves are independent on direction, and
have types defined in :doc:`derivs.cxx<_breathe_autogen/file/derivs_8cxx>`

::

    typedef BoutReal (*deriv_func)(stencil &); // f
    typedef BoutReal (*upwind_func)(stencil &, stencil &); // v, f

These operate on ``stencil`` objects. This class is in :doc:`stencils.hxx<_breathe_autogen/file/stencils_8hxx>`

::

    class stencil {
      public:
        int jx, jy, jz;  // Central location
        BoutReal c, p, m, pp, mm; // stencil 2 each side of the centre
        Overloaded operators
          =,+,-,*,/
        Functions
          min, max, abs
    };

The main purpose of this class is to store a 5-element stencil. To
simplify some code this class also has a bunch of overloaded operators
on BoutReals and other stencil objects. There are also some functions to
calculate things like absolute, minimum, and maximum values.

Lookup tables
-------------

To convert between short variable names (“C2”), long descriptions
(“2nd order Central Differencing”), ``DIFF_METHOD`` enums used to
specify methods at runtime (DIFF\_C2, defined in
:doc:`bout_types.hxx<_breathe_autogen/file/bout__types_8hxx>`), and
function pointers (``DDX_C2``), taking into account whether variables
are shifted or not, BOUT++ uses a set of lookup tables.

To find function pointers, tables of the following type are used:

::

    /// Translate between DIFF_METHOD codes, and functions
    struct DiffLookup {
      DIFF_METHOD method;
      deriv_func func;     // Single-argument differencing function
      upwind_func up_func; // Upwinding function
    };

Because the ``DiffLookup`` type contains a ``deriv_func`` and
``upwind_func`` pointer, it is used for all function lookup tables.
There is a separate table for each type of differencing method, so for
example the table of non-staggered upwinding methods is

::

    /// Upwinding functions lookup table
    static DiffLookup UpwindTable[] = { {DIFF_U1, NULL, VDDX_U1},
                        {DIFF_C2, NULL, VDDX_C2},
                        {DIFF_U4, NULL, VDDX_U4},
                        {DIFF_W3, NULL, VDDX_WENO3},
                        {DIFF_C4, NULL, VDDX_C4},
                        {DIFF_DEFAULT}};

The ``DIFF_DEFAULT`` at the end is used to terminate the array. These
tables are used by functions

::

    deriv_func lookupFunc(DiffLookup* table, DIFF_METHOD method);
    upwind_func lookupUpwindFunc(DiffLookup* table, DIFF_METHOD method);

which return the function pointer corresponding to the given method. If
the method isn’t in the table, then the first entry in the table is
used. These functions can be used at run-time to allow a user to specify
the method to use for specific operators.

When reading settings from the input file, they are specified as short
strings like “C2”, and a longer description of the method chosen should
be written to the output log. To do this, there is a name lookup table:

::

    /// Translate between short names, long names and DIFF_METHOD codes
    struct DiffNameLookup {
      DIFF_METHOD method;
      const char* label; // Short name
      const char* name;  // Long name
    };

    static DiffNameLookup DiffNameTable[] = { 
      {DIFF_U1, "U1", "First order upwinding"},
      {DIFF_C2, "C2", "Second order central"},
      {DIFF_W2, "W2", "Second order WENO"},
      {DIFF_W3, "W3", "Third order WENO"},
      {DIFF_C4, "C4", "Fourth order central"},
      {DIFF_U4, "U4", "Fourth order upwinding"},
      {DIFF_FFT, "FFT", "FFT"},
      {DIFF_DEFAULT}}; // Use to terminate the list

To search this table, there is the function

::

    DIFF_METHOD lookupFunc(DiffLookup *table, const string &label)

During initialisation, the lookup therefore works in two stages, shown
in figure [fig:diffLookup]. First the short description is turned into a
``DIFF_METHOD`` enum code, then this code is turned into a function
pointer.

.. figure:: figs/diffLookup.pdf
   :alt: Lookup tables for differential method

   Lookup tables for mapping between differential method labels, codes,
   descriptions and function pointers

Staggered grids
---------------

By default, all quantities in BOUT++ are defined at cell centre, and all
derivative methods map cell-centred quantities to cell centres.
Switching on staggered grid support in BOUT.inp:

::

    StaggerGrids = true

allows quantities to be defined on cell boundaries. Functions such as
``DDX`` now have to handle all possible combinations of input and output
locations, in addition to the possible derivative methods.

Several things are not currently implemented, which probably should be:

-  Only 3D fields currently have a cell location attribute. The location
   (cell centre etc) of 2D fields is ignored at the moment. The
   rationale for this is that 2D fields are assumed to be slowly-varying
   equilibrium quantities for which it won’t matter so much. Still,
   needs to be improved in future

-  Twist-shift and X shifting still treat all quantities as
   cell-centred.

-  No boundary condition functions yet account for cell location.

Currently, BOUT++ does not support values at cell corners; values can
only be defined at cell centre, or at the lower X,Y, or Z boundaries.
This is

Once staggered grids are enabled, two types of stencil are needed: those
which map between the same cell location (e.g. cell-centred values to
cell-centred values), and those which map to different locations (e.g.
cell-centred to lower X).

.. figure:: figs/diffStencils.pdf
   :alt: Stencils with cell-centred and lower shifted values

   Stencils with cell-centred (solid) and lower shifted values (open).
   Processor boundaries marked by vertical dashed line

Central differencing using 4-point stencil:

.. math::

   \begin{aligned}
   y &=& \left(9y_{-1/2} + 9y_{1/2} - y_{-3/2} - y_{3/2}\right) / 16 \\
   {\ensuremath{\frac{\partial y}{\partial x}}} &=& \left( 27y_{1/2} - 27y_{-1/2} - y_{3/2} + y_{-3/2}\right) / 24\Delta x \\
   \frac{\partial^2 y}{\partial x^2} &=& \left(y_{3/2} + y_{-3/2} - y_{1/2} - y_{-1/2}\right) / 2\Delta x^2\end{aligned}

+----------+-------------------+----------------------------------------------------------------+
| Input    | Output            | Actions                                                        |
+==========+===================+================================================================+
|          | Central stencil   |                                                                |
+----------+-------------------+----------------------------------------------------------------+
| CENTRE   | XLOW              | Lower staggered stencil                                        |
+----------+-------------------+----------------------------------------------------------------+
| XLOW     | CENTRE            | Upper staggered stencil                                        |
+----------+-------------------+----------------------------------------------------------------+
| XLOW     | Any               | Staggered stencil to CENTRE, then interpolate                  |
+----------+-------------------+----------------------------------------------------------------+
| CENTRE   | Any               | Central stencil, then interpolate                              |
+----------+-------------------+----------------------------------------------------------------+
| Any      | Any               | Interpolate to centre, use central stencil, then interpolate   |
+----------+-------------------+----------------------------------------------------------------+

Table: DDX actions depending on input and output locations. Uses first
match.

Laplacian inversion
===================

The Laplacian inversion code solves the equation:

.. math:: d\nabla^2_\perp x + \frac{1}{c}\nabla_\perp c\cdot\nabla_\perp x + a x = b

where :math:`x` and :math:`b` are 3D variables, whilst :math:`a`,
:math:`c` and :math:`d` are 2D variables. Several different algorithms
are implemented for Laplacian inversion, and they differ between
serial and parallel versions. Serial inversion can currently either be
done using a tridiagonal solver (Thomas algorithm), or a band-solver
(allowing :math:`4^{th}`-order differencing).

To support multiple implementations, a base class ``Laplacian`` is
defined in ``include/invert_laplace.hxx``. This defines a set of
functions which all implementations must provide:

.. code-block:: cpp

    class Laplacian {
     public:
      virtual void setCoefA(const Field2D &val) = 0;
      virtual void setCoefC(const Field2D &val) = 0;
      virtual void setCoefD(const Field2D &val) = 0;
     
      virtual const FieldPerp solve(const FieldPerp &b) = 0;
    }

At minimum, all implementations must provide a way to set coefficients,
and a solve function which operates on a single FieldPerp (X-Y) object
at once. Several other functions are also virtual, so default code
exists but can be overridden by an implementation.

For convenience, the ``Laplacian`` base class also defines a function to
calculate coefficients in a Tridiagonal matrix

.. code-block:: cpp

      void tridagCoefs(int jx, int jy, int jz, dcomplex &a, dcomplex &b,
                       dcomplex &c, const Field2D *ccoef = NULL,
                       const Field2D *d=NULL);

For the user of the class, some static functions are defined:

.. code-block:: cpp

      static Laplacian* create(Options *opt = NULL);
      static Laplacian* defaultInstance();

The create function allows new Laplacian implementations to be created,
based on options. To use the options in the ``[laplace]`` section, just
use the default:

.. code-block:: cpp

      Laplacian* lap = Laplacian::create();

The code for the ``Laplacian`` base class is in
``src/invert/laplace/invert_laplace.cxx``. The actual creation of new
Laplacian implementations is done in the ``LaplaceFactory`` class,
defined in ``src/invert/laplace/laplacefactory.cxx``. This file includes
all the headers for the implementations, and chooses which one to create
based on the “type” setting in the input options. This factory therefore
provides a single point of access to the underlying Laplacian inversion
implementations.

Each of the implementations is in a subdirectory of
``src/invert/laplace/impls`` and is discussed below.

Serial tridiagonal solver
-------------------------

This is the simplest implementation, and is in
``src/invert/laplace/impls/serial_tri/``

Serial band solver
------------------

This is band-solver which performs a :math:`4^{th}`-order inversion.
Currently this is only available when ``NXPE=1``; when more than one
processor is used in :math:`x`, the Laplacian algorithm currently
reverts to :math:`3^{rd}`-order.

SPT parallel tridiagonal
------------------------

This is a reference code which performs the same operations as the
serial code. To invert a single XZ slice (``FieldPerp`` object), data
must pass from the innermost processor (``mesh->PE_XIND = 0``) to the
outermost ``mesh->PE_XIND = mesh->NXPE-1`` and back again.

Some parallelism is achieved by running several inversions
simultaneously, so while processor 1 is inverting Y=0, processor 0 is
starting on Y=1. This works ok as long as the number of slices to be
inverted is greater than the number of X processors
(``MYSUB > mesh->NXPE``). If ``MYSUB < mesh->NXPE`` then not all
processors can be busy at once, and so efficiency will fall sharply.
Figure [fig:par\_laplace] shows the useage of 4 processors inverting a
set of 3 poloidal slices (i.e. MYSUB=3)

.. figure:: figs/par_laplace.pdf
   :alt: Parallel Laplacian inversion

   Parallel Laplacian inversion with MYSUB=3 on 4 processors. Red
   periods are where a processor is idle - in this case about 40% of the
   time

PDD algorithm
-------------

This is the Parallel Diagonally Dominant (PDD) algorithm. It’s very
fast, but achieves this by neglecting some cross-processor terms. For
ELM simulations, it has been found that these terms are important, so
this method is not usually used.

Mesh
====

The mesh is used in pretty much all parts of the code, and deals with
things like the geometry of the mesh (metric tensors etc.), and how the
mesh is divided between processors (communications). The Mesh class
(``include/mesh.hxx``) defines an interface, and there are currently two
implementations:

-  ``BoutMesh`` (``src/mesh/boutmesh.cxx``) which is backwards
   compatible with the BOUT and BOUT-06 codes. This is a logically
   rectangular mesh so the number of radial points (x) can’t change in
   the poloidal direction (y).

-  ``QuiltMesh`` (``src/mesh/quiltmesh.cxx``) is a more general mesh
   under development (i.e. **not** recommended except for testing). The
   primary advantage is that it allows the number of points in x to vary
   between regions so the number of radial grid points in the core no
   longer needs to be the same as the number in the private flux
   regions.

Grid data sources
-----------------

All data sources inherit from ``GridDataSource``, defined in
:doc:`grid.hxx<_breathe_autogen/file/griddata_8hxx>` at line 43. They must supply a method to test if a variable
exists:

::

    bool GridDataSource::hasVar(const char *name);

a method to get the size of the variable

::

    vector<int> GridDataSource::getSize(const char *name);

To fetch data, first the (x,y,z) origin must be set:

::

    bool GridDataSource::setOrigin(int x = 0, int y = 0, int z = 0);

and then use methods to fetch integers or reals:

::

    bool GridDataSource::fetch(int *var, const string &name, int lx = 1, int ly = 0, int lz = 0);
    bool GridDataSource::fetch(BoutReal *var, const string &name, int lx = 1, int ly = 0, int lz = 0);

In addition, GridDataSource implementations can have methods which
should be called before and after variables are accessed:

::

    void GridDataSource::open(const char *name = NULL);
    void GridDataSource::close();

Loading a mesh
--------------

To load in a mesh from a file or other source, there are the commands:

::

    int addSource(GridDataSource);   // Add a data source
    int load();                      // Load from added data sources
    int load(GridDataSource);        // Load from specified data source

all of which return an error code (0 if successful). ``addSource`` is
used to add a set of input data sources which inherit from
``GridDataSource``. ``load()`` loads the mesh from these sources,
querying each data source in turn for the required variables (in the
order in which they were added). ``load(GridDataSource)`` loads the mesh
from only the supplied data source.

In :doc:`bout++.cxx<_breathe_autogen/file/bout_09_09_8cxx>`, this is used to initialise the mesh:

::

    mesh->addSource(new GridFile(data_format(grid_name), grid_name));
    if(mesh->load()) {
      output << "Failed to read grid. Aborting\n";
      return 1;
    }

which creates a ``GridFile`` object based on the data format of the grid
file name, then adds that as a source of data for Mesh.

For post-processing of the results, it’s useful to have mesh quantities
in the dump files along with the results. To do this, there’s the
function

::

    void outputVars(Datafile &file); // Add mesh vars to file

which is called during BOUT++ initialisation.

Implementation: BoutMesh
~~~~~~~~~~~~~~~~~~~~~~~~

BoutMesh class uses the BOUT indices (which trace back to UEDGE):

::

    int ixseps1, ixseps2, jyseps1_1, jyseps2_1, jyseps1_2, jyseps2_2;

``ixseps1`` and ``ixseps2`` give the X location of the separatrices, and
are equal in the case of single-null configurations. The indexing is
such that all points ``0 <= x < ixseps1`` are inside the separatrix,
whilst ``ixseps1 <= x < ngx`` are outside.

Implementation: QuiltMesh
~~~~~~~~~~~~~~~~~~~~~~~~~

Index ranges
------------

The Mesh class includes several public members which describe the size
of the mesh, and are used all over BOUT++ to loop over variables:

::

    /// Size of the mesh on this processor including guard/boundary cells
    int ngx, ngy, ngz;  
    /// Local ranges of data (inclusive), excluding guard cells
    int xstart, xend, ystart, yend;

Getting data
------------

The ``load()`` code above needs to read data for the mesh, and physics
codes usually need to read their initial profiles during initialisation.
To do this, Mesh provides an overloaded function ``get``:

::

    int get(var, const char *name); // Request data from mesh file

where ``var`` can be just about any BOUT++ datatype (``Field2D``,
``Vector3D`` etc.).

Implementation: BoutMesh
~~~~~~~~~~~~~~~~~~~~~~~~

For integers and BoutReals, the implementation is fairly trivial. Uses
the Mesh protected functions to find a data source and read data from
it.

::

    GridDataSource* s = findSource(name);  // Find a source of data
    s->open(name);                          // Open the source
    bool success = s->fetch(&ival, name);   // Get the data
    s->close();                             // Close the source

To read 2D and 3D fields, the branch-cuts need to be taken into account.

Communications
--------------

The most common type of communication is to just exchange all guard
cells with neighboring processors. Mesh provides the following commands
for doing this:

::

    int communicate(FieldData, ...); // Communicate one or more fields
    int communicate(FieldGroup);     // Communicate a group of fields
    int communicate(FieldData);      // Returns error code
    comm_handle send(FieldGroup);    // Send data
    int wait(comm_handle);           // Receive data

``communicate(FieldData)`` can (currently) be used to communicate up to
4 variables together, and makes the code quite clear. For example in
``examples/DriftInstability/2fluid.cxx`` around line 360:

::

    // Need to communicate jpar
    mesh->communicate(jpar);

Since this uses the ``FieldData`` interface like Datafile, this can be
used to communicate all BOUT++ field data types. The limit of 4 is
because the C-style ``varargs`` system doesn’t work with “non POD”
variables, i.e. classes. To communicate a larger number of variables,
create a ``FieldGroup`` object to group fields together, then
communicate them all together:

::

    FieldGroup comgrp;  // Group of variables for communication
    Field3D P;
    Vector3D V;

    comgrp.add(P); // Add the variables
    comgrp.add(V); // Usually done in physics_init

    mesh->communicate(comgrp); // Communicate in physics_run

If you want to overlap communications with calculations then use the
``send`` and ``wait`` functions instead of ``communicate``.

::

    comm_handle ch = mesh->send(comgrp); // Start the communications
    // Calculations which don't need variables in comgrp
    wait(ch); // Wait for all communications to finish

Implementation: BoutMesh
~~~~~~~~~~~~~~~~~~~~~~~~

In BoutMesh, the communication is controlled by the variables

::

    int UDATA_INDEST, UDATA_OUTDEST, UDATA_XSPLIT;
    int DDATA_INDEST, DDATA_OUTDEST, DDATA_XSPLIT;
    int IDATA_DEST, ODATA_DEST;

In the Y direction, each boundary region (**U**\ p and **D**\ own in Y)
can be split into two, with ``0 <= x < UDATA_XSPLIT`` going to the
processor index ``UDATA_INDEST``, and ``UDATA_INDEST <= x < ngx`` going
to ``UDATA_OUTDEST``. Similarly for the Down boundary. Since there are
no branch-cuts in the X direction, there is just one destination for the
**I**\ nner and **O**\ uter boundaries. In all cases a negative
processor number means that there’s a domain boundary.

X communications
----------------

For parallel Laplacian inversions, communication is needed in the X
direction only, and involves quantities which are not in Fields.

::

    bool firstX();  // True if at the inner X boundary
    bool lastX();   // True if at the outer X boundary
    int NXPE, PE_XIND; // Number of processors in X, and X processor index
    int sendXOut(BoutReal *buffer, int size, int tag);
    sendXIn(BoutReal *buffer, int size, int tag);
    comm_handle irecvXOut(BoutReal *buffer, int size, int tag);
    comm_handle irecvXIn(BoutReal *buffer, int size, int tag);

The variables ``NXPE`` and ``PE_XIND`` shouldn’t really be there, but
are currently needed because the SPT algorithm in :doc:`invert_laplace.cxx<_breathe_autogen/file/invert__laplace_8cxx>`
needs to know when it’s going to be next and so keep track of which
processor number is currently working. This logic to pass a problem
along a chain in X should really be moved into Mesh.

Y-Z surface communications
--------------------------

Some operations (like parallel inversions in
``bout++/src/invert/invert_parderiv.cxx``) need to be performed on Y-Z
surfaces, i.e. slices at constant X. This needs to be able to handle
open and closed surfaces, and that closed surfaces may need a shift in
the Z direction to match one end onto the other (a twist-shift
condition).

The simplest operation is to average a quantity over Y:

::

    const Field2D averageY(const Field2D &f); // Average in Y

Currently this is only implemented for 2D fields. More generally a set
of FieldData objects could be used.

To test if a particular surface is closed, there is the function

::

    bool surfaceClosed(int jx, BoutReal &ts); // Test if a surface is closed, and if so get the twist-shift angle

The most general way to access data on surfaces is to use an iterator,
which can be created using:

::

    SurfaceIter* iterateSurfaces();

This then allows looping over the surfaces in the usual way

::

    for(surf->first(); !surf->isDone(); surf->next()) {
      ...
    }

**NB**: This iterator splits the surfaces between processors, so each
individual processor will iterate over a different set of surfaces. This
is to allow automatic load balancing when gathering and scattering data
from an entire surface onto one processor using:

::

    surf->gather(FieldData, BoutReal *recvbuffer);
    surf->scatter(BoutReal *sendbuffer, Field result);

The buffer is assumed to be large enough to hold all the data. To get
the number of points in Y for this surface, use

::

    int ysize = surf->ysize();

To test if the surface is closed, there’s the test

::

    bool surf->closed(BoutReal &ts)

which returns true if the surface is closed, along with the twist-shift
angle.

Initial profiles
----------------

The initial profiles code needs to construct a solution which is smooth
everywhere, with a form of perturbation specified in the input file for
each direction. In order to do this, it needs a continuous function to
use as an index. This is supplied by the functions:

::

    BoutReal GlobalX(int jx); // Continuous X index between 0 and 1
    BoutReal GlobalY(int jy); // Continuous Y index (0 -> 1)

which take a local x or y index and return a globally continuous x or y
index.

Differencing
------------

The mesh spacing is given by the public members

::

    // These used for differential operators 
    Field2D dx, dy;
    Field2D d2x, d2y;    // 2nd-order correction for non-uniform meshes     
    BoutReal zlength, dz;    // Derived from options (in radians)

Metrics
-------

The contravariant and covariant metric tensor components are public
members of ``Mesh``:

::

    // Contravariant metric tensor (g^{ij})
    Field2D g11, g22, g33, g12, g13, g23; // These are read in grid.cxx

    // Covariant metric tensor
    Field2D g_11, g_22, g_33, g_12, g_13, g_23;

    int calcCovariant();     // Invert contravatiant metric to get covariant
    int calcContravariant(); // Invert covariant metric to get contravariant

If only one of these sets is modified by an external code, then
``calc_covariant`` and ``calc_contravariant`` can be used to calculate
the other (uses Gauss-Jordan currently).

From the metric tensor components, Mesh calculates several other useful
quantities:

::

    int jacobian(); // Calculate J and Bxy
    Field2D J; // Jacobian
    Field2D Bxy; // Magnitude of B = nabla z times nabla x

    /// Calculate differential geometry quantities from the metric tensor
    int geometry();

    // Christoffel symbol of the second kind (connection coefficients)
    Field2D G1_11, G1_22, G1_33, G1_12, G1_13;
    Field2D G2_11, G2_22, G2_33, G2_12, G2_23;
    Field2D G3_11, G3_22, G3_33, G3_13, G3_23;
      
    Field2D G1, G2, G3;

These quantities are public and accessible everywhere, but this is
because they are needed in a lot of the code. They shouldn’t change
after initialisation, unless the physics model starts doing fancy things
with deforming meshes.

Miscellaneous
-------------

There are some public members of Mesh which are there for some specific
task and don’t really go anywhere else (yet).

To perform radial derivatives in tokamak geometry, interpolation is
needed in the Z direction. This is done by shifting in Z by a phase
factor, performing the derivatives, then shifting back. The following
public variables are currently used for this:

::

    bool ShiftXderivs; // Use shifted X derivatives
    int  ShiftOrder;   // Order of shifted X derivative interpolation
    Field2D zShift;    // Z shift for each point (radians)
      
    Field2D ShiftTorsion; // d <pitch angle> / dx. Needed for vector differentials (Curl)
    Field2D IntShiftTorsion; // Integrated shear (I in BOUT notation)
    bool IncIntShear; // Include integrated shear (if shifting X)

::

    int  TwistOrder;   // Order of twist-shift interpolation

This determines what order method to use for the interpolation at the
twist-shift location, with ``0`` meaning FFT during communication. Since
this must be 0 at the moment it’s fairly redundant and should be
removed.

A (currently experimental) feature is

::

    bool StaggerGrids;    ///< Enable staggered grids (Centre, Lower). Otherwise all vars are cell centred (default).

Boundary conditions
===================

The boundary condition system needs to be very flexible in order to
handle:

-  Meshes which can divide up the boundary into an arbitrary number of
   regions, giving each one a label. For example in BoutMesh (specific
   to tokamaks), the boundary regions are labelled “core”, “sol”, “pf”
   and “target”.

-  Each variable can have a different boundary condition in each region.
   It should be possible to have a global setting “all variables have
   dirichlet conditions on all boundaries”, which is over-ridden by more
   specific settings such as “All variables have neumann conditions on
   the inner x boundaries”, and finally to “variable ’Ni’ has laplacian
   boundary conditions in the ’sol’ regions”

-  Boundary conditions can be modified to be “relaxing”. This means that
   rather than enforcing a strict boundary condition, it’s a mixture of
   zero-gradient in the time-derivative combined with a damping
   (relaxation) towards the desired boundary condition. This can help
   improve the numerics of turbulence simulations.

-  Users should be able to implement their own boundary conditions, and
   add them to the system at run-time without modifying the core code.

-  After ``physics_init``, a boundary condition must be applied to the
   variables. During a simulation (at the end of ``physics_run``),
   boundary conditions need to be applied to the time-derivatives. The
   boundary system should ensure that these conditions are consistent.

Variable initialisation
=======================

Solver
======

The solver is the interface between BOUT++ and the time-integration code
such as SUNDIALS. All solvers implement the ``Solver`` class interface
(see ``src/solver/generic_solver.hxx``).

First all the fields which are to be evolved need to be added to the
solver. These are always done in pairs, the first specifying the field,
and the second the time-derivative:

::

    void add(Field2D &v, Field2D &F_v, const char* name);

This is normally called in the ``physics_init`` initialisation routine.
Some solvers (e.g. IDA) can support constraints, which need to be added
in the same way as evolving fields:

::

    bool constraints();
    void constraint(Field2D &v, Field2D &C_v, const char* name);

The ``constraints()`` function tests whether or not the current solver
supports constraints. The format of ``constraint(...)`` is the same as
``add``, except that now the solver will attempt to make ``C_v`` zero.
If ``constraint`` is called when the solver doesn’t support them then an
error should occur.

If the physics model implements a preconditioner or Jacobian-vector
multiplication routine, these can be passed to the solver during
initialisation:

::

    typedef int (*PhysicsPrecon)(BoutReal t, BoutReal gamma, BoutReal delta);
    void setPrecon(PhysicsPrecon f); // Specify a preconditioner
    typedef int (*Jacobian)(BoutReal t);
    void setJacobian(Jacobian j); // Specify a Jacobian

If the solver doesn’t support these functions then the calls will just
be ignored.

Once the problem to be solved has been specified, the solver can be
initialised using:

::

    int init(rhsfunc f, int argc, char **argv, bool restarting, int nout, BoutReal tstep);

which returns an error code (0 on success). This is currently called in
:doc:`bout++.cxx<_breathe_autogen/file/bout_09_09_8cxx>`:

::

    if(solver.init(physics_run, argc, argv, restart, NOUT, TIMESTEP)) {
      output.write("Failed to initialise solver. Aborting\n");
      return(1);
    }

which passes the (physics module) RHS function ``physics_run`` to the
solver along with the number and size of the output steps.

To run the solver using the (already supplied) settings, there is the
function:

::

    typedef int (*MonitorFunc)(BoutReal simtime, int iter, int NOUT);
    int run(MonitorFunc f);

Implementation: PVODE
---------------------

Implementation: IDA
-------------------

Implementation: PETSc
---------------------

File I/O
========

BOUT++ needs to deal with binary format files to read the grid; read and
write restart restart files; and write dump files. The two parts of the
code which need to read and write data are therefore the grid routines
(:doc:`grid.hxx<_breathe_autogen/file/griddata_8hxx>`), and the ``Datafile`` class
(:doc:`datafile.hxx<_breathe_autogen/file/datafile_8hxx>` and :doc:`datafile.cxx<_breathe_autogen/file/datafile_8cxx>`). All other parts which need to
read or write data go through these methods.

Several different file formats are commonly used, such as HDF, HDF5, and
netCDF. For historical reasons (inherited from BOUT), BOUT++ originally
used the Portable Data Binary (PDB) format developed at LLNL.  [3]_ To
separate the basic file format functions from the higher level grid and
Datafile classes, these use an abstract class ``DataFormat``. Any class
which implements the functions listed in :doc:`dataformat.hxx<_breathe_autogen/file/dataformat_8hxx>` can
therefore be passed to grid or datafile. This makes implementing a new
file format, and switching between formats at run-time, relatively
straightforward.

Access to data in files is provided using a Bridge pattern: The
``Datafile`` class provides an interface to the rest of the code to read
and write variables, whilst file formats implement the ``Dataformat``
interface.

::

    class Datafile {
     public:
      Datafile();
      Datafile(DataFormat *format);
      ~Datafile();
      
      /// Set the file format by passing an interface class
      void setFormat(DataFormat *format);

      void setLowPrecision(); ///< Only output floats

      void add(var, const char *name, int grow = 0);

      int read(const char *filename, ...);
      int write(const char *filename, ...);
      int append(const char *filename, ...);
      bool write(const string &filename, bool append=false);

      /// Set this to false to switch off all data writing
      static bool enabled;
    };

The important bits of the DataFormat interface are:

::

    class DataFormat {
     public:
      bool openr(const char *name);
      bool openw(const char *name, bool append=false);
      
      bool is_valid();
      
      void close();
      
      const char* filename();

      const vector<int> getSize(const char *var);
      const vector<int> getSize(const string &var);

      // Set the origin for all subsequent calls
      bool setOrigin(int x = 0, int y = 0, int z = 0); 
      bool setRecord(int t); // negative -> latest
      
      // Read / Write simple variables up to 3D

      bool read(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0);
      bool read(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0);

      bool write(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0);
      bool write(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0);

      // Read / Write record-based variables

      bool read_rec(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0);
      bool read_rec(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0);

      bool write_rec(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0);
      bool write_rec(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0);

      // Optional functions
      
      void setLowPrecision();
    };

Miscellaneous
=============

Other small modules which don’t really fit into any system, but are
needed.

Printing messages
-----------------

Iterating over ranges
---------------------

The boundary of a processor’s domain may consist of a set of disjoint
ranges, so the mesh needs a clean way to tell any code which depends on
the boundary how to iterate over it. The ``RangeIterator`` class in
``include/bout/sys/range.hxx`` and ``src/sys/range.cxx`` provides this.

RangeIterator can represent a single continuous range, constructed by
passing the minimum and maximum values.

::

    RangeIterator it(1,4);  // Range includes both end points
    for(it.first(); !it.isDone(); it.next())
      cout << it.ind; // Prints 1234

A more canonical C++ style is also supported, using overloaded ``++``,
``*``, and ``!=`` operators:

::

    for(it.first(); it != RangeIterator::end(); it++)
      cout << *it; // Prints 1234

where ``it++`` is the same as ``it.next()``, and ``*it`` the same as
``it.ind``.

To iterate over several ranges, ``RangeIterator`` can be constructed
with the next range as an argument:

::

    RangeIterator it(1,4, RangeIterator(6,9));
    for(it.first(); it != RangeIterator::end(); it++)
      cout << *it; // Prints 12346789

and these can be chained together to an arbitrary depth.

To support statements like

::

    for(RangeIterator it = mesh->iterateBndryLowerY(); !it.isDone(); it++)
      ...

the initial call to ``first()`` is optional, and everything is
initialised in the constructor.

.. [2]
   generated using Al Danial’s cloc

.. [3]
   Support for PDB files was removed in BOUT++ 4.0.0
