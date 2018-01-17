.. Use bash as the default language for syntax highlighting in this file
.. highlight:: console

Compiling and running under AIX
===============================

Most development and running of BOUT++ is done under Linux, with the
occasional FreeBSD and OSX. The configuration scripts are therefore
heavily tested on these architectures. IBM’s POWER architecture however
runs AIX, which has some crucial differences which make compiling a
pain.

-  Under Linux/BSD, it’s usual for a Fortran routine ``foo`` to appear
   under C as ``foo_``, whilst under AIX the name is unchanged

-  MPI compiler scripts are usually given the names ``mpicc`` and either
   ``mpiCC`` or ``mpicxx``. AIX uses ``mpcc`` and ``mpCC``.

-  Like BSD, the ``make`` command isn’t compatible with GNU make, so you
   have to run ``gmake`` to compile everything.

-  The POWER architecture is big-endian, different to the little endian
   Intel and AMD chips. This can cause problems with binary file
   formats.

.. _sec-sundials:

SUNDIALS
--------

To compile SUNDIALS, use::

    $ export CC=cc
    $ export CXX=xlC
    $ export F77=xlf
    $ export OBJECT_MODE=64
    $ ./configure --prefix=$HOME/local/ --with-mpicc=mpcc --with-mpif77=mpxlf CFLAGS=-maix64

You may get an error message like::

    make: Not a recognized flag: w

This is because the AIX ``make`` is being used, rather than ``gmake``.
The easiest way to fix this is to make a link to ``gmake`` in your local
bin directory::

    $ ln -s /usr/bin/gmake $HOME/local/bin/make

Running ``which make`` should now point to this ``local/bin/make``, and
if not then you need to make sure that your bin directory appears first
in the ``PATH``::

    $ export PATH=$HOME/local/bin:$PATH

If you see an error like this::

    ar: 0707-126 ../../src/sundials/sundials_math.o is not valid with the current object file mode.
            Use the -X option to specify the desired object mode.

then you need to set the environment variable ``OBJECT_MODE``::

    $ export OBJECT_MODE=64

Configuring BOUT++, you may get the error::

    configure: error: C compiler cannot create executables

In that case, you can try using::

    $ ./configure CFLAGS="-maix64"

When compiling, you may see warnings::

    xlC_r: 1501-216 (W) command option -64 is not recognized - passed to ld

At this point, the main BOUT++ library should compile, and you can try
compiling one of the examples.::

    ld: 0711-317 ERROR: Undefined symbol: .NcError::NcError(NcError::Behavior)
    ld: 0711-317 ERROR: Undefined symbol: .NcFile::is_valid() const
    ld: 0711-317 ERROR: Undefined symbol: .NcError::~NcError()
    ld: 0711-317 ERROR: Undefined symbol: .NcFile::get_dim(const char*) const

This is probably because the NetCDF libraries are 32-bit, whilst BOUT++
has been compiled as 64-bit. You can try compiling BOUT++ as 32-bit::

    $ export OBJECT_MODE=32
    $ ./configure CFLAGS="-maix32"
    $ gmake

If you still get undefined symbols, then go back to 64-bit, and edit
make.config, replacing ``-lnetcdf_c++`` with -lnetcdf64\_c++, and
``-lnetcdf`` with -lnetcdf64. This can be done by running::

    $ sed 's/netcdf/netcdf64/g' make.config > make.config.new
    $ mv make.config.new make.config

