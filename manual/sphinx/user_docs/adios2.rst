.. _sec-adios2:

ADIOS2 support
==============

This section summarises the use of `ADIOS2 <https://adios2.readthedocs.io/>`_ in BOUT++.

Installation
------------

The easiest way to configure BOUT++ with ADIOS2 is to tell CMake to download and build it
with this flag::

  -DBOUT_DOWNLOAD_ADIOS=ON

The ``master`` branch will be downloaded from `Github <https://github.com/ornladios/ADIOS2>`_,
configured and built with BOUT++.

Alternatively, if ADIOS is already installed then the following flags can be used::

  -DBOUT_USE_ADIOS=ON -DADIOS2_ROOT=/path/to/adios2

Output files
------------

The output (dump) files are controlled with the root ``output`` options.
By default the output format is NetCDF, so to use ADIOS2 instead set
the output type in BOUT.inp::

  [output]
  type = adios

or on the BOUT++ command line set ``output:type=adios``. The default
prefix is "BOUT.dmp" so the ADIOS file will be called "BOUT.dmp.bp". To change this,
set the ``output:prefix`` option.

Restart files
-------------

The restart files are contolled with the root ``restart_files`` options,
so to read and write restarts from an ADIOS dataset, put in BOUT.inp::

  [restart_files]
  type = adios

