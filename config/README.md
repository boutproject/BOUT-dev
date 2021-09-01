Configuration scripts
=====================

The CMake and autotools (configure/make) scripts supplied with BOUT++
should be able to automatically find and configure BOUT++ in most
cases. Where a complex configuration is desired, for example including
many dependencies (esp. complex dependencies like PETSc), or compiling
for GPUs, configuration can be quite complex.

The files in this directory are intended to be convenient shortcuts for
configuration on particular machines. Where there are many scripts, these
are put into sub-directories (e.g. "cori" and "lassen"). 

Environment
-----------

Scripts which set up the environment, for example loading and unloading
modules, start with `setup` or `setup-env`. These are typically modifying
shell environments and so should be invoked with `source`.

BOUT++ configuration
--------------------

The wrappers around CMake (or configure) start with `config` or `config-bout`.
These are shell scripts which can be run without `source`. 

