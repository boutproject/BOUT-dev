Cori configuration scripts
==========================


Environment
-----------

The `setup-*` scripts configure the environment, ensuring that the correct
modules are loaded. Note that as the system is upgraded these may become
out of date and need to be updated.

To use these scripts a Bash shell, use `source`:

    $ source setup-env-cgpu.sh

or if using C shell:

    > source setup-env-cgpu.csh

BOUT++ configuration
--------------------

The `config-*` scripts pass arguments to cmake, which have been designed to work
in the environment setup above.

From the BOUT-dev root directory:

    ./config/cori/config-bout-cgpu.sh

