.. _sec-logging:

Logging output
==============

Logging should be used to report simulation progress, record information,
and warn about potential problems. BOUT++ includes a simple logging facility
which supports both C printf and C++ iostream styles. For example:

::

   output.write("This is an integer: %d, and this a real: %e\n", 5, 2.0)
   
   output << "This is an integer: " << 5 << ", and this a real: " << 2.0 << endl;

Messages sent to ``output`` on processor 0 will be printed to console and saved to
``BOUT.log.0``. Messages from all other processors will only go to their log files,
``BOUT.log.#`` where ``#`` is the processor number.

**Note**: If an error occurs on a processor other than processor 0, then the
error message will usually only be in the log file, not printed to console. If BOUT++
crashes but no error message is printed, try looking at the ends of all log files:

.. code-block:: bash

   $ tail BOUT.log.*


For finer control over which messages are printed, several outputs are available,
listed in the table below.

===================   =================================================================
Name                  Useage
===================   =================================================================
``output_debug``      For highly verbose output messages, that are normally not needed.
                      Needs to be enabled with a compile switch
``output_info``       For infos like what options are used
``output_progress``   For infos about the current progress
``output_warn``       For warnings
``output_error``      For errors
===================   =================================================================


Controlling logging level
-------------------------

By default all of the outputs except ``output_debug`` are saved to log and printed
to console (processor 0 only).

To reduce the volume of outputs the command line argument ``-q`` (quiet) reduces
the output level by one, and ``-v`` (verbose) increases it by one. Running with ``-q``
in the command line arguments suppresses the ``output_info`` messages, so that they
will not appear in the console or log file. Running with ``-q -q`` suppresses everything
except ``output_warn`` and ``output_error``. 

To enable the ``output_debug`` messages, first configure BOUT++ with debug messages enabled
by adding ``-DDEBUG_ENABLED`` to ``BOUT_FLAGS`` in ``make.config`` and then recompiling
with ``make clean; make``. When running BOUT++ add a "-v" flag to see ``output_debug`` messages.



