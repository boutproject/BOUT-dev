"""
Routines for watching files for changes

"""
from __future__ import print_function
from builtins import zip

import time
import os


def watch(files, timeout=None, poll=2):
    """Watch a given file or collection of files until one changes. Uses
    polling.

    Parameters
    ----------
    files : str or list of str
        Name of one or more files to watch
    timeout : int, optional
        Timeout in seconds (default is no timeout)
    poll : int, optional
        Polling interval in seconds (default: 2)

    Returns
    -------
    str
        The name of the first changed file,
        or None if timed out before any changes

    Examples
    --------

    To watch one file, timing out after 60 seconds:

    >>> watch('file1', timeout=60)

    To watch 2 files, never timing out:

    >>> watch(['file1', 'file2'])

    Author: Ben Dudson <benjamin.dudson@york.ac.uk>

    """

    # Get modification time of file(s)
    try:
        if hasattr(files, '__iter__'):
            # Iterable
            lastmod = [ os.stat(f).st_mtime for f in files ]
            iterable = True
        else:
            # Not iterable -> just one file
            lastmod = os.stat(files).st_mtime
            iterable = False
    except:
        print("Can't test modified time. Wrong file name?")
        raise

    start_time = time.time()
    running = True
    while running:
        sleepfor = poll
        if timeout:
            # Check if timeout will be reached before next poll
            if time.time() - start_time + sleepfor > timeout:
                # Adjust time so that finish at timeout
                sleepfor = timeout - (time.time() - start_time)
                running = False # Stop after next test
                
        time.sleep(sleepfor)

        if iterable:
            for last_t, f in zip(lastmod, files):
                # Get the new modification time
                t = os.stat(f).st_mtime
                if t > last_t + 1.0:  # +1 to reduce risk of false alarms
                    # File has been modified
                    return f
        else:
            t = os.stat(files).st_mtime
            if t > lastmod + 1.0:
                return files
    return None
