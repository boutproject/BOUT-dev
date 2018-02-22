#!/usr/bin/env python

"""Post processing routine which shows the data"""

from boutdata.collect import collect
from boututils.showdata import showdata

# All post processing functions called by bout_runners must accept the
# first argument from bout_runners (called 'folder' in
# __call_post_processing_function)
def show_the_data(paths, t=None, x=None, y=None, z=None, **kwargs):
    """Function which plots the data.

    Parameters
    ----------
    paths : tuple
        The paths of the runs
    t : slice
        The desired t slice of showdata
    x : slice
        The desired x slice of showdata
    y : slice
        The desired y slice of showdata
    z : slice
        The desired z slice of showdata
    **kwargs : key word arguments
        Not used here, but acts like a "dumpster" for additional keyword
        arguments
    """

    for path in paths:
        print("Showing data from {}".format(path))
        n = collect('n', xguards=False, yguards=False, path=path, info=False)

        # Show the data
        showdata(n[t,x,y,z])
