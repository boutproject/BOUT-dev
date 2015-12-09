#!/usr/bin/env python

"""Post processing routine which shows the data"""

from boutdata.collect import collect
from boututils.showdata import showdata

# All post processing functions called by bout_runners must accept the
# first argument from bout_runners (called 'folder' in
# __call_post_processing_function)
def show_the_data(path, t=None, x=None, y=None, z=None):
    """Function which plots the data.

    Input
    path     -   the path of a run
    t        -   the desired t slice of showdata
    x        -   the desired x slice of showdata
    y        -   the desired y slice of showdata
    z        -   the desired z slice of showdata
    """

    print("Showing data from " + path)
    n = collect('n', xguards=False, yguards=False, path=path, info=False)

    # Show the data
    showdata(n[t,x,y,z])
