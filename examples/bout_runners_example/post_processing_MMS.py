#!/usr/bin/env python

"""Post processing which performs MMS"""

from boutdata import collect
from boututils.showdata import showdata

# NOTE: Under construction

# All post processing functions called by bout_runners must accept the
# first argument from bout_runners (called 'folder' in
# __call_post_processing_function
def show_the_data(paths, order, **kwargs):
    """Function which makes an MMS test on the data

    Input
    paths   -   a list of all the paths for the MMS (sorted after
                ascending grid size)
    order   -   expected order of the test
    """


    print("Showing data from " + path)
    n = collect('n', xguards=False, yguards=False, path=path, info=False)

    # Pick out the key word arguments
    t = kwargs['t']
    x = kwargs['x']
    y = kwargs['y']
    z = kwargs['z']

    # Show the data
    showdata(n[t,x,y,z])
