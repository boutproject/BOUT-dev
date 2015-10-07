# -*- coding: utf-8 -*-

# Adapted from http://stackoverflow.com/a/15860757/2043465

from sys import stdout

def update_progress(progress, barLength=10, ascii=False):
    """Displays or updates a console progress bar

    Accepts a float between 0 and 1. Any int will be converted to a float.
    A value under 0 represents a 'halt'.
    A value at 1 or bigger represents 100%

    Inputs
    ------
    progress  - Number between 0 and 1
    barLength - Length of the progress bar [10]
    ascii     - Use '#' as the progress indicator, otherwise use a Unicode
                character [False]
    """

    if ascii:
        cursor = '#'
    else:
        cursor = u'█'

    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = u"\rPercent: [{prog:-<{len}}] {perc:6.2f}% {stat}".format(
        len=barLength, prog=cursor*block, perc=progress*100, stat=status)
    stdout.write(text)
    stdout.flush()
