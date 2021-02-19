#!/usr/bin/env python3

import os
import re
from boututils.run_wrapper import shell


class Requirements(object):
    """Checks if test cases meet their specified requirements. Test
    scripts can contain a line starting with #requires or #Requires,
    optionally followed by a colon, then a Python expression which
    combines boolean values.

    Examples:

    #requires netcdf

    #requires not make

    #requires not (travis and netcdf)

    The individual requirements (netcdf, make, travis, etc.)
    are gathered from the bout-config scipt, or from executable
    scripts in the "requirements" directory (where this code is).

    """

    def __init__(self, path=None, verbose=False):
        selflocation = os.path.realpath(__file__)
        selflocation = selflocation.rsplit("/", 1)[0]

        if path is None:
            path = selflocation

        # Get list of files in subdirectory, excluding common temporaries,
        # hidden files, and python .py and .pyc files
        requirements_list = [
            x
            for x in os.listdir(path)
            if not (("#" in x) or ("~" in x) or (x[0] == ".") or (".py" in x))
        ]

        if verbose:
            print("======= Requirement checks ========")

        self._verbose = verbose
        self._requirements = {}
        for requirement in requirements_list:
            status, out = shell(os.path.join(path, requirement), pipe=True)
            self.add(requirement, (status == 0))

        with open(selflocation + "/../../bin/bout-config") as configfile:
            config = configfile.read()
            matches = re.findall("^has_(.*)=(.*)", config, re.MULTILINE)
            for match in matches:
                key = match[0]
                value = match[1]
                yesno = {'"yes"': True, '"no"': False}
                try:
                    value = yesno[value]
                except KeyError:
                    print(
                        "Error parsing " + match + ' - %s is not "yes"/"no"' % match[1]
                    )
                else:
                    self.add(key, value)

        # Now have dictionary of requirements, true/false.

    def check(self, script):
        """Tests if the requirements of the given script (file name) are met.
        Returns the result (True/False) and the requirements expression.
        """

        with open(script, "r", encoding="utf-8") as filein:
            contents = filein.read()

            # Find all lines starting with '#requires' or '#Requires:'
            match = re.findall("^\s*\#\s?[Rr]equires:?(.*)", contents, re.MULTILINE)
            # Iterate over all expressions to evaluate
            for expr in match:
                try:
                    ret = eval(expr, self._requirements)
                except:
                    print("Parsing '%s' from file %s failed." % (expr, script))
                    print("A list of known variable follows:")
                    for key, value in self._requirements.items():
                        if "_" not in key:
                            print("  %s => %s" % (key, value))
                    raise
                else:
                    # Only return if this one causes the test to skip
                    if ret is False:
                        return ret, expr

        return True, ""  # Default, if no requirements specified

    def add(self, requirement, value, override=False):
        """Add a requirement value to the dictionary. If the value is
        already defined with another value, it will raise an
        exception, unless override is set to True."""
        if requirement in self._requirements:
            old = self._requirements[requirement]
            if old != value and override is False:
                raise RuntimeError(
                    "Overwriting %s with %s in key %s" % (value, old, requirement)
                )
        self._requirements[requirement] = value
        if self._verbose:
            print("{0} => {1}".format(requirement, value))

    def has(self, requirement):
        """Check if a requrirement is available"""
        return self._requirements[requirement]
