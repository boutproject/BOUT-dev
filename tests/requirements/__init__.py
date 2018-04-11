#!/usr/bin/env python3

import os
import re
from boututils.run_wrapper import shell

class Requirements(object):
    def __init__(self,path=None,verbose=False):
        if path == None:
            path=os.path.realpath(__file__)
            path=path.rsplit("/",1)[0]


        # Get list of files in subdirectory, excluding common temporaries
        requirements_list = [ x for x in os.listdir(path)
                              if not (("#" in x) or ("~" in x)
                                      or (x[0] == ".") or (x == "__init__.py")) ]

        if verbose:
            print("======= Requirement checks ========")

        self.requirements = {}
        for requirement in requirements_list:
            status,out = shell(os.path.join(path, requirement), pipe=True)
            self.requirements[requirement] = (status == 0)
            if verbose:
                print("{0} => {1}".format(requirement, (status == 0)))

        # Now have dictionary of requirements, true/false.

    def check(self,script):
        """Tests if the requirements of the given script (file name) are met.
        Returns the result (True/False) and the requirements expression.
        """

        with open(script, "r") as filein:
            contents = filein.read()

            # Find all lines starting with '#requires' or '#Requires:'
            match = re.findall("^\s*\#[Rr]equires:?(.*)", contents,re.MULTILINE)
            # Iterate over all expressions to evaluate
            for expr in match:
                    ret=eval(expr, self.requirements)
                    # Only return if this one causes the test to skip
                    if ret == False:
                        return ret, expr

        return True, "" # Default, if no requirements specified

    def add(self,requirement, value):
        """Add a requirement value to the dictionary"""
        self.requirements[requirement]=value

