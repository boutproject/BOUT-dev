#!/usr/bin/env python3

import os
import re
from boututils.run_wrapper import shell

class Requirements(object):
    def __init__(self,path=None,verbose=False):
        selflocation=os.path.realpath(__file__)
        selflocation=selflocation.rsplit("/",1)[0]

        if path == None:
            path=selflocation

        # Get list of files in subdirectory, excluding common temporaries
        requirements_list = [ x for x in os.listdir(path)
                              if not (("#" in x) or ("~" in x)
                                      or (x[0] == ".") or (x == "__init__.py")) ]

        if verbose:
            print("======= Requirement checks ========")

        self._requirements = {}
        for requirement in requirements_list:
            status,out = shell(os.path.join(path, requirement), pipe=True)
            self._requirements[requirement] = (status == 0)
            if verbose:
                print("{0} => {1}".format(requirement, (status == 0)))

        with open(selflocation+"/../../bin/bout-config") as configfile:
            config=configfile.read()
            matches=re.findall("^has_(.*)=(.*)",config,re.MULTILINE)
            for match in matches:
                key=match[0]
                value=match[1]
                yesno={'"yes"':True, '"no"': False}
                try:
                    value=yesno[value]
                except KeyError:
                    print("Error parsing "+match+" - %s is not \"yes\"/\"no\""%match[1])
                else:
                    self._requirements[key]=value

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
                    ret=eval(expr, self._requirements)
                    # Only return if this one causes the test to skip
                    if ret == False:
                        return ret, expr

        return True, "" # Default, if no requirements specified

    def add(self,requirement, value):
        """Add a requirement value to the dictionary"""
        self._requirements[requirement]=value


    def has(self,requirement):
        """Check if a requrirement is available"""
        return self._requirements[requirement]
