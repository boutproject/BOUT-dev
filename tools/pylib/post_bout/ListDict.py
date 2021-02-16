from __future__ import print_function
import sys
import os

try:
    boutpath = os.environ["BOUT_TOP"]
    pylibpath = boutpath + "/pylib"
    pbpath = pylibpath + "/post_bout"
    boutdatapath = pylibpath + "/boutdata"
    boututilpath = pylibpath + "/boututils"

    allpath = [boutpath, pylibpath, pbpath, boutdatapath, boututilpath]
    [sys.path.append(elem) for elem in allpath]
except:
    print("meh")

import numpy as np


def ListDictKey(input, key):
    # given a key and a list of dictionaries this method returns an ordered
    # list of requested key values
    output = []
    for x in input:
        try:
            # print x[key]
            output.append(x[key])
        except:
            print("Key not found")
            return 1

    return output


def ListDictFilt(input, key, valuelist):
    # given a key,value pair and a list of dictionaries this
    # method returns an ordered list of dictionaries where (dict(key)==value) = True
    # http://stackoverflow.com/questions/5762643/how-to-filter-list-of-dictionaries-with-matching-values-for-a-given-key
    try:
        x = copyf(input, key, valuelist)
        return x
    except:
        return []


def copyf(dictlist, key, valuelist):
    return [dictio for dictio in dictlist if dictio[key] in valuelist]


# def subset(obj):
#     #def __init__(self,alldb,key,valuelist,model=False):
#     selection = ListDictFilt(obj.mode_db,obj.key,valuelist)
#     if len(selection) !=0:
#         LinRes.__init__(obj,selection)
#         self.skey = key
#         if model==True:
#             self.model()
#         else:
#             LinRes.__init__(self,alldb)
#             if model==True:
#                 self.model()
