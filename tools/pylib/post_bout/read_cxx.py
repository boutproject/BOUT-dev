from builtins import range
from read_grid import read_grid
from ordereddict import OrderedDict
import numpy as np
import string
import re


def findlowpass(cxxstring):
    ##p1="lowPass\(.*\)"
    p1 = "lowPass\(.*\,(.*)\)"
    maxN = np.array(re.findall(p1, cxxstring))
    # print substrings

    # p2=("[0-9]")

    # maxN = np.array([re.findall(p2,elem) for elem in substrings]).flatten()

    if maxN.size == 0:
        return 100
    else:
        output = int(min(maxN))
        if output == 0:
            return 20
        else:
            return output


def no_comment_cxx(path=".", boutcxx="physics_code.cxx.ref"):
    # print 'no_comment'
    boutcxx = path + "/" + boutcxx
    # boutcxx = open(boutcxx,'r').readlines()
    f = open(boutcxx, "r")
    boutcxx = f.read()
    f.close()

    start = string.find(boutcxx, "/*")
    end = string.find(boutcxx, "*/") + 2

    s = boutcxx[0:start]
    for i in range(string.count(boutcxx, "/*")):
        start = string.find(boutcxx, "/*", end)
        s = s + boutcxx[end + 1 : start - 1]

        end = string.find(boutcxx, "*/", end) + 2

        s = s + boutcxx[end + 1 :]

    # pattern = "\n \s* \(//)* .* \n" #pattern for a section start [All],[Ni], etc
    pattern = "\n+.*;"  # everythin
    pattern = re.compile(pattern)
    result = re.findall(pattern, s)

    # print result

    nocomment = []
    for elem in result:
        # print elem
        elem = elem.lstrip()
        stop = elem.find("//")
        # print start,stop

        if stop > 0:
            nocomment.append(elem[0:stop])
        elif stop == -1:
            nocomment.append(elem)

    # result = pattern.match(val)
    # start = string.find(z,'\n  //')
    # end =string.find(boutcxx,'*/')+2
    # print nocomment

    return nocomment


def get_evolved_cxx(cxxfile=None):
    if cxxfile is None:
        cxxfile = no_comment_cxx()

    # s = cxxfile
    # section_0 = string.find(s,'int physics_run(BoutReal t)')
    # section_1 = string.find(s,'return',section_0)
    # s =  s[section_0:section_1]
    temp = []

    for x in cxxfile:
        i = x.find("bout_solve(")
        # print i,x
        if i != -1:
            comma_i = x[i::].find('"')
            comma_j = x[i::].rfind('"')
            # print x[i+comma_i:i+comma_j+1]
            temp.append(x[i + comma_i + 1 : i + comma_j])

    evolved = []
    [evolved.append(x) for x in set(temp)]
    return np.array(evolved)


def read_cxx(path=".", boutcxx="physics_code.cxx.ref", evolved=""):

    # print path, boutcxx
    boutcxx = path + "/" + boutcxx
    # boutcxx = open(boutcxx,'r').readlines()
    f = open(boutcxx, "r")
    boutcxx = f.read()
    f.close()

    # start by stripping out all comments
    # look at the 1st character of all list elements
    # for now use a gross loop, vectorize later

    start = string.find(boutcxx, "/*")
    end = string.find(boutcxx, "*/") + 2

    s = boutcxx[0:start]
    for i in range(string.count(boutcxx, "/*")):
        start = string.find(boutcxx, "/*", end)
        s = s + boutcxx[end + 1 : start - 1]

        end = string.find(boutcxx, "*/", end) + 2

        s = s + boutcxx[end + 1 :]

    section_0 = string.find(s, "int physics_run(BoutReal t)")
    section_1 = string.find(s, "return", section_0)
    s = s[section_0:section_1]

    tmp = open("./read_cxx.tmp", "w")
    tmp.write(s)
    tmp.close()
    tmp = open("./read_cxx.tmp", "r")

    cxxlist = ""

    for line in tmp:
        if line[0] != "//" and line.isspace() == False:
            cxxlist = cxxlist + line.split("//")[0]

    return cxxlist
