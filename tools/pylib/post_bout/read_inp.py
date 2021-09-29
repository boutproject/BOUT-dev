from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import str
from past.utils import old_div
from builtins import object
from read_grid import read_grid
from ordereddict import OrderedDict
import numpy as np
from boututils.file_import import file_import
from .read_cxx import *


def read_inp(path="", boutinp="BOUT.inp"):

    boutfile = path + "/" + boutinp
    boutinp = open(boutfile, "r").readlines()

    # start by stripping out all comments
    # look at the 1st character of all list elements
    # for now use a gross loop, vectorize later
    boutlist = []

    for i, val in enumerate(boutinp):
        if val[0] != "#" and val.isspace() == False:
            boutlist.append(val.split("#")[0])

    return boutlist


def parse_inp(boutlist):

    import re
    from ordereddict import OrderedDict

    if not boutlist:
        return 0

    # boutdic={} unordered standard dict
    boutdic = OrderedDict()

    # regex is messy see http://docs.python.org/howto/regex.html#regex-howto
    pattern = "\[\S*\]"  # pattern for a section start [All],[Ni], etc

    pattern = re.compile(pattern)

    boutdic["[main]"] = {}
    current = "[main]"

    for i, val in enumerate(boutlist):
        # print i,val
        result = pattern.match(val)
        # while the current value is not a new section name add everything to the current section

        if result is None:
            # print val
            key, value = val.split("=")
            value = value.replace('"', "")
            # print current, key,value

            boutdic[current][key.strip()] = value.strip()
        else:
            boutdic[result.group()] = {}
            current = result.group()

    return boutdic


def read_log(path=".", logname="status.log"):

    print("in read_log")
    import re
    from ordereddict import OrderedDict

    # logfile = path+'/'+logname
    logfile = logname
    print(logfile)
    logcase = open(logfile, "r").readlines()

    # start by stripping out all comments
    # look at the 1st character of all list elements
    # for now use a gross loop, vectorize later
    loglist = []

    for i, val in enumerate(logcase):
        if val[0] != "#" and val.isspace() == False:
            loglist.append(val.split("#")[0])

    if not loglist:
        return 0

    logdict = OrderedDict()
    logdict["runs"] = []
    # print len(loglist)
    print(loglist)
    # print loglist[len(loglist)-1] == 'last one\n'

    # last = loglist.pop().rstrip()

    # logdict['done'] = last == 'done'

    logdict["current"] = loglist.pop().rstrip()
    for i, val in enumerate(loglist):
        print(val)
        logdict["runs"].append(val.rstrip())

    logdict["runs"].append(logdict["current"])

    # print logdict
    return logdict


def metadata(inpfile="BOUT.inp", path=".", v=False):
    filepath = path + "/" + inpfile
    print(filepath)
    inp = read_inp(path=path, boutinp=inpfile)
    inp = parse_inp(inp)  # inp file
    print(path)
    outinfo = file_import(path + "/BOUT.dmp.0.nc")  # output data

    try:
        print(path)
        cxxinfo = no_comment_cxx(path=path, boutcxx="physics_code.cxx.ref")
        # evolved = get_evolved_cxx(cxxinfo)
        fieldkeys = get_evolved_cxx(cxxinfo)
        fieldkeys = ["[" + elem + "]" for elem in fieldkeys]
    except:
        print("cant find the cxx file")

    # gridoptions = {'grid':grid,'mesh':mesh}
    if "[mesh]" in list(inp.keys()):
        # IC = outinfo
        IC = read_grid(path + "/BOUT.dmp.0.nc")  # output data again
    elif "grid" in inp["[main]"]:
        gridname = inp["[main]"]["grid"]
        try:
            IC = read_grid(gridname)  # have to be an ansoulte file path for now
            print("IC: ", type(IC))
        # print IC.variables
        # print gridname
        except:
            # print gridname
            print("Fail to load the grid file")
    # print IC

    # print gridname
    # print len(IC)
    # print IC

    evolved = []
    collected = []
    ICscale = []

    # fieldkeys = ['[Ni]','[Te]','[Ti]','[Vi]','[rho]',
    #              '[Ajpar]','[Apar]','[vEBx]','[vEBy]','[vEBz]',
    #              '[jpar]','[phi]']

    # construct fieldkeys from cxx info
    # fieldkeys = ['['+x+']' for x in evolved]
    # fieldkeys = evolved

    # just look ahead and see what 3D fields have been output
    available = np.array([str(x) for x in outinfo])
    a = np.array([(len(outinfo[x].shape) == 4) for x in available])
    available = available[a]

    defaultIC = float(inp["[All]"].get("scale", 0.0))

    # print inp.keys()

    # figure out which fields are evolved
    print(fieldkeys)

    for section in list(inp.keys()):  # loop over section keys
        print("section: ", section)
        if section in fieldkeys:  # pick the relevant sections
            print(section)
            # print inp[section].get('evolve','True')
            # rint (inp[section].get('evolve','True')).lower().strip()
            if (
                inp[section].get("evolve", "True").lower().strip() == "true"
            ):  # and section[1:-1] in available :
                print("ok reading")
                evolved.append(section.strip("[]"))
                ICscale.append(float(inp[section].get("scale", defaultIC)))

        if inp[section].get("collect", "False").lower().strip() == "true":
            collected.append(section.strip("[]"))

    try:
        if inp["[physics]"].get("transport", "False").lower().strip() == "true":
            vEBstr = ["vEBx", "vEBy", "vEBz", "vEBrms"]
            [collected.append(item) for item in vEBstr]
    except:
        print("no [physics] key")

    meta = OrderedDict()

    class ValUnit(object):
        def __init__(self, value=0, units=""):
            self.u = units
            self.v = value

        def todict(self):
            return {"u": self.u, "v": self.v}

    # def decode_valunit(d):

    def ToFloat(metaString):
        try:
            return float(metaString)
        except ValueError:
            return metaString

    # meta['evolved'] = ValUnit(evolved,'')
    meta["evolved"] = evolved
    meta["collected"] = collected
    meta["IC"] = np.array(ICscale)
    d = {}

    print("evolved: ", evolved)

    # read meta data from .inp file, this is whre most metadata get written
    for section in list(inp.keys()):
        if ("evolve" not in inp[section]) and (
            "first" not in inp[section]
        ):  # hacky way to exclude some less relevant metadata
            for elem in list(inp[section].keys()):
                meta[elem] = ValUnit(ToFloat(inp[section][elem]))
                d[elem] = np.array(ToFloat(inp[section][elem]))

    # read in some values from the grid(IC) and scale them as needed
    norms = {
        "Ni0": ValUnit(1.0e14, "cm^-3"),
        "bmag": ValUnit(1.0e4, "gauss"),
        "Ni_x": ValUnit(1.0e14, "cm^-3"),
        "Te_x": ValUnit(1.0, "eV"),
        "Ti_x": ValUnit(1.0, "eV"),
        "Rxy": ValUnit(1, "m"),
        "Bxy": ValUnit(1.0e4, "gauss"),
        "Bpxy": ValUnit(1.0e4, "gauss"),
        "Btxy": ValUnit(1.0e4, "gauss"),
        "Zxy": ValUnit(1, "m"),
        "dlthe": ValUnit(1, "m"),
        "dx": ValUnit(1, "m"),
        "hthe0": ValUnit(1, "m"),
    }

    availkeys = np.array([str(x) for x in outinfo])
    tmp1 = np.array([x for x in availkeys])
    # b = np.array([x if x not in available for x in a])
    tmp2 = np.array([x for x in tmp1 if x not in available])
    static_fields = np.array([x for x in tmp2 if x in list(norms.keys())])
    # static_fields = tmp2

    # print availkeys
    # print meta.keys()
    # print IC.variables.keys()
    # print tmp1
    # print tmp2

    for elem in static_fields:
        print("elem: ", elem)
        meta[elem] = ValUnit(IC.variables[elem][:] * norms[elem].v, norms[elem].u)
        d[elem] = np.array(IC.variables[elem][:] * norms[elem].v)

    for elem in IC.variables:
        if elem not in meta:
            if elem in list(norms.keys()):
                meta[elem] = ValUnit(
                    IC.variables[elem][:] * norms[elem].v, norms[elem].u
                )
                d[elem] = np.array(IC.variables[elem][:] * norms[elem].v)
            else:
                meta[elem] = IC.variables[elem][:]
                d[elem] = IC.variables[elem][:]

    # print d.keys()

    # if case some values are missing
    default = {
        "bmag": 1,
        "Ni_x": 1,
        "NOUT": 100,
        "TIMESTEP": 1,
        "MZ": 32,
        "AA": 1,
        "Zeff": ValUnit(1, ""),
        "ZZ": 1,
        "zlowpass": 0.0,
        "transport": False,
    }
    diff = set(default.keys()).difference(set(d.keys()))

    for elem in diff:
        # print 'diff: ',elem
        meta[elem] = default[elem]
        d[elem] = np.array(default[elem])

    # print meta.keys()
    # print d.keys()

    # print meta['zlowpass']

    if meta["zlowpass"] != 0:
        print(meta["MZ"].v, meta["zlowpass"].v)
        meta["maxZ"] = int(np.floor(meta["MZ"].v * meta["zlowpass"].v))
    else:
        meta["maxZ"] = 5

    # meta['nx'] = nx
    # meta['ny']= ny
    meta["dt"] = meta["TIMESTEP"]

    # nx,ny  = d['Rxy'].shape

    # print meta['AA'].v

    meta["rho_s"] = ValUnit(
        1.02e2 * np.sqrt(d["AA"] * d["Te_x"]) / (d["ZZ"] * d["bmag"]), "cm"
    )  # ion gyrorad at T_e, in cm
    meta["rho_i"] = ValUnit(
        1.02e2 * np.sqrt(d["AA"] * d["Ti_x"]) / (d["ZZ"] * d["bmag"]), "cm"
    )
    meta["rho_e"] = ValUnit(2.38 * np.sqrt(d["Te_x"]) / (d["bmag"]), "cm")

    meta["fmei"] = ValUnit(1.0 / 1836.2 / d["AA"])

    meta["lambda_ei"] = 24.0 - np.log(old_div(np.sqrt(d["Ni_x"]), d["Te_x"]))
    meta["lambda_ii"] = 23.0 - np.log(
        d["ZZ"] ** 3 * np.sqrt(2.0 * d["Ni_x"]) / (d["Ti_x"] ** 1.5)
    )  #

    meta["wci"] = 1.0 * 9.58e3 * d["ZZ"] * d["bmag"] / d["AA"]  # ion gyrofrteq
    meta["wpi"] = (
        1.32e3 * d["ZZ"] * np.sqrt(old_div(d["Ni_x"], d["AA"]))
    )  # ion plasma freq

    meta["wce"] = 1.78e7 * d["bmag"]  # electron gyrofreq
    meta["wpe"] = 5.64e4 * np.sqrt(d["Ni_x"])  # electron plasma freq

    meta["v_the"] = 4.19e7 * np.sqrt(d["Te_x"])  # cm/s
    meta["v_thi"] = 9.79e5 * np.sqrt(old_div(d["Ti_x"], d["AA"]))  # cm/s
    meta["c_s"] = 9.79e5 * np.sqrt(5.0 / 3.0 * d["ZZ"] * d["Te_x"] / d["AA"])  #
    meta["v_A"] = 2.18e11 * np.sqrt(old_div(1.0, (d["AA"] * d["Ni_x"])))

    meta["nueix"] = 2.91e-6 * d["Ni_x"] * meta["lambda_ei"] / d["Te_x"] ** 1.5  #
    meta["nuiix"] = (
        4.78e-8
        * d["ZZ"] ** 4.0
        * d["Ni_x"]
        * meta["lambda_ii"]
        / d["Ti_x"] ** 1.5
        / np.sqrt(d["AA"])
    )  #
    meta["nu_hat"] = meta["Zeff"].v * meta["nueix"] / meta["wci"]

    meta["L_d"] = 7.43e2 * np.sqrt(old_div(d["Te_x"], d["Ni_x"]))
    meta["L_i_inrt"] = (
        2.28e7 * np.sqrt(old_div(d["AA"], d["Ni_x"])) / d["ZZ"]
    )  # ion inertial length in cm
    meta["L_e_inrt"] = 5.31e5 * np.sqrt(d["Ni_x"])  # elec inertial length in cm

    meta["Ve_x"] = 4.19e7 * d["Te_x"]

    meta["R0"] = old_div((d["Rxy"].max() + d["Rxy"].min()), 2.0)

    print(d["Rxy"].mean(1))
    print(d["ZMAX"])
    print(d["ZMIN"])
    meta["L_z"] = (
        1e2 * 2 * np.pi * d["Rxy"].mean(1) * (d["ZMAX"] - d["ZMIN"])
    )  # in cm toroidal range
    meta["dz"] = d["ZMAX"] - d["ZMIN"]

    # meta['lbNorm']=meta['L_z']*(d['Bpxy']/d['Bxy']).mean(1)     #-binormal coord range [cm]
    meta["lbNorm"] = meta["L_z"] * (old_div(d["Bxy"], d["Bpxy"])).mean(1)

    # meta['zPerp']=np.array(meta['lbNorm']).mean*np.array(range(d['MZ']))/(d['MZ']-1)
    # let's calculate some profile properties
    dx = np.gradient(d["Rxy"])[0]
    meta["L"] = (
        1.0
        * 1e2
        * dx
        * (meta["Ni0"].v)
        / np.gradient(meta["Ni0"].v)[0]
        / meta["rho_s"].v
    )

    meta["w_Ln"] = old_div(
        meta["c_s"], (np.min(abs(meta["L"])) * meta["wci"] * meta["rho_s"].v)
    )  # normed to wci

    AA = meta["AA"].v
    ZZ = d["ZZ"]
    Te_x = d["Te_x"]
    Ti_x = d["Ti_x"]
    fmei = meta["fmei"].v

    meta["lpar"] = (
        1e2 * ((old_div(d["Bxy"], d["Bpxy"])) * d["dlthe"]).sum(1) / meta["rho_s"].v
    )  # -[normed], average over flux surfaces, parallel length

    # yes dlthe is always the vertical displacement
    # dlthe = (hthe0*2 pi)/nz
    # meta['lpar']=1e2*(d['Bxy']/d['Bpxy']).mean(1)*d['dlthe'].mean(1) #function of x
    meta["sig_par"] = old_div(1.0, (fmei * 0.51 * meta["nu_hat"]))
    # meta['heat_nue'] = ((2*np.pi/meta['lpar'])**2)/(fmei*meta['nu_hat'])
    # kz_e = kz_i*(rho_e/rho_i)
    # kz_s = kz_i*(rho_s/rho_i)
    # kz_i = (TWOPI/L_z)*(indgen((*current_str).fft.nz+1))*rho_i

    # knorm = (TWOPI/lbNorm)*(indgen((*current_str).fft.nz+1))*rho_s

    # for now just translate
    for elem in meta:
        if type(meta[elem]).__name__ == "ValUnit":
            meta[elem] = {"u": meta[elem].u, "v": meta[elem].v}

    print("meta: ", type(meta))
    return meta

    # meta['DZ'] =inp['[main]']['ZMAX']#-b['[main]']['ZMIN']
    # AA = inp['[2fluid]']['AA']
    # Ni0 = IC.variables['Ni0'][:]*1.e14
    # bmag = IC.variables['bmag'][:]*1.e4 #to cgs
    # Ni_x = IC.variables['Ni_x'][:]*1.e14 # cm^-3
    # Te_x

    # rho_s = 1.02e2*sqrt(AA.v*Te_x.v)/ZZ.v/bmag.v
    # rho_i
    # rho_e


# for i,val in enumerate(boutlist):
