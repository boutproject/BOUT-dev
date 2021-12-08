# note - these commands are only run by default in interactive mode
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import filter
from builtins import range
from past.utils import old_div
from builtins import object
import os
import sys

try:
    boutpath = os.environ["BOUT_TOP"]
    pylibpath = boutpath + "tools/pylib"
    pbpath = pylibpath + "/post_bout"
    boutdatapath = pylibpath + "/boutdata"
    boututilpath = pylibpath + "/boututils"

    allpath = [boutpath, pylibpath, pbpath, boutdatapath, boututilpath]
    [sys.path.append(elem) for elem in allpath]

except:
    print("unable to append needed .py files")

sys.path.append("/usr/local/pylib")

import post_bout as post_bout
from .ListDict import ListDictKey, ListDictFilt
from .read_inp import parse_inp, read_inp, read_log
from .basic_info import weighted_avg_and_std
from .read_cxx import read_cxx, findlowpass


import os
import numpy as np
import pickle
import subprocess


def corral(
    cached=True, refresh=False, debug=False, IConly=1, logname="status.log", skew=False
):

    print("in corral")
    log = read_log(logname=logname)
    # done = log['done']
    runs = log["runs"]  # a list of all directories, we need this,
    # only need 'runs' if the simulation is done

    current = log["current"]  # always return the last data_dir

    print(log)
    print("current:", current)

    if refresh == True:

        for i, path in enumerate(runs):
            print(i, path)
            a = post_bout.save(path=path, IConly=IConly)  # re post-process a run

    elif (
        cached == False
    ):  # if all the ind. simulation pkl files are in place skip this part

        a = post_bout.save(path=current)  # save to current dir
        # here is really where you shoudl write to status.log
        # write_log('status.log',
        cached = True

    # if done:

    all_ave = []
    all_modes = []
    print("last_one: ")
    for i, val in enumerate(runs):
        print(val)
        mode_db, ave_db = post_bout.read(path=val)
        # alldata.append(array)
        all_modes.append(mode_db)
        all_ave.append(ave_db)

        # build the end database

        # remove the read in pickle

    # return alldb
    def islist(input):
        return isinstance(input, list)

    all_modes = list(filter(islist, all_modes))
    # all_ave = filter(islist,all_ave)

    # alldb = sum(alldb,[])
    # alldata = np.array(alldata)
    all_modes = sum(all_modes, [])

    nt = []
    for mode in all_modes:
        nt.append(len(mode["amp"]))
        nt = [max(nt)]

    nt = nt[0]
    t = list(range(nt))
    i = 0

    if debug:
        return all_modes, all_ave
    else:
        return LinRes(all_modes)


class LinRes(object):
    def __init__(self, all_modes):

        self.mode_db = all_modes
        self.db = all_modes
        # self.ave_db = all_ave

        alldb = self.db
        # self.modekeys = data[0]['fields']['Ni']['modes'][0].keys()
        # print len(alldb)

        self.meta = np.array(ListDictKey(self.db, "meta"))[0]

        self.keys = list((self.mode_db)[0].keys())
        # self.avekeys = data[0]['fields']['Ni']['ave'].keys()

        # self.nrun = len(alldb) #number of runs

        self.path = np.array(ListDictKey(self.db, "path"))
        self.cxx = []
        self.maxN = []

        self.ave = np.array(ListDictKey(alldb, "ave"))

        # [self.cxx.append(read_cxx(path=elem,boutcxx='2fluid.cxx.ref')) for elem in self.path]
        # [self.maxN.append(findlowpass(elem)) for elem in self.cxx]
        [
            self.cxx.append(read_cxx(path=elem, boutcxx="physics_code.cxx.ref"))
            for elem in self.path
        ]

        self.maxZ = np.array(ListDictKey(alldb, "maxZ"))
        self.maxN = self.maxZ
        # self.maxN = findlowpass(self.cxx) #low pass filt from .cxx

        self.nx = np.array(ListDictKey(alldb, "nx"))[0]
        self.ny = np.array(ListDictKey(alldb, "ny"))[0]
        self.nz = np.array(ListDictKey(alldb, "nz"))[0]

        # self.nt = int(data[0]['meta']['NOUT']['v']+1)
        self.Rxy = np.array(ListDictKey(alldb, "Rxy"))
        self.Rxynorm = np.array(ListDictKey(alldb, "Rxynorm"))
        self.nt = np.array(ListDictKey(alldb, "nt"))

        self.dt = np.array(ListDictKey(alldb, "dt"))
        self.nfields = np.array(ListDictKey(alldb, "nfields"))

        self.field = np.array(ListDictKey(alldb, "field"))

        self.k = np.array(ListDictKey(alldb, "k"))
        self.k_r = np.array(ListDictKey(alldb, "k_r"))

        self.mn = np.array(ListDictKey(alldb, "mn"))

        # return ListDictKey(alldb,'phase')

        # self.phase = np.array(ListDictKey(alldb,'phase'))
        self.phase = ListDictKey(alldb, "phase")

        # self.amp= np.array(ListDictKey(alldb,'amp'))
        # self.amp_n=np.array(ListDictKey(alldb,'amp_n'))
        # self.dc= []
        # self.freq = np.array(ListDictKey(alldb,'k'))
        # self.gamma = np.array(ListDictKey(alldb,'gamma'))

        self.amp = ListDictKey(alldb, "amp")
        self.amp_n = ListDictKey(alldb, "amp_n")
        self.dc = []
        # self.freq = np.array(ListDictKey(alldb,'k'))
        self.gamma = np.array(ListDictKey(alldb, "gamma"))
        self.gamma_i = np.array(ListDictKey(alldb, "gamma_i"))

        self.freq = np.array(ListDictKey(alldb, "freq"))

        self.IC = np.array(ListDictKey(alldb, "IC"))
        self.dz = np.array(ListDictKey(alldb, "dz"))
        self.meta["dz"] = np.array(list(set(self.dz).union()))

        self.nmodes = self.dz.size

        self.MN = np.array(ListDictKey(alldb, "MN"))
        # self.MN = np.float32(self.mn)
        # self.MN[:,1] = self.mn[:,1]/self.dz
        self.nrun = len(set(self.path).union())
        self.L = np.array(ListDictKey(alldb, "L"))
        # self.C_s =
        self.modeid = np.array(ListDictKey(alldb, "modeid"))

        self.trans = np.array(ListDictKey(alldb, "transform"))

        if np.any(self.trans):
            self.phase_r = ListDictKey(alldb, "phase_r")
            self.gamma_r = np.array(ListDictKey(alldb, "gamma_r"))
            self.amp_r = ListDictKey(alldb, "amp_r")
            self.freq_r = np.array(ListDictKey(alldb, "freq_r"))

        # try:
        #     self.model(haswak=False) #
        # except:
        #     self.M = 0

        # try:
        try:  # analytic model based on simple matrix
            self.models = []
            # self.models.append(_model(self)) #create a list to contain models
            self.models.append(
                _model(self, haswak=True, name="haswak")
            )  # another model
        # self.models.append(_model(self,haswak=True,name='haswak_0',m=0))
        # for Ln in range(10):
        #     Lval = 10**((.2*Ln -1)/10)
        #     #Lval = 10**(Ln-1)
        #     #Lval =
        #     print Lval
        #     self.models.append(_model(self,varL=True,name='varL'+str(Lval),Lval=Lval,haswak=True))

        except:
            self.M = 0

        try:  # analytic models based on user defined complex omega
            self.ref = []
            self.ref.append(_ref(self))
        # self.ref.append(_ref(self,haswak=False,name='drift'))
        # demand a complex omega to compare
        # self.ref.append(_ref(self,haswas=True,name='haswak'))
        except:
            self.ref = 0

    # self.models.append(_model(self,haswak2=True,name='haswak2'))

    # except:
    #     print 'FAIL'

    def _amp(self, tind, xind):
        # first select modes that actually have valid (tind,xind)
        # indecies
        # s = subset(self.db,'modeid',modelist)
        return np.array([self.amp[i][tind, xind] for i in range(self.nmodes)])

    # def model(self,field='Ni',plot=False,haswak=False):

    #    #enrich the object
    #    allk = self.k_r[:,1,self.nx/2] #one location for now
    #    allkpar = self.k_r[:,0,self.nx/2] #one location for now

    #    self.M = []
    #    self.eigsys = []
    #    self.gammaA = []
    #    self.omegaA = []
    #    self.eigvec = []
    #    self.gammamax = []
    #    self.omegamax = []

    #    #allk = np.arange(0.1,100.0,.1)
    #    #allk=  np.sort(list(set(allk).union()))

    #    for i,k in enumerate(allk):
    #       #print i
    #       #M =np.matrix(np.random.rand(3,3),dtype=complex)
    #       M = np.zeros([3,3],dtype=complex)
    #       M[0,0] = 0
    #       M[0,1] = k/(self.L[i,self.nx/2,self.ny/2])
    #       M[1,0] = (2*np.pi/self.meta['lpar'][self.nx/2])**2 * self.meta['sig_par'][0]*complex(0,k**-2)
    #       M[1,1]= -(2*np.pi/self.meta['lpar'][self.nx/2])**2 * self.meta['sig_par'][0]*complex(0,k**-2)

    #       if haswak:
    #           M[0,0] = M[0,0] + M[1,1]*complex(0,k**2)
    #           M[0,1] = M[0,1] + M[1,0]*complex(0,k**2)

    #       #if rho_conv:

    #       #M[1,0] = (allkpar[i])**2 * self.meta['sig_par'][0]*complex(0,k**-2)
    #       #M[1,1]= -(allkpar[i])**2 * self.meta['sig_par'][0]*complex(0,k**-2)

    #       eigsys= np.linalg.eig(M)
    #       gamma = (eigsys)[0].imag
    #       omega =(eigsys)[0].real
    #       eigvec = eigsys[1]

    #       self.M.append(M)
    #       self.eigsys.append(eigsys)
    #       self.gammaA.append(gamma)
    #       self.gammamax.append(max(gamma))
    #       where = ((gamma == gamma.max()) & (omega != 0))
    #       self.omegamax.append(omega[where[0]])
    #       self.eigvec.append(eigvec)
    #       self.omegaA.append(omega)

    class __model__(object):
        def __init__(self):
            self.M = 0


class subset(LinRes):
    def __init__(self, alldb, key, valuelist, model=False):
        selection = ListDictFilt(alldb, key, valuelist)
        if len(selection) != 0:
            LinRes.__init__(self, selection)
            self.skey = key
            if model == True:
                self.model()
        else:
            LinRes.__init__(self, alldb)
            if model == True:
                self.model()


# class subset(originClass):
#     def __init__(self,alldb,key,valuelist,model=False):
#       selection = ListDictFilt(alldb,key,valuelist)
#       if len(selection) !=0:
#          originClass.__init__(self,selection,input_obj.ave_db)
#          self.skey = key
#          if model==True:
#             self.model()
#       else:
#          origin.__init__(self,alldb)
#          if model==True:
#             self.model()

# not sure if this is the best way . . .


# class subset(object):
#     def __init__(self,input_obj,key,valuelist,model=False):
#         selection = ListDictFilt(input_obj.mode_db,key,valuelist)
#         if len(selection) !=0:
#             import copy
#             self = copy.copy(input_obj)
#             self.__init__(selection,input_obj.ave_db)
#             self.skey = key
#             if model==True:
#                 self.model()
#         else:
#             self = input_obj
#             if model==True:
#                 self.model()


class _ref(object):  # NOT a derived obj, just takes one as a var
    def __init__(self, input_obj, name="haswak", haswak=True):
        allk = input_obj.k_r[:, 1, old_div(input_obj.nx, 2)]  # one location for now
        allkpar = input_obj.k_r[:, 0, old_div(input_obj.nx, 2)]  # one location for now
        self.name = name

        self.gamma = []
        self.omega = []
        self.soln = {}
        self.soln["gamma"] = []
        self.soln["freq"] = []

        for i, k in enumerate(allk):
            omega_star = old_div(
                -(k),
                (input_obj.L[i, old_div(input_obj.nx, 2), old_div(input_obj.ny, 2)]),
            )

            nu = (
                2 * np.pi / input_obj.meta["lpar"][old_div(input_obj.nx, 2)]
            ) ** 2 * input_obj.meta["sig_par"][0]

            if haswak:
                omega = old_div(-omega_star, (1 + (k) ** 2))
                gamma = old_div(((k ** 2) * omega_star ** 2), (nu * (1 + k ** 2) ** 3))
            else:
                # omega = -np.sqrt(nu*omega_star)/(np.sqrt(2)*k) + nu**(3/2)/(8*np.sqrt(2*omega_star)*k**3)
                # gamma = np.sqrt(nu*omega_star)/(np.sqrt(2)*k) - nu/(2* k**2) + nu**(3/2)/(8*np.sqrt(2*omega_star)*k**3)
                omega = -omega_star + old_div((2 * k ** 4 * omega_star ** 3), nu ** 2)
                gamma = old_div((k * omega_star) ** 2, nu) - (
                    5 * (k ** 6 * omega_star * 4 / nu ** 3)
                )
            self.gamma.append(gamma)
            self.omega.append(omega)

        self.soln["freq"] = np.transpose(np.array(self.omega))
        self.soln["gamma"] = np.transpose(np.array(self.gamma))


class _model(object):  # NOT a derived class,but one that takes a class as input
    def __init__(
        self,
        input_obj,
        name="drift",
        haswak=False,
        rho_conv=False,
        haswak2=False,
        varL=False,
        Lval=1.0,
        m=1,
    ):
        allk = input_obj.k_r[:, 1, old_div(input_obj.nx, 2)]  # one location for now
        allkpar = input_obj.k_r[:, 0, old_div(input_obj.nx, 2)]  # one location for now

        # numerical value to compare against

        numgam = input_obj.gamma[:, 0, old_div(input_obj.nx, 2)]
        numfreq = input_obj.freq[:, 0, old_div(input_obj.nx, 2)]
        self.name = name

        self.M = []
        self.eigsys = []
        self.gammaA = []
        self.omegaA = []
        self.eigvec = []
        self.gammamax = []
        self.omegamax = []
        self.k = []
        self.m = m

        self.soln = {}
        self.soln["freq"] = []
        self.soln["gamma"] = []
        self.soln["gammamax"] = []
        self.soln["freqmax"] = []

        self.chi = {}
        self.chi["freq"] = []
        self.chi["gamma"] = []

        for i, k in enumerate(allk):
            # print i
            # M =np.matrix(np.random.rand(3,3),dtype=complex)
            M = np.zeros([4, 4], dtype=complex)
            M[0, 0] = 0
            # k = k/np.sqrt(10)
            # L = (input_obj.L)*np.sqrt(10)

            if k == 0:
                k = 1e-5

            # print k {n,phi,v,ajpar}
            M[0, 1] = old_div(
                k, (input_obj.L[i, old_div(input_obj.nx, 2), old_div(input_obj.ny, 2)])
            )
            M[1, 0] = (
                (2 * m * np.pi / input_obj.meta["lpar"][old_div(input_obj.nx, 2)]) ** 2
                * input_obj.meta["sig_par"][0]
                * complex(0, (k) ** -2)
            )
            M[1, 1] = (
                -(
                    (2 * m * np.pi / input_obj.meta["lpar"][old_div(input_obj.nx, 2)])
                    ** 2
                )
                * input_obj.meta["sig_par"][0]
                * complex(0, (k) ** -2)
            )

            # parallel dynamics
            # M[2,2] = k/(input_obj.L[i,input_obj.nx/2,input_obj.ny/2])
            # M[2,0] = -(2*m*np.pi/input_obj.meta['lpar'][input_obj.nx/2])
            # M[0,2] = -(2*m*np.pi/input_obj.meta['lpar'][input_obj.nx/2])

            # M[1,0] = (2*m*np.pi/input_obj.meta['lpar'][input_obj.nx/2])**2 * input_obj.meta['sig_par'][0]
            # M[1,1]= -(2*m*np.pi/input_obj.meta['lpar'][input_obj.nx/2])**2 * input_obj.meta['sig_par'][0]

            # ajpar dynamics - effectively parallel electron dynamics instead of

            if haswak:
                M[0, 0] = (
                    -(
                        (
                            2
                            * m
                            * np.pi
                            / input_obj.meta["lpar"][old_div(input_obj.nx, 2)]
                        )
                        ** 2
                    )
                    * input_obj.meta["sig_par"][0]
                    * complex(0, 1)
                )
                M[0, 1] = (
                    2 * m * np.pi / input_obj.meta["lpar"][old_div(input_obj.nx, 2)]
                ) ** 2 * input_obj.meta["sig_par"][0] * complex(0, 1) + M[0, 1]

            if varL:
                M[0, 1] = Lval * M[0, 1]

            if rho_conv:  # not used
                M[1, 0] = (
                    (allkpar[i]) ** 2
                    * input_obj.meta["sig_par"][0]
                    * complex(0, (k) ** -2)
                )
                M[1, 1] = (
                    -((allkpar[i]) ** 2)
                    * input_obj.meta["sig_par"][0]
                    * complex(0, (k) ** -2)
                )

            eigsys = np.linalg.eig(M)
            gamma = (eigsys)[0].imag
            omega = (eigsys)[0].real
            eigvec = eigsys[1]
            self.k.append(k)

            self.M.append(M)
            self.eigsys.append(eigsys)

            self.gammaA.append(gamma)
            self.soln["gamma"].append(gamma)

            self.gammamax.append(max(gamma))
            self.soln["gammamax"].append(max(gamma))

            where = (gamma == gamma.max()) & (omega != 0)
            # if len(where) > 1:
            #     where = where[0]
            self.omegamax.append(omega[where])
            self.soln["freqmax"].append(omega[where])

            # print k,gamma,where,M,omega
            chigam = old_div(((numgam - max(gamma)) ** 2), max(gamma))
            chifreq = old_div(((numfreq - omega[where]) ** 2), omega[where])

            self.eigvec.append(eigvec)
            self.omegaA.append(omega)
            self.soln["freq"].append(omega)

            self.chi["freq"].append(chifreq[i])

            self.chi["gamma"].append(chigam[i])

        self.dim = M.shape[0]
        self.soln["freq"] = np.transpose(np.array(self.soln["freq"]))
        self.soln["gamma"] = np.transpose(np.array(self.soln["gamma"]))
        self.chi["freq"] = np.transpose(np.array(self.chi["freq"]))
        self.chi["gamma"] = np.transpose(np.array(self.chi["gamma"]))

        # self.soln = {}
        # self.soln['freq'] = self.omegaA
        # self.soln['gamma'] = self.gammaA
