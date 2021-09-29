# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 15:55:04 2016

@author: yolen
"""
import boutdata
import matplotlib.pyplot as plt

path = "."
E_n_ = boutdata.collect("E_N", xguards=False, yguards=False, path=path)
# S_n_ = boutdata.collect('S_N',xguards=False, yguards = False, path = path)
n = boutdata.collect("N", xguards=False, path=path)

print("shape E_n = {}".format(E_n_.shape))
# n_source_ = boutdata.collect('n_source_',xguards=False)

plt.subplot(2, 2, 1)
plt.plot(E_n_[1, 0, 0, :], "ro")
plt.title("error")
# plt.colorbar()
plt.subplot(2, 2, 2)
plt.plot(n[1, 0, 0, :], "ro")
plt.title("solution")
# plt.colorbar()
# plt.subplot(2,2,3)
# plt.plot(S_n_[1,1,0,:],'ro')
# plt.title('source')
# plt.colorbar()

plt.show()
