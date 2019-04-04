#!/usr/bin/env python3
from boututils.datafile import DataFile
import matplotlib.pyplot as plt

path = "data"
monitors = [
    ["PROBES", "T_up"],
    ["PROBES", "n_up"],
    ["slow", "T"],
]

for pack in monitors:
    filename, data_name = pack
    t = DataFile(path+'/'+filename+'.dmp.0.nc').read('t_array')
    data = DataFile(path+'/'+filename+'.dmp.0.nc').read(data_name).flatten()
    plt.plot(t, data, label="{} {}".format(filename, data_name))

time = DataFile(path+'/BOUT.dmp.0.nc').read('t_array')
data = DataFile(path+'/BOUT.dmp.0.nc').read("T")[:, 2, 2, 0]

plt.plot(time, data, marker='+', label="BOUT++ T")

plt.xlabel("Time")
plt.legend()
plt.show()
