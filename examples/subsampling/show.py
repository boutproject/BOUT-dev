from boututils.datafile import  DataFile
import matplotlib.pyplot as plt
path="data"
todo=[["slow","T"],
      ["PROBES","T_up"],
      ["PROBES","n_up"],
]
for pack in todo:
    filename,data = pack
    t=DataFile(path+'/'+filename+'.dmp.0.nc').read('t_array').flatten()
    data=DataFile(path+'/'+filename+'.dmp.0.nc').read(data).flatten()
    plt.plot(t,data)
plt.show()
