
import matplotlib.pyplot as plt

from boutdata.collect import collect

f = collect("f", path="data")
yup = collect("yup", path="data")
ydown = collect("ydown", path="data")

plt.plot(f[0,4,4,:], label="f")
plt.plot(yup[4,4,:], label="f.yup")
plt.plot(ydown[4,4,:], label="f.ydown")

plt.legend()

plt.savefig("plot_interp.pdf")
plt.show()
