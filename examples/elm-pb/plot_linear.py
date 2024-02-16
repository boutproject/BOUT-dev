# Plots an analysis of the linear growth rate
#
# Input argument is the directory containing data files.
#
# Example:
#     $ python plot_linear.py data/
#

from boutdata import collect
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

if len(sys.argv) != 2:
    raise ValueError(f"Usage: {sys.argv[0]} path")

# Path to the data
path = sys.argv[1]

# Read pressure at last time point, to find peak
p = collect("P", path=path, tind=-1).squeeze()
prms = np.sqrt(np.mean(p**2, axis=-1))

pyprof = np.amax(prms, axis=0)
yind = np.argmax(pyprof)
pxprof = prms[:, yind]
xind = np.argmax(pxprof)
print(f"Peak amplitude at x = {xind}, y = {yind}")

# Read pressure time history at index of peak amplitude
p = collect("P", path=path, xind=xind, yind=yind).squeeze()

# p = p[:,:-1] # Remove point in Z

prms = np.sqrt(np.mean(p**2, axis=-1))

t = collect("t_array", path=path)
dt = t[1] - t[0]

gamma = np.gradient(np.log(prms)) / dt
growth_rate = np.mean(gamma[len(gamma) // 2 :])
growth_rate_std = np.std(gamma[len(gamma) // 2 :])

print(f"Mean growth rate: {growth_rate} +/- {growth_rate_std}")

fig, axs = plt.subplots(2, 2)

ax = axs[0, 0]
ax.plot(pyprof)
ax.set_xlabel("Y index")
ax.set_ylabel("RMS pressure")
ax.axvline(yind, linestyle="--", color="k")
ax.text(yind, 0.5 * pyprof[yind], f"y = {yind}")

ax = axs[0, 1]
ax.plot(pxprof)
ax.set_xlabel("X index")
ax.set_ylabel("RMS pressure")
ax.axvline(xind, linestyle="--", color="k")
ax.text(xind, 0.5 * pxprof[xind], f"x = {xind}")

ax = axs[1, 0]
ax.plot(t, prms)
ax.set_xlabel(r"Time [$\tau_A$]")
ax.set_ylabel("RMS pressure")
ax.set_yscale("log")
ax.plot(t, prms[-1] * np.exp(growth_rate * (t - t[-1])), "--k")

ax = axs[1, 1]
ax.plot(t, gamma)
ax.set_xlabel(r"Time [$\tau_A$]")
ax.set_ylabel("Growth rate [$\omega_A$]")
ax.axhline(growth_rate, linestyle="--", color="k")
ax.text(
    (t[-1] + t[0]) / 4,
    growth_rate * 1.2,
    rf"$\gamma = {growth_rate:.3f}\pm {growth_rate_std:.3f} \omega_A$",
)

plt.savefig(os.path.join(path, "linear_growth.pdf"))
plt.savefig(os.path.join(path, "linear_growth.png"))

plt.tight_layout()
plt.show()
