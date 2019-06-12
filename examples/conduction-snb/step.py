#
# SNB model example
#
# Uses a step in the temperature, intended for comparison to VFP results

length = 6e-4 # Domain length in m

qe = 1.602176634e-19

import numpy as np
import matplotlib.pyplot as plt

from boututils.run_wrapper import shell_safe, launch_safe
from boutdata.collect import collect

path = "step"

shell_safe("make > make.log")
# Run the case
s, out = launch_safe("./conduction-snb -d " + path, nproc=1, mthread=1, pipe=True)

Te = collect("Te", path=path).ravel()

ny = len(Te)
dy = length / ny

position = (np.arange(ny) + 0.5) * length / ny

# Read divergence of heat flux
div_q = collect("Div_Q", path=path).ravel()
div_q_SH = collect("Div_Q_SH", path=path).ravel()

# Integrate the divergence of flux to get heat flux W/m^2
q = np.cumsum(div_q) * qe * dy
q_SH = np.cumsum(div_q_SH) * qe * dy

# Subtract the minimum value of each heat flux
q -= np.amin(q)
q_SH -= np.amin(q_SH)

fig, ax1 = plt.subplots()

color='tab:red'
ax1.plot(position, Te * 1e-3, color=color)
ax1.set_xlabel("position [m]")
ax1.set_ylabel("Electron temperature [keV]", color=color)
ax1.set_ylim(0,1)
ax1.tick_params(axis='y', colors=color)

ax2 = ax1.twinx()
ax2.plot(position, q_SH * 1e-4, '--k', label="Spitzer-Harm")
ax2.plot(position, q * 1e-4, label="SNB")
ax2.set_ylabel("Heat flux W/cm^2")
ax2.set_ylim(bottom=0.0)

plt.legend()
fig.tight_layout()

plt.xlim(0,3.5e-4)

plt.savefig("snb-step.png")
plt.show()

