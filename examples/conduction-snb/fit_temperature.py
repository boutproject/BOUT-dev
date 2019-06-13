import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

Te_ref = np.loadtxt("temperature.csv", delimiter=",")
Te_ref[:,0] *= 1e-4 # Convert X axis to m

def te_function(ypos, mid, wwid, w0, w1, w2, Tmax, Tmin, clip=False):
    width = w0 + ((ypos - mid)*w1 + (ypos - mid)**2 * w2) * np.exp(-((ypos - mid)/wwid)**2)

    if clip:
        width = np.clip(width, 1e-10, None)
    
    return Tmax - 0.5 * (1 + np.tanh((ypos - mid)/width)) * (Tmax - Tmin)

popt, pcov = optimize.curve_fit(te_function, Te_ref[:,0], Te_ref[:,1],
                                p0 = [2.2e-4, 1e-4, 1e-4, 0.0, 0.0, 0.960, 0.190])

print(popt)

xfit = np.linspace(Te_ref[0,0], Te_ref[-1,0], 100)

plt.plot(xfit, te_function(xfit, *popt, clip=True), '-k')
plt.plot(Te_ref[:,0], Te_ref[:,1], 'or')
plt.show()
