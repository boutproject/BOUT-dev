
import matplotlib.pyplot as plt
from boutdata import collect
import numpy as np

run = True  # Run the simulations?
nproc = 2

# Note: Data from fci-wave-logn examples commented out.

data_noboundary = [
    ("div", "Model 1 (density, point interpolation)")
    ,("div-integrate", "Model 2 (density, area integration)")
    ,("logn", "Model 3 (log density, area integration)")
    #,("../fci-wave-logn/div-integrate", "Model 5 (velocity, log density, area integration)")
]

data_boundary = [
    ("boundary", "Model 2 (density, momentum)")
    ,("boundary-logn", "Model 3 (log density, momentum)")
    #,("../fci-wave-logn/boundary", "Model 5 (log density, velocity)")
    ]

# Change this to select no boundary or boundary cases
data = data_noboundary

if run:
    from boututils.run_wrapper import shell_safe, launch_safe
    shell_safe("make > make.log")

    for path,label in data:
        launch_safe("./fci-wave -d "+path, nproc=nproc, pipe=False)

# Collect the results into a dictionary 
sum_n_B = {}

for path,label in data:
    n = collect("n", path=path)
    Bxyz = collect("Bxyz", path=path)

    time = collect("t_array", path=path)
    
    nt, nx, ny, nz = n.shape
    
    n_B = np.ndarray(nt)
    for t in range(nt):
        n_B[t] = np.sum(n[t,:,:,:] / Bxyz)

    sum_n_B[path] = (time, n_B)

    # Plot the density at the final time
    
    plt.figure()
    plt.contourf(n[-1,:,0,:].T, 100)
    plt.colorbar()
    plt.xlabel("Major radius")
    plt.ylabel("Height")
    plt.title("Density n, "+label)
    plt.savefig(path+".pdf")
    plt.show()

# Make a plot comparing total sum density / B
    
plt.figure()
for path,label in data:
    time, n_B = sum_n_B[path]
    plt.plot(time, n_B, label=label)
plt.legend()
plt.xlabel("Time")
plt.ylabel("Sum(n / B)")
plt.savefig("compare-density.pdf")

plt.show()

