
import matplotlib.pyplot as plt
from boutdata import collect
import numpy as np

run = True  # Run the simulations?

paths = ["div-integrate", "div"]

if run:
    from boututils.run_wrapper import shell_safe, launch_safe, getmpirun
    
    shell_safe("make > make.log")
    MPIRUN=getmpirun()
    
    for path in paths:
        launch_safe("./fci-wave -d "+path, runcmd=MPIRUN, nproc=nproc, pipe=False)

# Collect the results into a dictionary 
sum_n_B = {}

for path in paths:
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
    plt.title("Density n, "+path)
    plt.savefig(path+".pdf")
    plt.show()

# Make a plot comparing total sum density / B
    
plt.figure()
for path in paths:
    time, n_B = sum_n_B[path]
    plt.plot(time, n_B, label=path)
plt.legend()
plt.xlabel("Time")
plt.ylabel("Sum(n / B)")
plt.savefig("compare-density.pdf")

plt.show()

