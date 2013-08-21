# Launch a parallel job
# 
# If run on a known machine, will generate a run script and submit a job for
# execution. By default just uses mpirun

from boututils import shell, determineNumberOfCPUs

def launch(command, runcmd="mpirun -np", nproc=None, output=None, pipe=False, verbose=False):
    """Launch parallel MPI jobs
    
    status = launch(command, nproc, output=None)

    runcmd     Command for running parallel job; defaults to "mpirun -np"    
    command    The command to run (string)
    nproc      Number of processors (integer)
    output     Optional name of file for output
    """
    
    if nproc == None:
        # Determine number of CPUs on this machine
        nproc = determineNumberOfCPUs()

    cmd = runcmd + " " + str(nproc) + " " + command
    
    if output != None:
        cmd = cmd + " > "+output

    if verbose == True:
         print cmd    

    return shell(cmd, pipe=pipe)
