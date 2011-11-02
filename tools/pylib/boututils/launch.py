# Launch a parallel job
# 
# If run on a known machine, will generate a run script and submit a job for
# execution. By default just uses mpirun

from boututils import shell, determineNumberOfCPUs

def launch(command, nproc=None, output=None, pipe=False):
    """Launch parallel MPI jobs
    
    status = launch(command, nproc, output=None)
    
    command    The command to run (string)
    nproc      Number of processors (integer)
    output     Optional name of file for output
    """
    
    if nproc == None:
        # Determine number of CPUs on this machine
        nproc = determineNumberOfCPUs()

    cmd = "mpirun -np " + str(nproc) + " " + command
    
    if output != None:
        cmd = cmd + " > "+output
    
    return shell(cmd, pipe=pipe)

