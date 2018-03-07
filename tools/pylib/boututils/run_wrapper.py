"""Collection of functions which can be used to make a BOUT++ run"""

import os
import re
import subprocess
from os import getenv
try:
    # Python 2.4 onwards
    from subprocess import call, Popen, STDOUT, PIPE
    lib = "call"
except:
    # Use os.system (depreciated)
    from os import popen4, system
    lib = "system"
try:
  from builtins import str
except:
  pass


def getmpirun( default="mpirun -np" ):
  """
   getmpirun: return environment variable named MPIRUN, if it exists
              else return a default mpirun command
  """
  MPIRUN = getenv("MPIRUN")

  if MPIRUN is None:
    MPIRUN = default
    print("getmpirun: using the default " + str(default))

  return MPIRUN


def shell(command, pipe=False):
    """Run a shell command"""
    output = None
    status = 0
    if lib == "system":
        if pipe:
            handle = popen4(command)
            output = handle[1].read()
        else:
            status = system(command)
    else:
        if pipe:
            child = Popen(command, stderr=STDOUT, stdout=PIPE, shell=True)
            # This returns a b'string' which is casted to string in
            # python 2. However, as we want to use f.write() in our
            # runtest, we cast this to utf-8 here
            output = child.stdout.read().decode("utf-8")
            # Wait for the process to finish. Note that child.wait()
            # would have deadlocked the system as stdout is PIPEd, we
            # therefore use communicate, which in the end also waits for
            # the process to finish
            child.communicate()
            status = child.returncode
        else:
            status = call(command, shell=True)

    return status, output


def  determineNumberOfCPUs():
    """
    Number of virtual or physical CPUs on this system

    i.e.
    user/real as output by time(1) when called with an optimally scaling
    userspace-only program

    Taken from a post on stackoverflow:
    http://stackoverflow.com/questions/1006289/how-to-find-out-the-number-of-cpus-in-python
    """

    # Python 2.6+
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except (ImportError,NotImplementedError):
        pass

    # POSIX
    try:
        res = int(os.sysconf('SC_NPROCESSORS_ONLN'))

        if res > 0:
            return res
    except (AttributeError,ValueError):
        pass

    # Windows
    try:
        res = int(os.environ['NUMBER_OF_PROCESSORS'])

        if res > 0:
            return res
    except (KeyError, ValueError):
        pass

    # jython
    try:
        from java.lang import Runtime
        runtime = Runtime.getRuntime()
        res = runtime.availableProcessors()
        if res > 0:
            return res
    except ImportError:
        pass

    # BSD
    try:
        sysctl = subprocess.Popen(['sysctl', '-n', 'hw.ncpu'],
                                      stdout=subprocess.PIPE)
        scStdout = sysctl.communicate()[0]
        res = int(scStdout)

        if res > 0:
            return res
    except (OSError, ValueError):
        pass

    # Linux
    try:
        res = open('/proc/cpuinfo').read().count('processor\t:')

        if res > 0:
            return res
    except IOError:
        pass

    # Solaris
    try:
        pseudoDevices = os.listdir('/devices/pseudo/')
        expr = re.compile('^cpuid@[0-9]+$')

        res = 0
        for pd in pseudoDevices:
            if expr.match(pd) is not None:
                res += 1

        if res > 0:
            return res
    except OSError:
        pass

    # Other UNIXes (heuristic)
    try:
        try:
            dmesg = open('/var/run/dmesg.boot').read()
        except IOError:
            dmesgProcess = subprocess.Popen(['dmesg'], stdout=subprocess.PIPE)
            dmesg = dmesgProcess.communicate()[0]

        res = 0
        while '\ncpu' + str(res) + ':' in dmesg:
            res += 1

        if res > 0:
            return res
    except OSError:
        pass

    raise Exception('Can not determine number of CPUs on this system')


def launch(command, runcmd="mpirun -np", nproc=None, mthread=None, output=None, pipe=False, verbose=False):
    """Launch parallel MPI jobs

    status = launch(command, nproc, output=None)

    runcmd     Command for running parallel job; defaults to "mpirun -np"
    command    The command to run (string)
    nproc      Number of processors (integer)
    mthread      Number of omp threads (integer)
    output     Optional name of file for output
    """

    if nproc is None:
        # Determine number of CPUs on this machine
        nproc = determineNumberOfCPUs()

    cmd = runcmd + " " + str(nproc) + " " + command

    if output is not None:
        cmd = cmd + " > "+output

    if mthread is not None:
        cmd = "OMP_NUM_THREADS={j} ".format(j=mthread)+cmd
        
    if verbose == True:
         print(cmd)

    return shell(cmd, pipe=pipe)

def shell_safe(command,*args, **kwargs):
    """`Safe` version of shell.

    raises an RuntimeError exception if the command is not successfull.
    """
    s, out = shell(command,*args,**kwargs)
    if s:
        raise RuntimeError("Run failed with %d.\nCommand was:\n%s\n\n"
                           "Output was\n\n%s"%
                           (s,command,out))
    return s, out

def launch_safe(command,*args, **kwargs):
    """`Safe` version of launch.

    raises an RuntimeError exception if the command is not successfull.
    """
    s, out = launch(command,*args,**kwargs)
    if s:
        raise RuntimeError("Run failed with %d.\nCommand was:\n%s\n\n"
                           "Output was\n\n%s"%
                           (s,command,out))
    return s, out
