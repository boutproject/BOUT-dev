"""Collection of functions which can be used to make a BOUT++ run"""

from builtins import str
import os
import re
import subprocess

try:
    # Python 2.4 onwards
    from subprocess import call, Popen, STDOUT, PIPE
    lib = "call"
except ImportError:
    # FIXME: drop support for python < 2.4!
    # Use os.system (depreciated)
    from os import popen4, system
    lib = "system"

if os.name == "nt":
    # Default on Windows
    DEFAULT_MPIRUN = "mpiexec.exe -n"
else:
    DEFAULT_MPIRUN = "mpirun -np"


def getmpirun(default=DEFAULT_MPIRUN):
  """Return environment variable named MPIRUN, if it exists else return
   a default mpirun command

  Parameters
  ----------
  default : str, optional
      An mpirun command to return if ``MPIRUN`` is not set in the environment

  """
  MPIRUN = os.getenv("MPIRUN")

  if MPIRUN is None or MPIRUN == "":
    MPIRUN = default
    print("getmpirun: using the default " + str(default))

  return MPIRUN


def shell(command, pipe=False):
    """Run a shell command

    Parameters
    ----------
    command : list of str
        The command to run, split into (shell) words
    pipe : bool, optional
        Grab the output as text, else just run the command in the
        background

    Returns
    -------
    tuple : (int, str)
        The return code, and either command output if pipe=True else None
    """
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
            output = child.stdout.read().decode("utf-8", "ignore")
            # Wait for the process to finish. Note that child.wait()
            # would have deadlocked the system as stdout is PIPEd, we
            # therefore use communicate, which in the end also waits for
            # the process to finish
            child.communicate()
            status = child.returncode
        else:
            status = call(command, shell=True)

    return status, output


def determineNumberOfCPUs():
    """Number of virtual or physical CPUs on this system

    i.e. user/real as output by time(1) when called with an optimally
    scaling userspace-only program

    Taken from a post on stackoverflow:
    http://stackoverflow.com/questions/1006289/how-to-find-out-the-number-of-cpus-in-python

    Returns
    -------
    int
        The number of CPUs
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


def launch(command, runcmd=None, nproc=None, mthread=None,
           output=None, pipe=False, verbose=False):
    """Launch parallel MPI jobs

    >>> status = launch(command, nproc, output=None)

    Parameters
    ----------
    command : str
        The command to run
    runcmd : str, optional
        Command for running parallel job; defaults to what getmpirun() returns"
    nproc : int, optional
        Number of processors (default: all available processors)
    mthread : int, optional
        Number of omp threads (default: the value of the
        ``OMP_NUM_THREADS`` environment variable
    output : str, optional
        Name of file to save output to
    pipe : bool, optional
        If True, return the output of the command
    verbose : bool, optional
        Print the full command to be run before running it

    Returns
    -------
    tuple : (int, str)
        The return code, and either command output if pipe=True else None

    """

    if runcmd is None:
        runcmd = getmpirun()

    if nproc is None:
        # Determine number of CPUs on this machine
        nproc = determineNumberOfCPUs()

    cmd = runcmd + " " + str(nproc) + " " + command

    if output is not None:
        cmd = cmd + " > "+output

    if mthread is not None:
        if os.name == "nt":
            # We're on windows, so we have to do it a little different
            cmd = 'cmd /C "set OMP_NUM_THREADS={} && {}"'.format(mthread, cmd)
        else:
            cmd = "OMP_NUM_THREADS={} {}".format(mthread, cmd)
        
    if verbose == True:
         print(cmd)

    return shell(cmd, pipe=pipe)


def shell_safe(command, *args, **kwargs):
    """'Safe' version of shell.

    Raises a `RuntimeError` exception if the command is not
    successful

    Parameters
    ----------
    command : str
        The command to run
    *args, **kwargs
        Optional arguments passed to `shell`

    """
    s, out = shell(command,*args,**kwargs)
    if s:
        raise RuntimeError("Run failed with %d.\nCommand was:\n%s\n\n"
                           "Output was\n\n%s"%
                           (s,command,out))
    return s, out


def launch_safe(command, *args, **kwargs):
    """'Safe' version of launch.

    Raises an RuntimeError exception if the command is not successful

    Parameters
    ----------
    command : str
        The command to run
    *args, **kwargs
        Optional arguments passed to `shell`

    """
    s, out = launch(command,*args,**kwargs)
    if s:
        raise RuntimeError("Run failed with %d.\nCommand was:\n%s\n\n"
                           "Output was\n\n%s"%
                           (s,command,out))
    return s, out
