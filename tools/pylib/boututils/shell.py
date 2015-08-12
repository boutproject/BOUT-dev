# Run a shell command
#
#

try:
    # Python 2.4 onwards
    from subprocess import call, Popen, STDOUT, PIPE
    lib = "call"
except:
    # Use os.system (depreciated)
    from os import popen4, system
    lib = "system"


def shell(command, pipe=False):
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
