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
            output = child.stdout.read()
            status = child.returncode
        else:
            status = call(command, shell=True)
        
    return status, output
