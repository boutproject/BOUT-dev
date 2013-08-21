# 
# Get environment variable
#

def getmpirun( default="mpirun -np" ):
  """
   getmpirun: return environment variable named MPIRUN, if it exists
              else return a default mpirun command
  """
  from os import getenv
  MPIRUN = getenv("MPIRUN")

  if MPIRUN == None:
    MPIRUN = default
    print "getmpirun: using the default ", default

  return MPIRUN   
