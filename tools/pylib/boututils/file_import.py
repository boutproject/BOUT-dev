#
# Import an entire NetCDF file into memory
#

from boututils import DataFile

def file_import(name):
    f = DataFile(name)   # Open file
    varlist = f.list() # Get list of all variables in file
    data = {}          # Create empty dictionary
    for v in varlist:
        data[v] = f.read(v)
    f.close()
    return data
