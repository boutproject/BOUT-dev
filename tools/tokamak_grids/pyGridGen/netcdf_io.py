import numpy

from boututils import DataFile

def file_open(output, format='NETCDF3_CLASSIC'):
    return DataFile(output, write=True, create=True,format=format)
    
def file_write(out,name,var):
    out.write(name, var)
        
def file_close(out):
    out.close()
    
