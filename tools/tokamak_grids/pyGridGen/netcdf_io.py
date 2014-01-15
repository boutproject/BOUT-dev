import numpy
from netCDF4 import Dataset


def file_open(output):
    out = Dataset(output,'w', format='NETCDF3_CLASSIC')
    return out
    

def file_write(out,name,var):

    format={int:'i', numpy.int64:'i4', float:'f', numpy.float64:'f8'}

    ns=numpy.ndim(var) 
    if ns==0 :
            out.createDimension(name,1)
             
            if type(var) in format.keys() : 
                typ=format.get(type(var))
            else:
                print 'Unknown format'    

            v0=out.createVariable(name,typ,(name))
            v0[:]=var
            
    elif ns==1:
        
            out.createDimension(name, numpy.size(var))
            
            if type(var[0]) in format.keys() : 
                typ=format.get(type(var))
            else:
                print 'Unknown format'    

            v0=out.createVariable(name,typ,(name))
            v0[:]=var

    elif ns==2:
            ni=0
            for dimobj in out.dimensions.items():
                if dimobj[0]=='x' : ni=ni+1
                if dimobj[0]=='y' : ni=ni+1
            if ni==0 :
                out.createDimension('x',numpy.shape(var)[0])
                out.createDimension('y',numpy.shape(var)[1])

            if type(var[0,0]) in format.keys() : 
                typ=format.get(type(var))
            else:
                print 'Unknown format'    
             
            v1 = out.createVariable(name,typ,('x','y'))
            v1[:,:] = var
        
        
        
def file_close(out):
    out.close()
    
