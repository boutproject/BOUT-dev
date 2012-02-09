# File I/O class
# A wrapper around the NetCDF4 library, used by
# BOUT++ routines. This allows easily changing
# methods later.
#
# NOTE: NetCDF includes unlimited dimensions,
# but this library is just for very simple
# I/O operations. Educated guesses are made
# for the dimensions.

try:
    import numpy as np
except ImportError:
    print "ERROR: NumPy module not available"
    raise

library = None # Record which library to use

try:
    from netCDF4 import Dataset
    library = "netCDF4"
except ImportError:
    #print "netcdf4-python module not found"
    
    try:
        from Scientific.IO.NetCDF import NetCDFFile as Dataset
        from Scientific.N import Int, Float
        library = "Scientific"
        #print "  => Using Scientific.IO.NetCDF instead"
    except ImportError:
        print "DataFile: Neither netCDF4 nor Scientific.IO.NetCDF modules available"
        raise

import time

def getUserName():
    try:
	import os, pwd, string
    except ImportError:
	return 'unknown user'
    pwd_entry = pwd.getpwuid(os.getuid())
    name = string.strip(string.splitfields(pwd_entry[4], ',')[0])
    if name == '':
	name = pwd_entry[0]
    return name

class DataFile:
    handle = None

    def open(self, filename, write=False, create=False,
             format='NETCDF3_CLASSIC'):
        if (not write) and (not create):
            self.handle = Dataset(filename, "r")
        elif create:
            if library == "Scientific":
                self.handle = Dataset(filename, "w", 
                                      'Created ' + time.ctime(time.time())
                                      + ' by ' + getUserName())
            else:
                self.handle = Dataset(filename, "w", format=format)
        else:
            self.handle = Dataset(filename, "a")
    
    def close(self):
        if self.handle != None:
            self.handle.close()
        self.handle = None
    
    def __init__(self, filename=None):
        if filename != None:
            self.open(filename)
    
    def __del__(self):
        self.close()

    def read(self, name, ranges=None):
        """Read a variable from the file."""
        if self.handle == None: return None
        
        try:
            var = self.handle.variables[name]
        except KeyError:
            # Not found. Try to find using case-insensitive search
            var = None
            for n in self.handle.variables.keys():
                if n.lower() == name.lower():
                    print "WARNING: Reading '"+n+"' instead of '"+name+"'"
                    var = self.handle.variables[n]
            if var == None:
                return None
        ndims = len(var.dimensions)
        if ndims == 0:
            data = var.getValue()
            return data #[0]
        else:
            if ranges != None:
                if len(ranges) != 2*ndims:
                    print "Incorrect number of elements in ranges argument"
                    return None
                
                if library == "Scientific":
                    # Passing ranges to var[] doesn't seem to work
                    data = var[:]
                    if ndims == 1:
                        data = data[ranges[0]:ranges[1]]
                    elif ndims == 2:
                        data = data[ranges[0]:ranges[1], 
                                    ranges[2]:ranges[3]]
                    elif ndims == 3:
                        data = data[ranges[0]:ranges[1], 
                                    ranges[2]:ranges[3],
                                    ranges[4]:ranges[5]]
                    elif ndims == 4:
                        data = data[(ranges[0]):(ranges[1]),
                                    (ranges[2]):(ranges[3]),
                                    (ranges[4]):(ranges[5]),
                                    (ranges[6]):(ranges[7])]
                else:
                    if ndims == 1:
                        data = var[ranges[0]:ranges[1]]
                    elif ndims == 2:
                        data = var[ranges[0]:ranges[1], 
                                   ranges[2]:ranges[3]]
                    elif ndims == 3:
                        data = var[ranges[0]:ranges[1], 
                                   ranges[2]:ranges[3],
                                   ranges[4]:ranges[5]]
                    elif ndims == 4:
                        print "Ranges = ", ranges
                        data = var[(ranges[0]):(ranges[1]),
                                   (ranges[2]):(ranges[3]),
                                   (ranges[4]):(ranges[5]),
                                   (ranges[6]):(ranges[7])]
                return data
            else:
                return var[:]

    def list(self):
        """List all variables in the file."""
        if self.handle == None: return []
        return self.handle.variables.keys()

    def dimensions(self, varname):
        """Array of dimension names"""
        if self.handle == None: return None
        try:
            var = self.handle.variables[varname]
        except KeyError:
            return None
        return var.dimensions
        
    def ndims(self, varname):
        """Number of dimensions for a variable."""
        if self.handle == None: return None
        try:
            var = self.handle.variables[varname]
        except KeyError:
            return None
        return len(var.dimensions)
    
    def size(self, varname):
        """List of dimension sizes for a variable."""
        if self.handle == None: return []
        try:
            var = self.handle.variables[varname]
        except KeyError:
            return []
        
        def dimlen(d):
            dim = self.handle.dimensions[d]
            if dim != None:
                return len(dim)
            return 0
        return map(lambda d: dimlen(d), var.dimensions)

    def write(self, name, data):
        """Writes a variable to file, making guesses for the dimensions"""
        s = np.shape(data)

        # Get the variable type
        t = type(data).__name__
        
        if t == 'NoneType':
            print "DataFile: None passed as data to write. Ignoring"
            return
        
        if t == 'ndarray':
            # Numpy type
            t = data.dtype.str
        
        try:
            # See if the variable already exists
            var = self.handle.variables[name]

            # Check the shape of the variable
            if var.shape != s:
                print "Datafile: Variable already exists with different size: "+ name
                raise
        except KeyError:
            # Not found, so add.

            # Get dimensions
            defdims = [(),
                       ('x',),
                       ('x','y'),
                       ('x','y','z'),
                       ('t','x','y','z')]

            def find_dim(dim):
                # Find a dimension with given name and size
                size, name = dim

                # See if it exists already
                try:
                    d = self.handle.dimensions[name]

                    # Check if it's the correct size
                    if len(d) == size:
                        return name

                    # Find another with the correct size
                    for dn, d in self.handle.dimensions.iteritems():
                        if len(d) == size:
                            return dn

                    # None found, so create a new one
                    i = 2
                    while True:
                        dn = name + str(i)
                        try:
                            d = self.handle.dimensions[dn]
                            # Already exists, so keep going
                        except KeyError:
                            # Not found. Create
                            print "Defining dimension "+ dn + " of size %d" % size
                            self.handle.createDimension(dn, size)
                            return dn
                        i = i + 1
                    
                except KeyError:
                    # Doesn't exist, so add
                    print "Defining dimension "+ name + " of size %d" % size
                    self.handle.createDimension(name, size)
                return name
                
            # List of (size, 'name') tuples
            dlist = zip(s, defdims[len(s)])
            # Get new list of variables, and turn into a tuple
            dims = tuple( map(find_dim, dlist) )
            
            # Create the variable
            if library == "Scientific":
                if t == 'int':
                    tc = Int
                else:
                    tc = Float
                var = self.handle.createVariable(name, tc, dims)
            else:
                var = self.handle.createVariable(name, t, dims)

        # Write the data
        var.assignValue(data)
        
            
