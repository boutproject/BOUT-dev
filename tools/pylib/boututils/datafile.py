"""
File I/O class
A wrapper around various NetCDF libraries, used by
BOUT++ routines. Creates a consistent interface
across machines

NOTE: NetCDF includes unlimited dimensions,
but this library is just for very simple
I/O operations. Educated guesses are made
for the dimensions.

Supported libraries:
-------------------

netCDF4

Scientific.IO.NetCDF

scipy.io.netcdf
  old version (create_dimension, create_variable)
  new version (createDimension, createVariable)
"""

from __future__ import print_function
try:
    from builtins import map
    from builtins import zip
    from builtins import str
    from builtins import object
except:
    pass

import numpy as np
import time, getpass

# Record which library to use
library = None

try:
    from netCDF4 import Dataset
    library = "netCDF4"
    has_netCDF = True
except ImportError:
    try:
        from Scientific.IO.NetCDF import NetCDFFile as Dataset
        from Scientific.N import Int, Float, Float32
        library = "Scientific"
        has_netCDF = True
    except ImportError:
        try:
            from scipy.io.netcdf import netcdf_file as Dataset
            library = "scipy"
            has_netCDF = True
        except:
            raise ImportError("DataFile: No supported NetCDF modules available")

try:
    import h5py
    has_h5py = True
except ImportError:
    has_h5py = False

class DataFile:
    impl = None
    def __init__(self, filename=None, write=False, create=False, format='NETCDF3_CLASSIC'):
        if filename is not None:
            if filename.split('.')[-1] in ('hdf5','hdf','h5'):
                self.impl = DataFile_HDF5(filename=filename, write=write, create=create, format=format)
            else:
                self.impl = DataFile_netCDF(filename=filename, write=write, create=create, format=format)
        elif format == 'HDF5':
            self.impl = DataFile_HDF5(filename=filename, write=write, create=create, format=format)
        else:
            self.impl = DataFile_netCDF(filename=filename, write=write, create=create, format=format)

    def open(self, filename, write=False, create=False,
             format='NETCDF3_CLASSIC'):
      self.impl.open(filename, write=write, create=create,
                  format=format)
    def close(self):
        self.impl.close()

    def __del__(self):
        if self.impl is not None:
            self.impl.__del__()

    def __enter__(self):
        return self.impl.__enter__()

    def __exit__(self, type, value, traceback):
        self.impl.__exit__(type, value, traceback)

    def read(self, name, ranges=None):
        """Read a variable from the file."""
        return self.impl.read(name, ranges=ranges)

    def list(self):
        """List all variables in the file."""
        return self.impl.list()

    def keys(self):
        return self.list()

    def dimensions(self, varname):
        """Array of dimension names"""
        return self.impl.dimensions(varname)

    def ndims(self, varname):
        """Number of dimensions for a variable."""
        return self.impl.ndims(varname)

    def size(self, varname):
        """List of dimension sizes for a variable."""
        return self.impl.size(varname)

    def write(self, name, data):
        """Writes a variable to file, making guesses for the dimensions"""
        return self.impl.write(name, data)

    def __getitem__(self, name):
        return self.impl.__getitem__(name)

    def __setitem__(self, key, value):
        self.impl.__setitem__(key, value)

    def attributes(self, varname):
        """Return a dictionary of attributes"""
        return self.impl.attributes(varname);
        
class DataFile_netCDF(DataFile):
    handle = None
    # Print warning if netcdf is used without the netcdf library
    if library != "netCDF4":
        print("WARNING: netcdf4-python module not found")
        print("         expect poor performance")
        if library == "Scientific":
            print("  => Using Scientific.IO.NetCDF instead")
        elif library == "scipy":
            print("  => Using scipy.io.netcdf instead")


    def open(self, filename, write=False, create=False,
             format='NETCDF3_CLASSIC'):
        if (not write) and (not create):
            if library == "scipy":
                self.handle = Dataset(filename, "r", mmap=False)
            else:
                self.handle = Dataset(filename, "r")
        elif create:
            if library == "Scientific":
                self.handle = Dataset(filename, "w",
                                      'Created ' + time.ctime(time.time())
                                      + ' by ' + getpass.getuser())
            elif library == "scipy":
                self.handle = Dataset(filename, "w")
            else:
                self.handle = Dataset(filename, "w", format=format)
        else:
            if library == "scipy":
                raise Exception("scipy.io.netcdf doesn't support appending");
            else:
                self.handle = Dataset(filename, "a")
        # Record if writing
        self.writeable = write or create

    def close(self):
        if self.handle is not None:
            self.handle.close()
        self.handle = None

    def __init__(self, filename=None, write=False, create=False,
                 format='NETCDF3_CLASSIC'):
        if not has_netCDF:
            message = "DataFile: No supported NetCDF python-modules available"
            raise ImportError(message)
        if filename is not None:
            self.open(filename, write=write, create=create, format=format)

    def __del__(self):
        self.close()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def read(self, name, ranges=None):
        """Read a variable from the file."""
        if self.handle is None: return None

        try:
            var = self.handle.variables[name]
        except KeyError:
            # Not found. Try to find using case-insensitive search
            var = None
            for n in list(self.handle.variables.keys()):
                if n.lower() == name.lower():
                    print("WARNING: Reading '"+n+"' instead of '"+name+"'")
                    var = self.handle.variables[n]
            if var is None:
                return None
        ndims = len(var.dimensions)
        if ndims == 0:
            data = var.getValue()
            return data #[0]
        else:
            if ranges is not None:
                if len(ranges) != 2*ndims:
                    print("Incorrect number of elements in ranges argument")
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
                        data = var[(ranges[0]):(ranges[1]),
                                   (ranges[2]):(ranges[3]),
                                   (ranges[4]):(ranges[5]),
                                   (ranges[6]):(ranges[7])]
                return data
            else:
                return var[:]

    def __getitem__(self, name):
        var = self.read(name)
        if var is None:
            raise KeyError("No variable found: "+name)
        return var

    def __setitem__(self, key, value):
        self.write(key, value)

    def list(self):
        """List all variables in the file."""
        if self.handle is None: return []
        return list(self.handle.variables.keys())

    def keys(self):
        """List all variables in the file."""
        return self.list()

    def dimensions(self, varname):
        """Array of dimension names"""
        if self.handle is None: return None
        try:
            var = self.handle.variables[varname]
        except KeyError:
            raise ValueError("No such variable")
        return var.dimensions

    def ndims(self, varname):
        """Number of dimensions for a variable."""
        if self.handle is None:
           raise ValueError("File not open")
        try:
            var = self.handle.variables[varname]
        except KeyError:
            raise ValueError("No such variable")
        return len(var.dimensions)

    def size(self, varname):
        """List of dimension sizes for a variable."""
        if self.handle is None: return []
        try:
            var = self.handle.variables[varname]
        except KeyError:
            return []

        def dimlen(d):
            dim = self.handle.dimensions[d]
            if dim is not None:
                t = type(dim).__name__
                if t == 'int':
                    return dim
                return len(dim)
            return 0
        return [dimlen(d) for d in var.dimensions]

    def write(self, name, data, info = False):
        """Writes a variable to file, making guesses for the dimensions"""

        if not self.writeable:
            raise Exception("File not writeable. Open with write=True keyword")

        s = np.shape(data)

        # Get the variable type
        t = type(data).__name__

        if t == 'NoneType':
            print("DataFile: None passed as data to write. Ignoring")
            return

        if t == 'ndarray':
            # Numpy type. Get the data type
            t = data.dtype.str

        if t == 'list':
            # List -> convert to numpy array
            data = np.array(data)
            t = data.dtype.str

        if (t == 'int') or (t == '<i8') or (t == 'int64') :
            # NetCDF 3 does not support type int64
            data = np.int32(data)
            t = data.dtype.str

        try:
            # See if the variable already exists
            var = self.handle.variables[name]

            # Check the shape of the variable
            if var.shape != s:
                print("DataFile: Variable already exists with different size: "+ name)
                # Fallthrough to the exception
                raise
        except:
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
                    if type(d).__name__ == 'int':
                        if d == size:
                            return name;
                    else:
                        if len(d) == size:
                            return name

                    # Find another with the correct size
                    for dn, d in list(self.handle.dimensions.items()):
                        # Some implementations need len(d) here, some just d
                        if type(d).__name__ == 'int':
                            if d == size:
                                return dn
                        else:
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
                            if info:
                                print("Defining dimension "+ dn + " of size %d" % size)
                            try:
                                self.handle.createDimension(dn, size)
                            except AttributeError:
                                # Try the old-style function
                                self.handle.create_dimension(dn, size)
                            return dn
                        i = i + 1

                except KeyError:
                    # Doesn't exist, so add
                    if info:
                        print("Defining dimension "+ name + " of size %d" % size)
                    try:
                        self.handle.createDimension(name, size)
                    except AttributeError:
                        self.handle.create_dimension(name, size)

                return name

            # List of (size, 'name') tuples
            dlist = list(zip(s, defdims[len(s)]))
            # Get new list of variables, and turn into a tuple
            dims = tuple( map(find_dim, dlist) )

            # Create the variable
            if library == "Scientific":
                if t == 'int' or t == '<i4' or t == 'int32':
                    tc = Int
                elif t=='<f4':
                    tc = Float32
                else:
                    tc = Float
                var = self.handle.createVariable(name, tc, dims)

            elif library == "scipy":
                try:
                    # New style functions
                    var = self.handle.createVariable(name, t, dims)
                except AttributeError:
                    # Old style functions
                    var = self.handle.create_variable(name, t, dims)
            else:
                var = self.handle.createVariable(name, t, dims)

            if var is None:
                raise Exception("Couldn't create variable")

        # Write the data

        try:
            # Some libraries allow this for arrays
            var.assignValue(data)
        except:
            # And some others only this
            var[:] = data

    def attributes(self, varname):
        """Return a dictionary of variable attributes"""
        if self.handle is None: return None
        try:
            var = self.handle.variables[varname]
        except KeyError:
            # Not found. Try to find using case-insensitive search
            var = None
            for n in list(self.handle.variables.keys()):
                if n.lower() == varname.lower():
                    print("WARNING: Reading '"+n+"' instead of '"+varname+"'")
                    var = self.handle.variables[n]
            if var is None:
                return None
            
        result = {} # Map of attribute names to values
        
        try:
            # This code tested with NetCDF4 library
            attribs = var.ncattrs() # List of attributes
            for attrname in attribs:
                result[attrname] = var.getncattr(attrname) # Get all values and insert into map
        except:
            print("Error reading attributes")
            # Result will be an empty map
            
        return result
        

class DataFile_HDF5(DataFile):
    handle = None

    def open(self, filename, write=False, create=False, format=None):
        if (not write) and (not create):
            self.handle = h5py.File(filename,mode="r")
        elif create:
            self.handle = h5py.File(filename,mode="w")
        else:
            self.handl = h5py.File(filename,mode="a")
        # Record if writing
        self.writeable = write or create

    def close(self):
        if self.handle is not None:
            self.handle.close()
        self.handle = None

    def __init__(self, filename=None, write=False, create=False,
                 format=None):
        if not has_h5py:
            message = "DataFile: No supported HDF5 python-modules available"
            raise ImportError(message)
        if filename is not None:
            self.open(filename, write=write, create=create, format=format)

    def __del__(self):
        self.close()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def read(self, name, ranges=None):
        """Read a variable from the file."""
        if self.handle is None: return None

        try:
            var = self.handle[name]
        except KeyError:
            # Not found. Try to find using case-insensitive search
            var = None
            for n in self.handle.variables.keys():
                if n.lower() == name.lower():
                    print("WARNING: Reading '"+n+"' instead of '"+name+"'")
                    var = self.handle[name]
            if var is None:
                return None
        ndims = len(var.shape)
        if ndims == 1 and var.shape[0] == 1:
            data = var
            return data[0]
        else:
            if ranges is not None:
                if len(ranges) != 2*ndims:
                    print("Incorrect number of elements in ranges argument")
                    return None

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
                    data = var[(ranges[0]):(ranges[1]),
                                (ranges[2]):(ranges[3]),
                                (ranges[4]):(ranges[5]),
                                (ranges[6]):(ranges[7])]
                return data
            else:
                return var[...]

    def __getitem__(self, name):
        var = self.read(name)
        if var is None:
            raise KeyError("No variable found: "+name)
        return var

    def __setitem__(self, key, value):
        self.write(key, value)
    
    def list(self):
        """List all variables in the file."""
        if self.handle is None: return []
        names = []
        self.handle.visit(names.append)
        return names

    def keys(self):
        """List all variables in the file."""
        return self.list()

    def dimensions(self, varname):
        """Array of dimension names"""
        var = self.handle[varname]
        vartype = var.attrs['type']
        if vartype == 'Field3D_t':
            return ('t','x','y','z')
        elif vartype == 'Field2D_t':
            return ('t','x','y')
        elif vartype == 'scalar_t':
            return ('t')
        elif vartype == 'Field3D':
            return ('x','y','z')
        elif vartype == 'Field2D':
            return ('x','y')
        elif vartype == 'scalar':
            return ()
        else:
            raise ValueError("Variable type not recognized")

    def ndims(self, varname):
        """Number of dimensions for a variable."""
        if self.handle is None: return None
        try:
            var = self.handle[varname]
        except KeyError:
            raise ValueError("Variable not found")
        return len(var.shape)

    def size(self, varname):
        """List of dimension sizes for a variable."""
        if self.handle is None: return None
        try:
            var = self.handle[varname]
        except KeyError:
            return None
        return var.shape

    def write(self, name, data):
        """Writes a variable to file"""

        if not self.writeable:
            raise Exception("File not writeable. Open with write=True keyword")

        self.handle.create_dataset(name, data=data)
        
    def attributes(self, varname):
        """Return a map of variable attributes"""
        return {} # Empty for now
    
