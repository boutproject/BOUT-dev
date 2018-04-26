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
import time
import getpass
from boututils.boutwarnings import alwayswarn
from boututils.boutarray import BoutArray

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
            raise ImportError(
                "DataFile: No supported NetCDF modules available")

try:
    import h5py
    has_h5py = True
except ImportError:
    has_h5py = False


class DataFile:
    impl = None

    def __init__(self, filename=None, write=False, create=False, format='NETCDF3_CLASSIC'):
        if filename is not None:
            if filename.split('.')[-1] in ('hdf5', 'hdf', 'h5'):
                self.impl = DataFile_HDF5(
                    filename=filename, write=write, create=create, format=format)
            else:
                self.impl = DataFile_netCDF(
                    filename=filename, write=write, create=create, format=format)
        elif format == 'HDF5':
            self.impl = DataFile_HDF5(
                filename=filename, write=write, create=create,
                format=format)
        else:
            self.impl = DataFile_netCDF(
                filename=filename, write=write, create=create, format=format)

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

    def read(self, name, ranges=None, asBoutArray=True):
        """Read a variable from the file."""
        return self.impl.read(name, ranges=ranges, asBoutArray=asBoutArray)

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

    def bout_type(self, varname):
        """BOUT++ type of a variable"""
        return self.attributes(varname)["bout_type"]


    def write(self, name, data, info=False):
        """Writes a variable to file, making guesses for the dimensions"""
        return self.impl.write(name, data, info)

    def __getitem__(self, name):
        return self.impl.__getitem__(name)

    def __setitem__(self, key, value):
        self.impl.__setitem__(key, value)

    def attributes(self, varname):
        """Return a dictionary of attributes"""
        return self.impl.attributes(varname)


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
                raise Exception("scipy.io.netcdf doesn't support appending")
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
        self._attributes_cache = {}

    def __del__(self):
        self.close()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def read(self, name, ranges=None, asBoutArray=True):
        """Read a variable from the file."""
        if self.handle is None:
            return None

        try:
            var = self.handle.variables[name]
            n = name
        except KeyError:
            # Not found. Try to find using case-insensitive search
            var = None
            for n in list(self.handle.variables.keys()):
                if n.lower() == name.lower():
                    print(
                        "WARNING: Reading '" + n + "' instead of '" + name + "'")
                    var = self.handle.variables[n]
            if var is None:
                return None

        if asBoutArray:
            attributes = self.attributes(n)

        ndims = len(var.dimensions)
        if ndims == 0:
            data = var.getValue()
            if asBoutArray:
                data = BoutArray(data, attributes=attributes)
            return data  # [0]
        else:
            if ranges is not None:
                if len(ranges) != 2 * ndims:
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
                if asBoutArray:
                    data = BoutArray(data, attributes=attributes)
                return data
            else:
                data = var[:]
                if asBoutArray:
                    data = BoutArray(data, attributes=attributes)
                return data

    def __getitem__(self, name):
        var = self.read(name)
        if var is None:
            raise KeyError("No variable found: " + name)
        return var

    def __setitem__(self, key, value):
        self.write(key, value)

    def list(self):
        """List all variables in the file."""
        if self.handle is None:
            return []
        return list(self.handle.variables.keys())

    def keys(self):
        """List all variables in the file."""
        return self.list()

    def dimensions(self, varname):
        """Array of dimension names"""
        if self.handle is None:
            return None
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
        if self.handle is None:
            return []
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

    def _bout_type_from_dimensions(self, varname):
        dims = self.dimensions(varname)
        if dims == ('t', 'x', 'y', 'z'):
            return "Field3D_t"
        elif dims == ('t', 'x', 'y'):
            return "Field2D_t"
        elif dims == ('t',):
            return "scalar_t"
        elif dims == ('x', 'y', 'z'):
            return "Field3D"
        elif dims == ('x', 'y'):
            return "Field2D"
        elif dims == ():
            return "scalar"
        else:
            # Unknown bout_type, but still want to be able to read, so give it a value...
            return None

    def write(self, name, data, info=False):
        """Writes a variable to file, making guesses for the dimensions"""

        if not self.writeable:
            raise Exception("File not writeable. Open with write=True keyword")

        s = np.shape(data)

        # Get the variable type
        t = type(data).__name__

        if t == 'NoneType':
            print("DataFile: None passed as data to write. Ignoring")
            return

        if t == 'ndarray' or t == 'BoutArray':
            # Numpy type or BoutArray wrapper for Numpy type. Get the data type
            t = data.dtype.str

        if t == 'list':
            # List -> convert to numpy array
            data = np.array(data)
            t = data.dtype.str

        if (t == 'int') or (t == '<i8') or (t == 'int64'):
            # NetCDF 3 does not support type int64
            data = np.int32(data)
            t = data.dtype.str

        try:
            # See if the variable already exists
            var = self.handle.variables[name]

            # Check the shape of the variable
            if var.shape != s:
                print(
                    "DataFile: Variable already exists with different size: " + name)
                # Fallthrough to the exception
                raise
        except:
            # Not found, so add.

            # Get dimensions
            defdims = [(),
                       ('t',),
                       ('x', 'y'),
                       ('x', 'y', 'z'),
                       ('t', 'x', 'y', 'z')]

            def find_dim(dim):
                # Find a dimension with given name and size
                size, name = dim

                # See if it exists already
                try:
                    d = self.handle.dimensions[name]

                    # Check if it's the correct size
                    if type(d).__name__ == 'int':
                        if d == size:
                            return name
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
                                print(
                                    "Defining dimension " + dn + " of size %d" % size)
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
                        print(
                            "Defining dimension " + name + " of size %d" % size)
                    if name == 't':
                        size = None
                    try:
                        self.handle.createDimension(name, size)
                    except AttributeError:
                        self.handle.create_dimension(name, size)

                return name

            # List of (size, 'name') tuples
            dlist = list(zip(s, defdims[len(s)]))
            # Get new list of variables, and turn into a tuple
            dims = tuple(map(find_dim, dlist))

            # Create the variable
            if library == "Scientific":
                if t == 'int' or t == '<i4' or t == 'int32':
                    tc = Int
                elif t == '<f4':
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

        # Write attributes, if present
        try:
            for attrname in data.attributes:
                var.setncattr(attrname, data.attributes[attrname])
        except AttributeError:
            pass

    def attributes(self, varname):
        """Return a dictionary of variable attributes"""
        try:
            return self._attributes_cache[varname]
        except KeyError:
            # Need to build the attributes dictionary for this variable
            if self.handle is None:
                return None
            try:
                var = self.handle.variables[varname]
            except KeyError:
                # Not found. Try to find using case-insensitive search
                var = None
                for n in list(self.handle.variables.keys()):
                    if n.lower() == varname.lower():
                        print(
                            "WARNING: Reading '" + n + "' instead of '" + varname + "'")
                        var = self.handle.variables[n]
                if var is None:
                    return None

            attributes = {}  # Map of attribute names to values

            try:
                # This code tested with NetCDF4 library
                attribs = var.ncattrs()  # List of attributes
                for attrname in attribs:
                    attributes[attrname] = var.getncattr(
                        attrname)  # Get all values and insert into map
            except:
                print("Error reading attributes for " + varname)
                # Result will be an empty map

            if not "bout_type" in attributes:
                attributes["bout_type"] = self._bout_type_from_dimensions(varname)

            # Save the attributes for this variable to the cache
            self._attributes_cache[varname] = attributes

            return attributes


class DataFile_HDF5(DataFile):
    handle = None

    def open(self, filename, write=False, create=False, format=None):
        if (not write) and (not create):
            self.handle = h5py.File(filename, mode="r")
        elif create:
            self.handle = h5py.File(filename, mode="w")
        else:
            self.handle = h5py.File(filename, mode="a")
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
        self._attributes_cache = {}

    def __del__(self):
        self.close()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def read(self, name, ranges=None, asBoutArray=True):
        """Read a variable from the file."""
        if self.handle is None:
            return None

        try:
            var = self.handle[name]
            n = name
        except KeyError:
            # Not found. Try to find using case-insensitive search
            var = None
            for n in self.handle:
                if n.lower() == name.lower():
                    print(
                        "WARNING: Reading '" + n + "' instead of '" + name + "'")
                    var = self.handle[n]
            if var is None:
                return None

        if asBoutArray:
            attributes = self.attributes(n)

        ndims = len(var.shape)
        if ndims == 1 and var.shape[0] == 1:
            data = var
            if asBoutArray:
                data = BoutArray(data, attributes=attributes)
            return data[0]
        else:
            if ranges is not None:
                if len(ranges) != 2 * ndims:
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
                if asBoutArray:
                    data = BoutArray(data, attributes=attributes)
                return data
            else:
                data = var[...]
                if asBoutArray:
                    data = BoutArray(data, attributes=attributes)
                return data

    def __getitem__(self, name):
        var = self.read(name)
        if var is None:
            raise KeyError("No variable found: " + name)
        return var

    def __setitem__(self, key, value):
        self.write(key, value)

    def list(self):
        """List all variables in the file."""
        if self.handle is None:
            return []
        names = []
        self.handle.visit(names.append)
        return names

    def keys(self):
        """List all variables in the file."""
        return self.list()

    def dimensions(self, varname):
        """Array of dimension names"""
        bout_type = self.bout_type(varname)
        if bout_type == 'Field3D_t':
            return ('t', 'x', 'y', 'z')
        elif bout_type == 'Field2D_t':
            return ('t', 'x', 'y')
        elif bout_type == 'scalar_t':
            return ('t')
        elif bout_type == 'Field3D':
            return ('x', 'y', 'z')
        elif bout_type == 'Field2D':
            return ('x', 'y')
        elif bout_type == 'scalar':
            return ()
        else:
            raise ValueError("Variable bout_type not recognized")

    def _bout_type_from_array(self, data):
        """
        Get the bout_type from the array 'data'. If 'data' is a BoutArray, it knows
        its bout_type, otherwise we have to guess.
        """
        try:
            # If data is a BoutArray, it should have a type attribute that we can use
            bout_type = data.attributes["bout_type"]
            return bout_type
        except AttributeError:
            # Otherwise data is a numpy.ndarray and we have to guess the bout_type
            pass

        try:
            ndim = len(data.shape)
        except AttributeError:
            ndim = 0
        if ndim == 4:
            return 'Field3D_t'
        elif ndim == 3:
            # not ideal, 3d field might be time-evolving 2d field,
            # 'Field2D_t', but can't think of a good way to distinguish
            alwayswarn("Warning: assuming bout_type of 3d array is Field3D. If it "
                       "should be a time-evolving Field2D, this may cause errors in "
                       "dimension sizes.")
            return 'Field3D'
        elif ndim == 2:
            return 'Field2D'
        elif ndim == 1:
            return 'scalar_t'
        elif ndim == 0:
            return 'scalar'
        else:
            raise ValueError("Unrecognized variable bout_type, ndims=" + str(ndim))

    def ndims(self, varname):
        """Number of dimensions for a variable."""
        if self.handle is None:
            return None
        try:
            var = self.handle[varname]
        except KeyError:
            raise ValueError("Variable not found")
        if var.size == 1:
            # variable is a scalar, but h5py always (?) returns numpy.ndarray,
            # so var.shape=(1,)
            return 0
        else:
            return len(var.shape)

    def size(self, varname):
        """List of dimension sizes for a variable."""
        if self.handle is None:
            return None
        try:
            var = self.handle[varname]
        except KeyError:
            return None
        return var.shape

    def write(self, name, data, info=False):
        """Writes a variable to file"""

        if not self.writeable:
            raise Exception("File not writeable. Open with write=True keyword")

        try:
            bout_type = data.attributes["bout_type"]
        except AttributeError:
            bout_type = self._bout_type_from_array(data)

        if info:
            print("Creating variable '" + name +
                  "' with bout_type '" + bout_type + "'")

        if bout_type in ["Field3D_t", "Field2D_t", "scalar_t"]:
            # time evolving fields
            shape = list(data.shape)
            # set time dimension to None to make unlimited
            shape[0] = None
            self.handle.create_dataset(name, data=data, maxshape=shape)
        elif bout_type == 'scalar':
            # Need to create scalars as one element arrays to be compatible
            # with BOUT++ assumptions (maybe it would be better to read/write
            # scalars in BOUT++?)
            self.handle.create_dataset(name, data=np.array([data]))
        else:
            self.handle.create_dataset(name, data=data)

        # Need encodes in the following to make sure we pass a byte-string to
        # attrs and not a regular python string.

        # Check if the bout_type of the variable will be written when copying
        # attributes from data, which it should be if data is a BoutArray.
        # Otherwise, need to write it explicitly
        try:
            if (not "bout_type" in data.attributes):
                raise AttributeError("'bout_type' not found in attributes")
        except AttributeError:
            self.handle[name].attrs.create(
                'bout_type', bout_type.encode(encoding='utf-8'))

        try:
            for attrname in data.attributes:
                attrval = data.attributes[attrname]
                if type(attrval == str):
                    attrval = attrval.encode(encoding='utf-8')
                self.handle[name].attrs.create(attrname, attrval)
        except AttributeError:
            # data is not a BoutArray, so doesn't have attributes to write
            pass

    def attributes(self, varname):
        """Return a map of variable attributes"""

        try:
            return self._attributes_cache[varname]
        except KeyError:
            # Need to add attributes for this variable to the cache
            attributes = {}
            var = self.handle[varname]
            for attrname in var.attrs:
                attribute = var.attrs[attrname]
                if type(attribute) in [bytes, np.bytes_]:
                    attribute = str(attribute, encoding="utf-8")
                attributes[attrname] = attribute

            if not "bout_type" in attributes:
                # bout_type is a required attribute for BOUT++ outputs, so it should
                # have been found
                raise ValueError(
                    "Error: bout_type not found in attributes of "+varname)

            # Save the attributes for this variable to the cache
            self._attributes_cache[varname] = attributes

            return attributes
