"""File I/O class

A wrapper around various NetCDF libraries and h5py, used by BOUT++
routines. Creates a consistent interface across machines

Supported libraries:

- ``h5py`` (for HDF5 files)
- ``netCDF4`` (preferred NetCDF library)

NOTE
----
NetCDF and HDF5 include unlimited dimensions, but this library is just
for very simple I/O operations. Educated guesses are made for the
dimensions.

TODO
----
- Don't raise ``ImportError`` if no NetCDF libraries found, use HDF5
  instead?
- Cleaner handling of different NetCDF libraries
- Support for h5netcdf?

"""

from __future__ import print_function
from builtins import map, zip, str, object

import numpy as np
import time
import getpass
from boututils.boutwarnings import alwayswarn
from boututils.boutarray import BoutArray

try:
    from netCDF4 import Dataset
    has_netCDF = True
except ImportError:
    raise ImportError(
        "DataFile: No supported NetCDF modules available -- requires netCDF4")

try:
    import h5py
    has_h5py = True
except ImportError:
    has_h5py = False


class DataFile(object):
    """File I/O class

    A wrapper around various NetCDF libraries and h5py, used by BOUT++
    routines. Creates a consistent interface across machines

    Parameters
    ----------
    filename : str, optional
        Name of file to open. If no filename supplied, you will need
        to call :py:obj:`~DataFile.open` and supply `filename` there
    write : bool, optional
        If True, open the file in read-write mode (existing files will
        be appended to). Default is read-only mode
    create : bool, optional
        If True, open the file in write mode (existing files will be
        truncated). Default is read-only mode
    format : str, optional
        Name of a filetype to use (e.g. ``NETCDF3_CLASSIC``,
        ``NETCDF3_64BIT``, ``NETCDF4``, ``HDF5``)

    TODO
    ----
    - `filename` should not be optional!
    - Take a ``mode`` argument to be more in line with other file types
    - `format` should be checked to be a sensible value
    - Make sure ``__init__`` methods are first
    - Make `impl` and `handle` private

    """
    impl = None

    def __init__(self, filename=None, write=False, create=False, format='NETCDF3_64BIT', **kwargs):
        """

        NetCDF formats are described here: http://unidata.github.io/netcdf4-python/
        - NETCDF3_CLASSIC   Limited to 2.1Gb files
        - NETCDF3_64BIT_OFFSET or NETCDF3_64BIT is an extension to allow larger file sizes
        - NETCDF3_64BIT_DATA adds 64-bit integer data types and 64-bit dimension sizes
        - NETCDF4 and NETCDF4_CLASSIC use HDF5 as the disk format
        """
        if filename is not None:
            if filename.split('.')[-1] in ('hdf5', 'hdf', 'h5'):
                self.impl = DataFile_HDF5(
                    filename=filename, write=write, create=create, format=format)
            else:
                self.impl = DataFile_netCDF(
                    filename=filename, write=write, create=create, format=format, **kwargs)
        elif format == 'HDF5':
            self.impl = DataFile_HDF5(
                filename=filename, write=write, create=create,
                format=format)
        else:
            self.impl = DataFile_netCDF(
                filename=filename, write=write, create=create, format=format, **kwargs)

    def open(self, filename, write=False, create=False,
             format='NETCDF3_CLASSIC'):
        """Open the file

        Parameters
        ----------
        filename : str, optional
            Name of file to open
        write : bool, optional
            If True, open the file in read-write mode (existing files will
            be appended to). Default is read-only mode
        create : bool, optional
            If True, open the file in write mode (existing files will be
            truncated). Default is read-only mode
        format : str, optional
            Name of a filetype to use (e.g. ``NETCDF3_CLASSIC``,
            ``NETCDF4``, ``HDF5``)

        TODO
        ----
        - Return the result of calling open to be more like stdlib's
          open
        - `keys` should be more pythonic (return generator)

        """
        self.impl.open(filename, write=write, create=create,
                       format=format)

    def close(self):
        """Close a file and flush data to disk

        """
        self.impl.close()

    def __del__(self):
        if self.impl is not None:
            self.impl.__del__()

    def __enter__(self):
        self.impl.__enter__()
        return self

    def __exit__(self, type, value, traceback):
        self.impl.__exit__(type, value, traceback)

    def read(self, name, ranges=None, asBoutArray=True):
        """Read a variable from the file

        Parameters
        ----------
        name : str
            Name of the variable to read
        ranges : list of slice objects, optional
            Slices of variable to read, can also be converted from lists or
            tuples of (start, stop, stride). The number of elements in `ranges`
            should be equal to the number of dimensions of the variable you
            wish to read. See :py:obj:`~DataFile.size` for how to get the
            dimensions
        asBoutArray : bool, optional
            If True, return the variable as a
            :py:obj:`~boututils.boutarray.BoutArray` (the default)

        Returns
        -------
        ndarray or :py:obj:`~boututils.boutarray.BoutArray`
            The variable from the file
            (:py:obj:`~boututils.boutarray.BoutArray` if `asBoutArray`
            is True)

        """
        if ranges is not None:
            for x in ranges:
                if isinstance(x, (list, tuple)):
                    x = slice(*x)
        return self.impl.read(name, ranges=ranges, asBoutArray=asBoutArray)

    def list(self):
        """List all variables in the file

        Returns
        -------
        list of str
            A list containing all the names of the variables

        """
        return self.impl.list()

    def keys(self):
        """A synonym for :py:obj:`~DataFile.list`

        TODO
        ----
        - Make a generator to be more like python3 dict keys

        """
        return self.list()

    def dimensions(self, varname):
        """Return the names of all the dimensions of a variable

        Parameters
        ----------
        varname : str
            The name of the variable

        Returns
        -------
        tuple of str
            The names of the variable's dimensions

        """
        return self.impl.dimensions(varname)

    def ndims(self, varname):
        """Return the number of dimensions for a variable

        Parameters
        ----------
        varname : str
            The name of the variable

        Returns
        -------
        int
            The number of dimensions

        """
        return self.impl.ndims(varname)

    def sync(self):
        """Write pending changes to disk.

        """
        self.impl.sync()

    def size(self, varname):
        """Return the size of each dimension of a variable

        Parameters
        ----------
        varname : str
            The name of the variable

        Returns
        -------
        tuple of int
            The size of each dimension

        """
        return self.impl.size(varname)

    def bout_type(self, varname):
        """Return the name of the BOUT++ type of a variable

        Possible values are:

        - scalar
        - Field2D
        - Field3D

        If the variable is an evolving variable (i.e. has a time
        dimension), then it is appended with a "_t"

        Parameters
        ----------
        varname : str
            The name of the variable

        Returns
        -------
        str
            The name of the BOUT++ type

        """
        return self.attributes(varname)["bout_type"]

    def write(self, name, data, info=False):
        """Write a variable to file

        If the variable is not a :py:obj:`~boututils.boutarray.BoutArray` with
        the ``bout_type`` attribute, a guess will be made for the
        dimensions

        Parameters
        ----------
        name : str
            Name of the variable to use in the file
        data : :py:obj:`~boututils.boutarray.BoutArray` or ndarray
            An array containing the variable data
        info : bool, optional
            If True, print information about what is being written to
            file

        Returns
        -------
        None

        """
        return self.impl.write(name, data, info)

    def __getitem__(self, name):
        return self.impl.__getitem__(name)

    def __setitem__(self, key, value):
        self.impl.__setitem__(key, value)

    def attributes(self, varname):
        """Return a dictionary of attributes

        Parameters
        ----------
        varname : str
            The name of the variable

        Returns
        -------
        dict
            The attribute names and their values

        """
        return self.impl.attributes(varname)


class DataFile_netCDF(DataFile):
    handle = None

    def open(self, filename, write=False, create=False,
             format='NETCDF3_CLASSIC'):
        if (not write) and (not create):
            self.handle = Dataset(filename, "r")
        elif create:
            self.handle = Dataset(filename, "w", format=format)
        else:
            self.handle = Dataset(filename, "a")
        # Record if writing
        self.writeable = write or create

    def close(self):
        if self.handle is not None:
            self.handle.close()
        self.handle = None

    def __init__(self, filename=None, write=False, create=False,
                 format='NETCDF3_CLASSIC', **kwargs):
        self._kwargs = kwargs
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
            if ranges:
                if len(ranges) == 2 * ndims:
                    # Reform list of pairs of ints into slices
                    ranges = [slice(a, b) for a, b in
                              zip(ranges[::2], ranges[1::2])]
                elif len(ranges) != ndims:
                    raise ValueError("Incorrect number of elements in ranges argument "
                                     "(got {}, expected {} or {})"
                                     .format(len(ranges), ndims, 2 * ndims))

                data = var[ranges[:ndims]]
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
        if self.handle is None:
            return []
        return list(self.handle.variables.keys())

    def keys(self):
        return self.list()

    def dimensions(self, varname):
        if self.handle is None:
            return None
        try:
            var = self.handle.variables[varname]
        except KeyError:
            raise ValueError("No such variable")
        return var.dimensions

    def ndims(self, varname):
        if self.handle is None:
            raise ValueError("File not open")
        try:
            var = self.handle.variables[varname]
        except KeyError:
            raise ValueError("No such variable")
        return len(var.dimensions)

    def sync(self):
        self.handle.sync()

    def size(self, varname):
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
        dims_dict = {
            ('t', 'x', 'y', 'z'): "Field3D_t",
            ('t', 'x', 'y'): "Field2D_t",
            ('t',): "scalar_t",
            ('x', 'y', 'z'): "Field3D",
            ('x', 'y'): "Field2D",
            (): "scalar",
        }

        return dims_dict.get(dims, None)

    def write(self, name, data, info=False):

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
                                print("Defining dimension {} of size {}"
                                      .format(dn, size))

                            self.handle.createDimension(dn, size)
                            return dn
                        i = i + 1

                except KeyError:
                    # Doesn't exist, so add
                    if info:
                        print(
                            "Defining dimension " + name + " of size %d" % size)
                    if name == 't':
                        size = None

                    self.handle.createDimension(name, size)

                return name

            # List of (size, 'name') tuples
            dlist = list(zip(s, defdims[len(s)]))
            # Get new list of variables, and turn into a tuple
            dims = tuple(map(find_dim, dlist))

            # Create the variable
            var = self.handle.createVariable(name, t, dims, **self._kwargs)

            if var is None:
                raise Exception("Couldn't create variable")

        # Write the data
        if t == 'str':
            var[0] = data
        else:
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
                attribs = var.ncattrs()  # List of attributes
                for attrname in attribs:
                    attributes[attrname] = var.getncattr(
                        attrname)  # Get all values and insert into map
            except:
                print("Error reading attributes for " + varname)
                # Result will be an empty map

            if "bout_type" not in attributes:
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
            if ranges:
                if len(ranges) == 2 * ndims:
                    # Reform list of pairs of ints into slices
                    ranges = [slice(a, b) for a, b in
                              zip(ranges[::2], ranges[1::2])]
                elif len(ranges) != ndims:
                    raise ValueError("Incorrect number of elements in ranges argument "
                                     "(got {}, expected {} or {})"
                                     .format(len(ranges), ndims, 2 * ndims))
                # Probably a bug in h5py, work around by passing tuple
                data = var[tuple(ranges[:ndims])]
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
        if self.handle is None:
            return []
        names = []
        self.handle.visit(names.append)
        return names

    def keys(self):
        return self.list()

    def dimensions(self, varname):
        bout_type = self.bout_type(varname)
        dims_dict = {
            "Field3D_t": ('t', 'x', 'y', 'z'),
            "Field2D_t": ('t', 'x', 'y'),
            "scalar_t": ('t',),
            "Field3D": ('x', 'y', 'z'),
            "Field2D": ('x', 'y'),
            "scalar": (),
        }
        try:
            return dims_dict[bout_type]
        except KeyError:
            raise ValueError("Variable bout_type not recognized (got {})"
                             .format(bout_type))

    def _bout_type_from_array(self, data):
        """Get the bout_type from the array 'data'

        If 'data' is a BoutArray, it knows its bout_type, otherwise we
        have to guess.

        Parameters
        ----------
        data : :py:obj:`~boututils.boutarray.BoutArray` or ndarray
            An array with between 0 and 4 dimensions

        Returns
        -------
        str
            Either the actual bout_type or our best guess

        See Also
        --------
        - `DataFile.bout_type`

        TODO
        ----
        - Make standalone function

        """
        try:
            # If data is a BoutArray, it should have a type attribute that we can use
            return data.attributes["bout_type"]
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

    def sync(self):
        self.handle.flush()

    def size(self, varname):
        if self.handle is None:
            return None
        try:
            var = self.handle[varname]
        except KeyError:
            return None
        return var.shape

    def write(self, name, data, info=False):

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
