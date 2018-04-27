# Provides a class BoutData which makes access to code
# inputs and outputs easier. Creates a tree of maps,
# inspired by approach used in OMFIT
#
# 

import os
import sys
import glob
import numpy

from boutdata.collect import findVar as collect_findVar, findFiles as collect_findFiles

try:
    from boututils.datafile import DataFile
except ImportError:
    print("ERROR: boututils.datafile.DataFile couldn't be loaded")
    raise
    
class BoutOptions(object):
    """
    This class represents a tree structure.
    Each node (BoutOptions object) can have several
    sub-nodes (sections), and several key-value pairs.

    Example
    -------
    
    optRoot = BoutOptions()  # Create a root
    
    # Specify value of a key in a section "test" 
    # If the section does not exist then it is created
    
    optRoot.getSection("test")["key"] = value
    
    # Get the value of a key in a section "test"
    # If the section does not exist then a KeyError is raised
    
    print optRoot["test"]["key"]
    
    # To pretty print the options
    
    print optRoot
    
    """
    def __init__(self, name="root", parent=None):
        self._sections = {}
        self._keys = {}
        self._name = name
        self._parent = parent

    def getSection(self, name):
        """
        Return a section object. If the section
        does not exist then it is created
        """
        name = name.lower()
        
        if name in self._sections:
            return self._sections[name]
        else:
            newsection = BoutOptions(name, self)
            self._sections[name] = newsection
            return newsection
        
    def __getitem__(self, key):
        """
        First check if it's a section, then a value
        """
        key = key.lower()
        if key in self._sections:
            return self._sections[key]
            
        if key not in self._keys:
            raise KeyError("Key '%s' not in section '%s'" % (key, self.path()))
        return self._keys[key]

    def __setitem__(self, key, value):
        """
        Set a key
        """
        if len(key) == 0:
            return
        self._keys[key.lower()] = value

    def path(self):
        """
        Returns the path of this section,
        joining together names of parents
        """
        
        if self._parent:
            return self._parent.path() + ":" + self._name
        return self._name

    def keys(self):
        """
        Returns all keys, including sections and values
        """
        return self._sections.keys() + self._keys.keys()

    def sections(self):
        """
        Return a list of sub-sections
        """
        return self._sections.keys()

    def values(self):
        """
        Return a list of values
        """
        return self._keys.keys()

    def __len__(self):
        return len(self._sections) + len(self._keys)

    def __iter__(self):
        """
        Iterates over all keys. First values, then sections
        """
        for k in self._keys:
            yield k
        for s in self._sections:
            yield s

    def __str__(self, indent=""):
        """
        Print a pretty version of the options tree
        """
        text = self._name + "\n"
        
        for k in self._keys:
            text += indent + " |- " + k + " = " + str(self._keys[k]) + "\n"
        
        for s in self._sections:
            text += indent + " |- " + self._sections[s].__str__(indent+" |  ")
        return text

class BoutOptionsFile(BoutOptions):
    """
    Parses a BOUT.inp configuration file, producing
    a tree of BoutOptions.
    
    Slight differences from ConfigParser, including allowing
    values before the first section header.

    Example
    -------
    
    opts = BoutOptionsFile("BOUT.inp")
    
    print opts   # Print all options in a tree
    
    opts["All"]["scale"] # Value "scale" in section "All"
    
    """
    def __init__(self, filename, name="root"):
        BoutOptions.__init__(self, name)
        # Open the file
        with open(filename, "r") as f:
            # Go through each line in the file
            section = self # Start with root section
            for linenr, line in enumerate(f.readlines()):
                # First remove comments, either # or ;
                startpos = line.find("#")
                if startpos != -1:
                    line = line[:startpos]
                startpos = line.find(";")
                if startpos != -1:
                    line = line[:startpos]
                
                # Check section headers
                startpos = line.find("[")
                endpos = line.find("]")
                if startpos != -1:
                    # A section heading
                    if endpos == -1:
                        raise SyntaxError("Missing ']' on line %d" % (linenr,))
                    line = line[(startpos+1):endpos].strip()
                    
                    section = self
                    while True:
                        scorepos = line.find(":")
                        if scorepos == -1:
                            break
                        sectionname = line[0:scorepos]
                        line = line[(scorepos+1):]
                        section = section.getSection(sectionname)
                    section = section.getSection(line)
                else:
                    # A key=value pair
                    
                    eqpos = line.find("=")
                    if eqpos == -1:
                        # No '=', so just set to true
                        section[line.strip()] = True
                    else:
                        value = line[(eqpos+1):].strip()
                        try:
                            # Try to convert to an integer
                            value = int(value)
                        except ValueError:
                            try:
                                # Try to convert to float
                                value = float(value)
                            except ValueError:
                                # Leave as a string
                                pass
                        
                        section[line[:eqpos].strip()] = value
                    
                        
            
class BoutOutputs(object):
    """
    Emulates a map class, represents the contents of a BOUT++ dmp files. Does
    not allow writing, only reading of data.  By default there is no cache, so
    each time a variable is read it is collected; if caching is set to True
    variables are stored once they are read.  Extra keyword arguments are
    passed through to collect.
    
    Example
    -------

    d = BoutOutputs(".")  # Current directory
    
    d.keys()     # List all valid keys

    d.dimensions["ne"] # Get the dimensions of the field ne
    
    d["ne"] # Read "ne" from data files

    d = BoutOutputs(".", prefix="BOUT.dmp", caching=True) # Turn on caching

    Options
    -------
    prefix - sets the prefix for data files (default "BOUT.dmp")

    caching - switches on caching of data, so it is only read into memory when
              first accessed (default False) If caching is set to a number, it
              gives the maximum size of the cache in GB, after which entries
              will be discarded in first-in-first-out order to prevent the
              cache getting too big.  If the variable being returned is bigger
              than the maximum cache size, then the variable will be returned
              without being added to the cache, and the rest of the cache will
              be left.

    **kwargs - keyword arguments that are passed through to _caching_collect()
    """
    def __init__(self, path=".", prefix="BOUT.dmp", caching=False, DataFileCaching=True, **kwargs):
        """
        Initialise BoutOutputs object
        """
        self._path = path
        self._prefix = prefix
        self._caching = caching
        self._DataFileCaching = DataFileCaching
        self._kwargs = kwargs
        
        # Label for this data
        self.label = path

        # Check that the path contains some data
        # Search for BOUT++ dump files
        self._file_list,self._parallel,self._suffix=collect_findFiles(path,prefix)
        if len(self._file_list) == 0:
            raise ValueError("ERROR: No data files found")
        
        # Available variables
        self.varNames = []
        self.dimensions = {}
        self.evolvingVariableNames = []

        # Private variables
        if self._caching:
            from collections import OrderedDict
            self._datacache = OrderedDict()
            if self._caching is not True:
                # Track the size of _datacache and limit it to a maximum of _caching
                try:
                    # Check that _caching is a number of some sort
                    float(self._caching)
                except ValueError:
                    raise ValueError("BoutOutputs: Invalid value for caching argument. Caching should be either a number (giving the maximum size of the cache in GB), True for unlimited size or False for no caching.")
                self._datacachesize = 0
                self._datacachemaxsize = self._caching*1.e9

        if self._DataFileCaching:
            self._DataFileCache = []
        
        with DataFile(self._file_list[0]) as f:
            # Get variable names
            self.varNames = f.keys()
            for name in f.keys():
                dimensions = f.dimensions(name)
                self.dimensions[name] = dimensions
                if name != "t_array" and "t" in dimensions:
                    self.evolvingVariableNames.append(name)
        
    def keys(self):
        """
        Return a list of available variable names
        """
        return self.varNames

    def evolvingVariables(self):
        """
        Return a list of names of time-evolving variables
        """
        return self.evolvingVariableNames
        
    def __len__(self):
        return len(self.varNames)
            
    def __getitem__(self, name):
        """
        Reads a variable using _caching_collect.
        Caches result and returns later if called again, if self._caching=True
        
        """
        if self._caching:
            if name not in self._datacache.keys():
                item = self._caching_collect(name, path=self._path, prefix=self._prefix, **self._kwargs)
                if self._caching is not True:
                    itemsize = item.nbytes
                    if itemsize>self._datacachemaxsize:
                        return item
                    self._datacache[name] = item
                    self._datacachesize += itemsize
                    while self._datacachesize > self._datacachemaxsize:
                        self._removeFirstFromCache()
                else:
                    self._datacache[name] = item
                return item
            else:
                return self._datacache[name]
        else:
            # Collect the data from the repository
            data = self._caching_collect(name, path=self._path, prefix=self._prefix, **self._kwargs)
            return data
    
    def _removeFirstFromCache(self):
        # pop the first item from the OrderedDict _datacache
        item = self._datacache.popitem(last=False)
        self._datacachesize -= item[1].nbytes

    def _getDataFile(self, i):
        """
        Get a DataFile, creating the cache if necessary
        """
        if self._DataFileCaching:
            try:
                f = self._DataFileCache[i]
            except IndexError:
                for filename in self._file_list:
                    self._DataFileCache.append(DataFile(filename))
                f = self._DataFileCache[i]
        else:
            f = DataFile(self._file_list[i])
        return f

    def _caching_collect(self,varname, xind=None, yind=None, zind=None, tind=None, path=".",yguards=False, xguards=True, info=True,prefix="BOUT.dmp",strict=False,tind_auto=False):
        """Collect a variable from a set of BOUT++ outputs.

        data = self._caching_collect(name)

        varname   Name of the variable (string)

        Optional arguments:

        xind = [min,max]   Range of X indices to collect
        yind = [min,max]   Range of Y indices to collect
        zind = [min,max]   Range of Z indices to collect
        tind = [min,max]   Range of T indices to collect

        path    = "."          Path to data files
        prefix  = "BOUT.dmp"   File prefix
        yguards = False        Collect Y boundary guard cells?
        xguards = True         Collect X boundary guard cells?
                               (Set to True to be consistent with the
                               definition of nx)
        info    = True         Print information about _caching_collect?
        strict  = False        Fail if the exact variable name is not found?
        tind_auto = False      Read all files, to get the shortest length of time_indices
                               useful if writing got interrupted.
        """
        if self._parallel:
            print("Single (parallel) data file")
            f = self._getDataFile(0) # Open the file

            data = f.read(varname)
            return data
        nfiles = len(self._file_list)

        # Read data from the first file
        f = self._getDataFile(0)

        try:
            dimens = f.dimensions(varname)
            #ndims = len(dimens)
            ndims = f.ndims(varname)
        except:
            if strict:
                raise
            else:
                # Find the variable
                varname = collect_findVar(varname, f.list())

                dimens = f.dimensions(varname)
                #ndims = len(dimens)
                ndims = f.ndims(varname)

        # ndims is 0 for reals, and 1 for f.ex. t_array
        if ndims < 2:
            # Just read from file
            if varname != 't_array':
                data = f.read(varname)
            elif (varname == 't_array') and (tind is None):
                data = f.read(varname)
            elif (varname == 't_array') and (tind is not None):
                data = f.read(varname, ranges=[tind[0],tind[1]+1])
            return data

        if ndims > 4:
            raise ValueError("ERROR: Too many dimensions")

        mxsub = f.read("MXSUB")
        if mxsub is None:
            raise ValueError("Missing MXSUB variable")
        mysub = f.read("MYSUB")
        mz    = f.read("MZ")
        myg   = f.read("MYG")
        t_array = f.read("t_array")
        if t_array is None:
            nt = 1
            t_array = numpy.zeros(1)
        else:
            nt = len(t_array)
            if tind_auto:
                for i in range(len(self._file_list)):
                    t_array_ = self._getDataFile(i).read("t_array")
                    nt = min(len(t_array_),nt)

        if info:
            print("mxsub = %d mysub = %d mz = %d\n" % (mxsub, mysub, mz))

        # Get the version of BOUT++ (should be > 0.6 for NetCDF anyway)
        try:
            version = f["BOUT_VERSION"]
        except KeyError:
            print("BOUT++ version : Pre-0.2")
            version = 0
        if version < 3.5:
            # Remove extra point
            nz = mz-1
        else:
            nz = mz

        # Fallback to sensible (?) defaults
        try:
            nxpe = f["NXPE"]
        except KeyError:
            nxpe = 1
            print("NXPE not found, setting to {}".format(nxpe))
        try:
            mxg  = f["MXG"]
        except KeyError:
            mxg = 0
            print("MXG not found, setting to {}".format(mxg))
        try:
            nype = f["NYPE"]
        except KeyError:
            nype = nfiles
            print("NYPE not found, setting to {}".format(nype))

        npe = nxpe * nype
        if info:
            print("nxpe = %d, nype = %d, npe = %d\n" % (nxpe, nype, npe))
            if npe < nfiles:
                print("WARNING: More files than expected (" + str(npe) + ")")
            elif npe > nfiles:
                print("WARNING: Some files missing. Expected " + str(npe))

        if xguards:
            nx = nxpe * mxsub + 2*mxg
        else:
            nx = nxpe * mxsub

        if yguards:
            ny = mysub * nype + 2*myg
        else:
            ny = mysub * nype

        # Check ranges

        def check_range(r, low, up, name="range"):
            r2 = r
            if r is not None:
                try:
                    n = len(r2)
                except:
                    # No len attribute, so probably a single number
                    r2 = [r2,r2]
                if (len(r2) < 1) or (len(r2) > 2):
                    print("WARNING: "+name+" must be [min, max]")
                    r2 = None
                else:
                    if len(r2) == 1:
                        r2 = [r2,r2]
                    if r2[0] < 0 and low >= 0:
                        r2[0]+=(up-low+1)
                    if r2[1] < 0 and low >= 0:
                        r2[1]+=(up-low+1)
                    if r2[0] < low:
                        r2[0] = low
                    if r2[0] > up:
                        r2[0] = up
                    if r2[1] < low:
                        r2[1] = low
                    if r2[1] > up:
                        r2[1] = up
                    if r2[0] > r2[1]:
                        tmp = r2[0]
                        r2[0] = r2[1]
                        r2[1] = tmp
            else:
                r2 = [low, up]
            return r2

        xind = check_range(xind, 0, nx-1, "xind")
        yind = check_range(yind, 0, ny-1, "yind")
        zind = check_range(zind, 0, nz-1, "zind")
        tind = check_range(tind, 0, nt-1, "tind")

        xsize = xind[1] - xind[0] + 1
        ysize = yind[1] - yind[0] + 1
        zsize = zind[1] - zind[0] + 1
        tsize = tind[1] - tind[0] + 1

        # Map between dimension names and output size
        sizes = {'x':xsize, 'y':ysize, 'z':zsize, 't':tsize}

        # Create a list with size of each dimension
        ddims = [sizes[d] for d in dimens]

        # Create the data array
        data = numpy.zeros(ddims)

        for i in range(npe):
            # Get X and Y processor indices
            pe_yind = int(i/nxpe)
            pe_xind = i % nxpe

            inrange = True

            if yguards:
                # Get local ranges
                ymin = yind[0] - pe_yind*mysub
                ymax = yind[1] - pe_yind*mysub

                # Check lower y boundary
                if pe_yind == 0:
                    # Keeping inner boundary
                    if ymax < 0: inrange = False
                    if ymin < 0: ymin = 0
                else:
                    if ymax < myg: inrange = False
                    if ymin < myg: ymin = myg

                # Upper y boundary
                if pe_yind == (nype - 1):
                    # Keeping outer boundary
                    if ymin >= (mysub + 2*myg): inrange = False
                    if ymax > (mysub + 2*myg - 1): ymax = (mysub + 2*myg - 1)
                else:
                    if ymin >= (mysub + myg): inrange = False
                    if ymax >= (mysub + myg): ymax = (mysub+myg-1)

                # Calculate global indices
                ygmin = ymin + pe_yind * mysub
                ygmax = ymax + pe_yind * mysub

            else:
                # Get local ranges
                ymin = yind[0] - pe_yind*mysub + myg
                ymax = yind[1] - pe_yind*mysub + myg

                if (ymin >= (mysub + myg)) or (ymax < myg):
                    inrange = False # Y out of range

                if ymin < myg:
                    ymin = myg
                if ymax >= mysub+myg:
                    ymax = myg + mysub - 1

                # Calculate global indices
                ygmin = ymin + pe_yind * mysub - myg
                ygmax = ymax + pe_yind * mysub - myg

            if xguards:
                # Get local ranges
                xmin = xind[0] - pe_xind*mxsub
                xmax = xind[1] - pe_xind*mxsub

                # Check lower x boundary
                if pe_xind == 0:
                    # Keeping inner boundary
                    if xmax < 0: inrange = False
                    if xmin < 0: xmin = 0
                else:
                    if xmax < mxg: inrange = False
                    if xmin < mxg: xmin = mxg

                # Upper x boundary
                if pe_xind == (nxpe - 1):
                    # Keeping outer boundary
                    if xmin >= (mxsub + 2*mxg): inrange = False
                    if xmax > (mxsub + 2*mxg - 1): xmax = (mxsub + 2*mxg - 1)
                else:
                    if xmin >= (mxsub + mxg): inrange = False
                    if xmax >= (mxsub + mxg): xmax = (mxsub+mxg-1)

                # Calculate global indices
                xgmin = xmin + pe_xind * mxsub
                xgmax = xmax + pe_xind * mxsub

            else:
                # Get local ranges
                xmin = xind[0] - pe_xind*mxsub + mxg
                xmax = xind[1] - pe_xind*mxsub + mxg

                if (xmin >= (mxsub + mxg)) or (xmax < mxg):
                    inrange = False # X out of range

                if xmin < mxg:
                    xmin = mxg
                if xmax >= mxsub+mxg:
                    xmax = mxg + mxsub - 1

                # Calculate global indices
                xgmin = xmin + pe_xind * mxsub - mxg
                xgmax = xmax + pe_xind * mxsub - mxg


            # Number of local values
            nx_loc = xmax - xmin + 1
            ny_loc = ymax - ymin + 1

            if not inrange:
                continue # Don't need this file

            if info:
                sys.stdout.write("\rReading from file " + str(i) + ": [" + \
                                     str(xmin) + "-" + str(xmax) + "][" + \
                                     str(ymin) + "-" + str(ymax) + "] -> [" + \
                                     str(xgmin) + "-" + str(xgmax) + "][" + \
                                     str(ygmin) + "-" + str(ygmax) + "]")

            f = self._getDataFile(i)

            if ndims == 4:
                d = f.read(varname, ranges=[tind[0],tind[1]+1,
                                            xmin, xmax+1,
                                            ymin, ymax+1,
                                            zind[0],zind[1]+1])
                data[:, (xgmin-xind[0]):(xgmin-xind[0]+nx_loc), (ygmin-yind[0]):(ygmin-yind[0]+ny_loc), :] = d
            elif ndims == 3:
                # Could be xyz or txy

                if dimens[2] == 'z': # xyz
                    d = f.read(varname, ranges=[xmin, xmax+1,
                                                ymin, ymax+1,
                                                zind[0],zind[1]+1])
                    data[(xgmin-xind[0]):(xgmin-xind[0]+nx_loc), (ygmin-yind[0]):(ygmin-yind[0]+ny_loc), :] = d
                else: # txy
                    d = f.read(varname, ranges=[tind[0],tind[1]+1,
                                                xmin, xmax+1,
                                                ymin, ymax+1])
                    data[:, (xgmin-xind[0]):(xgmin-xind[0]+nx_loc), (ygmin-yind[0]):(ygmin-yind[0]+ny_loc)] = d
            elif ndims == 2:
                # xy
                d = f.read(varname, ranges=[xmin, xmax+1,
                                            ymin, ymax+1])
                data[(xgmin-xind[0]):(xgmin-xind[0]+nx_loc), (ygmin-yind[0]):(ygmin-yind[0]+ny_loc)] = d

        # Force the precision of arrays of dimension>1
        if ndims>1:
            try:
                data = data.astype(t_array.dtype, copy=False)
            except TypeError:
                data = data.astype(t_array.dtype)

        # Finished looping over all files
        if info:
            sys.stdout.write("\n")
        return data

    def __iter__(self):
        """
        Iterate through all keys, starting with "options"
        then going through all variables for _caching_collect
        """
        for k in self.varNames:
            yield k
            
    def __str__(self, indent=""):
        """
        Print a pretty version of the tree
        """
        text = ""
        for k in self.varNames:
            text += indent+k+"\n"
        
        return text


def BoutData(path=".", prefix="BOUT.dmp", caching=False, **kwargs):
    """
    Returns a dictionary, containing the contents of a BOUT++ output directory.
    Does not allow writing, only reading of data.  By default there is no
    cache, so each time a variable is read it is collected; if caching is set
    to True variables are stored once they are read.
    
    Example
    -------

    d = BoutData(".")  # Current directory
    
    d.keys()     # List all valid keys

    print d["options"]  # Prints tree of options

    d["options"]["nout"]   # Value of nout in BOUT.inp file
    
    print d["outputs"]    # Print available outputs

    d["outputs"]["ne"] # Read "ne" from data files
    
    d = BoutData(".", prefix="BOUT.dmp", caching=True) # Turn on caching

    Options
    -------
    prefix - sets the prefix for data files (default "BOUT.dmp")

    caching - switches on caching of data, so it is only read into memory when
              first accessed (default False) If caching is set to a number, it
              gives the maximum size of the cache in GB, after which entries
              will be discarded in first-in-first-out order to prevent the
              cache getting too big.  If the variable being returned is bigger
              than the maximum cache size, then the variable will be returned
              without being added to the cache, and the rest of the cache will
              be left.
    
    **kwargs - keyword arguments that are passed through to collect()
    """
    
    data = {} # Map for the result
    
    data["path"] = path
    
    # Options from BOUT.inp file
    data["options"] = BoutOptionsFile(os.path.join(path, "BOUT.inp"), name="options")
    
    # Output from .dmp.* files
    data["outputs"] = BoutOutputs(path, prefix=prefix, caching=caching, **kwargs)
    
    return data
