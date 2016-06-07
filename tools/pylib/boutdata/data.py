# Provides a class BoutData which makes access to code
# inputs and outputs easier. Creates a tree of maps,
# inspired by approach used in OMFIT
#
# 

import os
import glob

from boutdata import collect

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
            text += " |- " + self._sections[s].__str__(indent+" |  ")
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
                        sectionname = line[0,scorepos]
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
    Emulates a map class, represents the contents of a BOUT++
    dmp files. Does not allow writing, only reading of data.
    Currently there is no cache, so each time a variable
    is read it is collected.
    
    Example
    -------

    d = BoutOutputs(".")  # Current directory
    
    d.keys()     # List all valid keys
    
    d["ne"] # Read "ne" from data files
    
    """
    def __init__(self, path=".", prefix="BOUT.dmp"):
        """
        Initialise BoutOutputs object
        """
        self._path = path
        self._prefix = prefix
        
        # Label for this data
        self.label = path

        # Check that the path contains some data
        file_list = glob.glob(os.path.join(path, prefix+"*.nc"))
        if len(file_list) == 0:
            raise ValueError("ERROR: No data files found")
        
        # Available variables
        self.varNames = []
        
        with DataFile(file_list[0]) as f:
            # Get variable names
            self.varNames = f.keys()
        
    def keys(self):
        """
        Return a list of available variable names
        """
        return self.varNames
        
    def __len__(self):
        return len(self.varNames)
            
    def __getitem__(self, name):
        """
        Reads a variable using collect.
        
        """

        # Collect the data from the repository
        data = collect(name, path=self._path, prefix=self._prefix)
        return data

    def __iter__(self):
        """
        Iterate through all keys, starting with "options"
        then going through all variables for collect
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
    

def BoutData(path=".", prefix="BOUT.dmp"):
    """
    Returns a dictionary, containing the contents of a BOUT++
    output directory. Does not allow writing, only reading of data.
    Currently there is no cache, so each time a variable
    is read it is collected.
    
    Example
    -------

    d = BoutData(".")  # Current directory
    
    d.keys()     # List all valid keys

    print d["options"]  # Prints tree of options

    d["options"]["nout"]   # Value of nout in BOUT.inp file
    
    print d["outputs"]    # Print available outputs

    d["outputs"]["ne"] # Read "ne" from data files
    
    """
    
    data = {} # Map for the result
    
    data["path"] = path
    
    # Options from BOUT.inp file
    data["options"] = BoutOptionsFile(os.path.join(path, "BOUT.inp"), name="options")
    
    # Output from .dmp.* files
    data["outputs"] = BoutOutputs(path)
    
    return data
