"""Import an entire BOUT++ DataFile into memory

"""

from boututils.datafile import DataFile


def file_import(name):
    """Read all variables from file into a dictionary

    Parameters
    ----------
    name : str
        Name of file to read

    Returns
    -------
    dict
        Dictionary containing all the variables in the file
    """
    f = DataFile(name)   # Open file
    varlist = f.list() # Get list of all variables in file
    data = {}          # Create empty dictionary
    for v in varlist:
        data[v] = f.read(v)
    f.close()
    return data
