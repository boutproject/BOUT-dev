"""Parse BOUT.inp settings file

"""


def get(filename, name, section=None):
    """Find and return a single value from a BOUT.inp settings file

    .. deprecated::3.0
            `settings.get` has been replaced with
            `boututils.options.BoutOptions`

    Parameters
    ----------
    filename : str
        Name of the settings file
    name : str
        The name of the setting
    section : str, optional
        The section to look in (default: the global section)

    Note that names and sections are case insensitive

    Returns
    -------
    str
        Value of the setting. If not found, raises a ValueError

    Examples
    --------

    >>> settings.get("BOUT.inp", "nout")
    '100'

    >>> settings.get("BOUT.inp", "compress", section="highbeta")
    'true'

    """
    with open(filename, "rt") as f:
        if section is not None:
            # First find the section
            found = False
            for line in f:
                # Strip spaces from left
                line = line.lstrip(' \t\n\r')
                if len(line) < 1:
                    continue  # Empty line
                    
                # if line starts with '[' then this is a section
                if line[0] == '[':
                    # Split on ']'
                    head, _ = line[1:].split(']', 1)
                    # head is now the section name
                    if head.lower() == section.lower():
                        found = True
                        break
            if not found:
                raise ValueError("Section '%s' not found" % (section))
        
        # Now in the correct section
        
        for line in f:
            # Strip spaces from left
            line = line.lstrip(' \t\n\r')
            if len(line) < 1:
                continue  # Empty line
                
            # if line starts with '[' then this is a section
            if line[0] == '[':
                raise ValueError("Name '%s' not found in section '%s'" % (name,section))
            # Check if this line contains an '='
            if '=' in line:
                # Check if contains comment
                comment = ''
                if '#' in line:
                    line, comment = line.split('#', 1)
                # Split on '='
                key, value = line.split('=',1)
                # Strip whitespace
                key   = key.strip(' \t\n\r')
                value = value.strip(' \t\n\r')
                
                # Strip out quotes if present
                if value[0] == '"' or value[0] == "'": 
                    value = value[1:]
                if value[-1] == '"' or value[-1] == "'":
                    value = value[:-1]
                
                #print("'%s' = '%s'" % (key, value))
                if key.lower() == name.lower(): # Case insensitive
                    return value

