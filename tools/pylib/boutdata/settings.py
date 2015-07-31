# Parse BOUT.inp settings file

def get(filename, name, section=None):
    """
    Find and return a single value from a BOUT.inp settings file
    
    Inputs
    ------
    
    filename     A string containing a BOUT.inp file name
    name         The name of the setting
    section      (optional) The section to look in. 
                 By default this is None (global)

    Note that names and sections are case insensitive
    
    Returns
    -------
    
    value of the setting as a string

    If not found, raises a ValueError
    
    Example
    -------
    
    nout = settings.get("BOUT.inp", "nout")
    
    compress = settings.get("BOUT.inp", "compress", section="highbeta")
    
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

