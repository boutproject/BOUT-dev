; Tests if a given name is a member of the structure

FUNCTION in_struct, str, name
  ON_ERROR, 2
  
  RETURN, in_list( TAG_NAMES(str), name )
END
