FUNCTION hypnotoad_version
  ; This function defines, and returns, the version number of Hypnotoad
  ; - The major version should increase for substantial changes that change the
  ;   format of the output grid files
  ; - The minor version should increase when new features are added
  ; - The patch number will now increase when bugs are fixed or minor tweaks are made
  ;
  ; Version history:
  ; 1.0.0 - original version of hypnotoad
  ; 1.1.0 - non-orthogonal grid generation added
  ; 1.1.1 - Hypnotoad version number added here, and now saved to grid files
  ; 1.1.2 - Fixed bug in calculation of qloop. Should be only in closed regions
  
  major_version = 1
  minor_version = 1
  patch_number = 2

  RETURN, LONG([major_version, minor_version, patch_number])

END
