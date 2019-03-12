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
  ; 1.1.3 - Handle case when break of contour is very close to X-point.
  ;         Add checkbox to write metrics for orthogonal coordinates (i.e.
  ;         metrics with integrated shear I=0); include attribute in output
  ;         file labelling the 'coordinates_type', either 'field_aligned' or
  ;         'orthogonal'.
  ;         Various small fixes to make non-orthogonal grid generation more
  ;         robust.
  ;         Better handling of closed contours
  ;         Make transitions of non-orthogonal coordinates between X-point and
  ;         wall more flexible: instead of forcing coordinates to be orthogonal
  ;         half-way through the poloidal region, have separate weights for the
  ;         boundary vectors at either end to make the transition between the
  ;         two continuous (with dominant weighting of the orthogonal direction
  ;         in the middle of the region).
  ;         Enable 'Detailed settings' dialog for non-orthogonal grids, add
  ;         setting to change the exponent of the power-law decrease of the
  ;         weight of non-orthogonal vectors.
  ; 1.1.4 - For non-orthogonal, always refine of starting locations (previously
  ;         was only done for grids with only closed field lines) - improves
  ;         robustness especially around X-points.
  ;         Pass 'simple' setting through when restarting grid generation.
  ;         Option to save y-boundary guard cells (defaults to 0 for backward
  ;         compatibility with versions of BOUT++ before 4.3).
  
  major_version = 1
  minor_version = 1
  patch_number = 4

  RETURN, LONG([major_version, minor_version, patch_number])

END
