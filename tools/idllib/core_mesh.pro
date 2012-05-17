;
; Extracts only the core part of data arrays or mesh
;

FUNCTION extract_2d, var, xr, yr
  nx = N_ELEMENTS(xr)
  ny = N_ELEMENTS(yr)
  
  result = DBLARR(nx, ny)
  
  FOR i=0, nx-1 DO BEGIN
    result[i,*] = var[xr[i], yr]
  ENDFOR
  RETURN, result
END

FUNCTION in_list, val, list
  FOR i=0, N_ELEMENTS(list)-1 DO BEGIN
    IF val EQ list[i] THEN RETURN, 1
  ENDFOR
  RETURN, 0
END

FUNCTION core_mesh, var, mesh

  IF SIZE(var, /type) EQ 8 THEN BEGIN
    ; Given a mesh as input. Process and return another mesh
    
    xrange = range(0, MIN([var.ixseps1, var.ixseps2])-1)
    yrange = [range(var.jyseps1_1+1, var.jyseps2_1), range(var.jyseps1_2+1, var.jyseps2_2)]
    nx = N_ELEMENTS(xrange)
    ny = N_ELEMENTS(yrange)
    
    ; Create a new mesh
    g = {nx:nx, ny:ny, $
         jyseps1_1:-1, jyseps2_1:FIX(ny/2), jyseps1_2:FIX(ny/2), jyseps2_2:ny-1, $
         ixseps1:nx, ixseps2:nx}
    
    
    ; Get a list of variables
    tags = TAG_NAMES(var)
    done = TAG_NAMES(g) ; Tags already set
    nt = N_ELEMENTS(tags)
    
    FOR i=0, nt-1 DO BEGIN
      IF in_list(tags[i], done) THEN CONTINUE ; Skip these
      
      ndim = SIZE(var.(i), /n_dim) ; Number of dimensions
      IF ndim EQ 2 THEN BEGIN
        g = CREATE_STRUCT(tags[i], extract_2d(var.(i), xrange, yrange), g)
      ENDIF ELSE IF ndim EQ 1 THEN BEGIN
        IF N_ELEMENTS(var.(i)) EQ var.nx THEN BEGIN
          g = CREATE_STRUCT(tags[i], (var.(i))[xrange], g)
        ENDIF ELSE g = CREATE_STRUCT(tags[i], var.(i), g)
      ENDIF ELSE BEGIN
        ; Just add everything else to g
        g = CREATE_STRUCT(tags[i], var.(i), g)
      ENDELSE
    ENDFOR
    RETURN, g
  ENDIF
  
  ; Should be an array
  ndim = SIZE(var, /N_DIM) ; Number of dimensions
  
  IF SIZE(mesh, /type) NE 8 THEN STOP

  xrange = range(0, MIN([mesh.ixseps1, mesh.ixseps2])-1)
  yrange = [range(mesh.jyseps1_1+1, mesh.jyseps2_1), range(mesh.jyseps1_2+1, mesh.jyseps2_2)]

  nx = N_ELEMENTS(xrange)
  ny = N_ELEMENTS(yrange)

  CASE ndim OF
    4: BEGIN
      nz = (SIZE(var, /dim))[2]
      nt = (SIZE(var, /dim))[3]
      result = DBLARR(nx, ny, nz, nt)
      FOR z=0, nz-1 DO BEGIN
        FOR t=0, nt-1 DO BEGIN
          result[*,*,z] = extract_2d(REFORM(var[*,*,z]), xrange, yrange)
        ENDFOR
      ENDFOR
      RETURN, result
    END
    3: BEGIN
      nz = (SIZE(var, /dim))[2]
      result = DBLARR(nx, ny, nz)
      FOR z=0, nz-1 DO BEGIN
        result[*,*,z] = extract_2d(REFORM(var[*,*,z]), xrange, yrange)
      ENDFOR
      RETURN, result
    END
    2: BEGIN
      RETURN, extract_2d(var, xrange, yrange)
    END
    ELSE: BEGIN
      PRINT, "Not sure how to handle variable. dimensions =", SIZE(var, /dim)
      RETURN, var
    END
  ENDCASE
END

