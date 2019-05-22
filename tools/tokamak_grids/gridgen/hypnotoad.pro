;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                  HYPNO-TOAD grid generator
;
; 
; This is a graphical interface to some grid generation routines.
; Aims to allow tokamak grids to be easily generated from
; a variety of input sources.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO plot_region, R, Z, ymin, ymax, _extra=_extra
  s = SIZE(R, /dimen)
  nx = s[0]
  ny = s[1]
  
  FOR j=0, nx-1 DO BEGIN
    OPLOT, R[j,ymin:ymax], Z[j,ymin:ymax], _extra=_extra
  ENDFOR
  
  FOR j=ymin, ymax DO BEGIN
    OPLOT, R[*,j], Z[*,j], _extra=_extra
  ENDFOR
END

PRO oplot_mesh, rz_mesh, flux_mesh
  ; Plot X-points and separatrices
  FOR i=0, flux_mesh.critical.n_xpoint-1 DO BEGIN
    ; plot the separatrix contour
    CONTOUR, rz_mesh.psi, rz_mesh.R, rz_mesh.Z, levels=[flux_mesh.critical.xpt_f[i]], c_colors=2, /overplot
    oplot, [INTERPOLATE(rz_mesh.R, flux_mesh.critical.xpt_ri[i], /DOUBLE)], [INTERPOLATE(rz_mesh.Z, flux_mesh.critical.xpt_zi[i], /DOUBLE)], psym=7, color=2
  ENDFOR

  ; Plot O-points
  FOR i=0, flux_mesh.critical.n_opoint-1 DO BEGIN
    oplot, [INTERPOLATE(rz_mesh.R, flux_mesh.critical.opt_ri[i], /DOUBLE)], [INTERPOLATE(rz_mesh.Z, flux_mesh.critical.opt_zi[i], /DOUBLE)], psym=7, color=3
  ENDFOR
  
  ypos = 0
  FOR i=0, N_ELEMENTS(flux_mesh.npol)-1 DO BEGIN
    plot_region, flux_mesh.rxy, flux_mesh.zxy, ypos, ypos+flux_mesh.npol[i]+flux_mesh.n_y_boundary_guards[i]-1, color=i+1
    ypos = ypos + flux_mesh.npol[i]+flux_mesh.n_y_boundary_guards[i]
  ENDFOR
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Event handling procedures

; For detailed settings popup
PRO popup_event, event
  ; Get the UVALUE
  widget_control, event.id, get_uvalue=uvalue
  
  ; Retrieve a copy of information stored in tlb
  widget_control, event.top, get_uvalue=info

  widget_control, info.top, get_uvalue=base_info
  
  IF N_ELEMENTS(uvalue) EQ 0 THEN RETURN ; Undefined
  
  CASE 1 OF
    uvalue EQ 'mesh' OR uvalue EQ 'mesh2': BEGIN
      IF base_info.rz_grid_valid EQ 0 THEN BEGIN
        PRINT, "ERROR: No valid equilibrium data. Read from file first"
        a = DIALOG_MESSAGE("No valid equilibrium data. Read from file first", /error)
        RETURN
      ENDIF
      
      boundary = TRANSPOSE([[(*base_info.rz_grid).rlim], [(*base_info.rz_grid).zlim]])
      
      ; retrieve the values from the entry fields
      
      nnrad = N_ELEMENTS(info.nrad_field)
      nrad = LONARR(nnrad)
      FOR i=0, nnrad-1 DO BEGIN
        widget_control, info.nrad_field[i], get_value=nr
        nrad[i] = nr
      ENDFOR

      ninpsi = N_ELEMENTS(info.in_psi_field)
      psi_inner = DBLARR(ninpsi)
      FOR i=0, ninpsi-1 DO BEGIN
        widget_control, info.in_psi_field[i], get_value=inp
        psi_inner[i] = inp
      ENDFOR

      noutpsi = N_ELEMENTS(info.out_psi_field)
      psi_outer = DBLARR(noutpsi)
      FOR i=0, noutpsi-1 DO BEGIN
        widget_control, info.out_psi_field[i], get_value=inp
        psi_outer[i] = inp
      ENDFOR

      nnpol = N_ELEMENTS(info.npol_field)
      npol = LONARR(nnpol)
      FOR i=0, nnpol-1 DO BEGIN
        widget_control, info.npol_field[i], get_value=np
        npol[i] = np
      ENDFOR

      widget_control, base_info.rad_peak_field, get_value=rad_peak
      
      widget_control, base_info.xpt_dist_field, get_value=xpt_mul

      widget_control, info.nonorthogonal_weight_decay_power, get_value=nonorthogonal_weight_decay_power

      widget_control, base_info.y_boundary_guards_field, get_value=y_boundary_guards

      PRINT, "HYP: psi_inner =", psi_inner

      settings = {nrad:nrad, npol:npol, psi_inner:psi_inner, psi_outer:psi_outer, $
                  rad_peaking:rad_peak, $
                  nonorthogonal_weight_decay_power:nonorthogonal_weight_decay_power}
      
      WIDGET_CONTROL, base_info.status, set_value="Generating mesh ..."
      
      ; Delete the window, as number of fields may change
      WIDGET_CONTROL, event.top, /destroy

      ; Check if a simplified boundary should be used
      IF base_info.simple_bndry THEN BEGIN
        ; Simplify the boundary to a square box
        boundary = TRANSPOSE([ [MIN(boundary[0,*]), MAX(boundary[0,*]), $
                                MAX(boundary[0,*]), MIN(boundary[0,*])], $
                               [MIN(boundary[1,*]), MIN(boundary[1,*]), $
                                MAX(boundary[1,*]), MAX(boundary[1,*])] ])
      ENDIF
    END
  ENDCASE
      
  CASE uvalue OF
    'mesh': BEGIN
      ; Orthogonal mesh button was pushed

      ; Create the mesh
      mesh = create_grid((*(base_info.rz_grid)).psi, (*(base_info.rz_grid)).r, (*(base_info.rz_grid)).z, $
                         settings, $
                         boundary=boundary, strict=base_info.strict_bndry, $
                         single_rad_grid=base_info.single_rad_grid, $
                         critical=(*(base_info.rz_grid)).critical, $
                         fast=base_info.fast, xpt_mul=xpt_mul, /simple, $
                         y_boundary_guards=y_boundary_guards)
    END
    'mesh2': BEGIN
      ; Non-orthogonal mesh button was pushed

      ; Create the mesh
      mesh = create_nonorthogonal((*(base_info.rz_grid)).psi, (*(base_info.rz_grid)).r, $
                         (*(base_info.rz_grid)).z, settings, $
                         boundary=boundary, strict=base_info.strict_bndry, $
                         /nrad_flexible, $
                         single_rad_grid=base_info.single_rad_grid, $
                         critical=(*(base_info.rz_grid)).critical, $
                         xpt_only=base_info.xpt_only, /simple, $
                         y_boundary_guards=y_boundary_guards)
    END
  ENDCASE
      
  CASE 1 OF
    uvalue EQ 'mesh' OR uvalue EQ 'mesh2': BEGIN
      IF mesh.error EQ 0 THEN BEGIN
        PRINT, "Successfully generated mesh"
        WIDGET_CONTROL, base_info.status, set_value="Successfully generated mesh. All glory to the Hypnotoad!"
        oplot_mesh, (*base_info.rz_grid), mesh
        
        base_info.flux_mesh_valid = 1
        base_info.flux_mesh = PTR_NEW(mesh)
        widget_control, info.top, set_UVALUE=base_info
      ENDIF ELSE BEGIN
        a = DIALOG_MESSAGE("Could not generate mesh", /error)
        WIDGET_CONTROL, base_info.status, set_value="  *** FAILED to generate mesh ***"
      ENDELSE
    END
  ENDCASE
END

; For the main window
PRO event_handler, event
  ; Get the UVALUE
  widget_control, event.id, get_uvalue=uvalue
  
  ; Get info stored in base
  widget_control, event.top, get_uvalue=info
  
  IF N_ELEMENTS(uvalue) EQ 0 THEN RETURN ; Undefined
  
  CASE uvalue OF
    'aandg': BEGIN
      PRINT, "Open G-eqdsk (neqdsk) file"
      filename = DIALOG_PICKFILE(dialog_parent=event.top, file="neqdsk", /read, path=info.path, get_path=newpath)
      info.path=newpath
      IF STRLEN(filename) EQ 0 THEN BEGIN
        WIDGET_CONTROL, info.status, set_value="   *** Cancelled open file ***"
        RETURN ;BREAK
      ENDIF
      PRINT, "Trying to read file "+filename
      g = read_neqdsk(filename)
      
      IF SIZE(g, /TYPE) EQ 8 THEN BEGIN
        ; Got a structure
        PRINT, "Successfully read equilibrium"
        WIDGET_CONTROL, info.status, set_value="Successfully read "+filename
        
        ; Analyse the equilibrium
        critical = analyse_equil(g.psi, REFORM(g.r[*,0]), REFORM(g.z[0,*]))

        IF critical.n_xpoint EQ 0 THEN BEGIN
          ; Analyse_equil will make a guess for where the separatrix
          ; is. Reset to the value in the G file
          critical.xpt_f[0] = g.sibdry
          PRINT, "WARNING: No X-points. Resetting outer psi to "+STR(g.sibdry)
        ENDIF

        ; Extract needed data from g-file struct
        
        rz_grid = {nr:g.nx, nz:g.ny, $  ; Number of grid points
                   r:REFORM(g.r[*,0]), z:REFORM(g.z[0,*]), $  ; R and Z as 1D arrays
                   simagx:g.simagx, sibdry:g.sibdry, $ ; Range of psi
                   psi:g.psi, $  ; Poloidal flux in Weber/rad on grid points
                   npsigrid:(FINDGEN(N_ELEMENTS(g.pres))/(N_ELEMENTS(g.pres)-1)), $ ; Normalised psi grid for fpol, pres and qpsi
                   fpol:g.fpol, $ ; Poloidal current function on uniform flux grid
                   pres:g.pres, $ ; Plasma pressure in nt/m^2 on uniform flux grid
                   qpsi:g.qpsi, $ ; q values on uniform flux grid
                   nlim:g.nlim, rlim:g.xlim, zlim:g.ylim, $ ; Wall boundary
                   ;nlim:g.nbdry, rlim:g.rbdry, zlim:g.zbdry, $
                   critical:critical} ; Critical point structure
        
        
        IF info.rz_grid_valid GT 0 THEN BEGIN
          ; Need to free existing data
          PTR_FREE, info.rz_grid
        ENDIF

        ; Put pointer to data into info struct
        info.rz_grid = PTR_NEW(rz_grid)
        info.rz_grid_valid = 1
        
        ; Plot the equilibrium
        plot_rz_equil, rz_grid
        
        ; Set info to new values
        widget_control, event.top, set_UVALUE=info
        
        IF rz_grid.nlim LT 3 THEN BEGIN
          PRINT, "WARNING: No boundary found!"
        ENDIF
      ENDIF ELSE BEGIN
        ; Couldn't read data
        PRINT, "ERROR: Failed to read grid file"
        WIDGET_CONTROL, info.status, set_value="   *** Failed to read G-EQDSK file "+filename+" ***"
      ENDELSE
    END
    'restorerz': BEGIN
      ; Restore a file containing rz_grid
      filename = DIALOG_PICKFILE(dialog_parent=event.top, /read, path=info.path, get_path=newpath)
      info.path = newpath

      CATCH, err
      IF err NE 0 THEN BEGIN
        PRINT, "ERROR: Failed to restore state"
        WIDGET_CONTROL, info.status, set_value="   *** Failed to restore file "+filename+" ***"
        PRINT, 'Error message: ', !ERROR_STATE.MSG
      ENDIF ELSE BEGIN
        RESTORE, filename
      ENDELSE
      CATCH, /cancel
      
      IF SIZE(rz_grid, /type) NE 8 THEN BEGIN
        PRINT, "Error: File does not contain a variable called 'rz_grid'"
        WIDGET_CONTROL, info.status, set_value="   *** Failed to restore file "+filename+" ***"
      ENDIF ELSE BEGIN
      
        IF info.rz_grid_valid GT 0 THEN BEGIN
          ; Need to free existing data
          PTR_FREE, info.rz_grid
        ENDIF
        
        ; Put pointer to data into info struct
        info.rz_grid = PTR_NEW(rz_grid)
        info.rz_grid_valid = 1
        
        ; Plot the equilibrium
        plot_rz_equil, rz_grid
        
        ; Set info to new values
        widget_control, event.top, set_UVALUE=info
        
        IF rz_grid.nlim LT 3 THEN BEGIN
          PRINT, "WARNING: No boundary found!"
        ENDIF

        WIDGET_CONTROL, info.status, set_value="Restored R-Z data from "+filename
      ENDELSE
    END
    'bndry': BEGIN
      IF info.rz_grid_valid NE 1 THEN BEGIN
        PRINT, "ERROR: Need to read an equilibrium first"
        WIDGET_CONTROL, info.status, set_value="   *** Need to read equilibrium first ***"
        BREAK
      ENDIF
      PRINT, "Read boundary from G-eqdsk (neqdsk) file"
      filename = DIALOG_PICKFILE(dialog_parent=event.top, file="neqdsk", /read, path=info.path, get_path=newpath)
      info.path = newpath
      IF STRLEN(filename) EQ 0 THEN BEGIN
        WIDGET_CONTROL, info.status, set_value="   *** Cancelled open file ***"
        BREAK
      ENDIF
      PRINT, "Trying to read file "+filename
      g = read_neqdsk(filename)
      IF SIZE(g, /TYPE) EQ 8 THEN BEGIN
        ; Got a structure
        PRINT, "Successfully read equilibrium"
        WIDGET_CONTROL, info.status, set_value="Successfully read "+filename
        
        rz_grid = {nr:(*info.rz_grid).nr, nz:(*info.rz_grid).nz, $
                   r:(*info.rz_grid).r, z:(*info.rz_grid).z, $
                   simagx:(*info.rz_grid).simagx, sibdry:(*info.rz_grid).sibdry, $
                   psi:(*info.rz_grid).psi, $
                   npsigrid:(*info.rz_grid).npsigrid, $
                   fpol:(*info.rz_grid).fpol, $
                   pres:(*info.rz_grid).pres, $
                   qpsi:(*info.rz_grid).qpsi, $
                   nlim:g.nlim, rlim:g.xlim, zlim:g.ylim, $ ; Wall boundary
                   critical:(*info.rz_grid).critical}
        IF rz_grid.nlim LT 3 THEN BEGIN
          PRINT, "WARNING: No boundary found!"
        ENDIF ELSE BEGIN
          PTR_FREE, info.rz_grid
          info.rz_grid = PTR_NEW(rz_grid)
          widget_control, event.top, set_UVALUE=info
          
          ; Plot the equilibrium
          plot_rz_equil, rz_grid
        ENDELSE
      ENDIF ELSE BEGIN
        ; Couldn't read data
        PRINT, "ERROR: Failed to read grid file"
        WIDGET_CONTROL, info.status, set_value="   *** Failed to read G-EQDSK file "+filename+" ***"
      ENDELSE
    END
    'mesh': BEGIN
      ; Create a mesh
      IF info.rz_grid_valid EQ 0 THEN BEGIN
        PRINT, "ERROR: No valid equilibrium data. Read from file first"
        a = DIALOG_MESSAGE("No valid equilibrium data. Read from file first", /error)
        RETURN
      ENDIF
      
      boundary = TRANSPOSE([[(*info.rz_grid).rlim], [(*info.rz_grid).zlim]])
      
      IF info.detail_set THEN BEGIN
        settings = {dummy:0}
      ENDIF ELSE BEGIN
        ; Get settings
        widget_control, info.nrad_field, get_value=nrad
        widget_control, info.npol_field, get_value=npol

        widget_control, info.psi_inner_field, get_value=psi_inner
        widget_control, info.psi_outer_field, get_value=psi_outer

        widget_control, info.rad_peak_field, get_value=rad_peak

        widget_control, info.parweight_field, get_value=parweight

        settings = {nrad:nrad, npol:npol, psi_inner:psi_inner, psi_outer:psi_outer, rad_peaking:rad_peak, parweight:parweight}
      ENDELSE
      
      widget_control, info.xpt_dist_field, get_value=xpt_mul
      PRINT, "xpt_mul = ", xpt_mul

      widget_control, info.y_boundary_guards_field, get_value=y_boundary_guards

      ; Check if a simplified boundary should be used
      IF info.simple_bndry THEN BEGIN
        ; Simplify the boundary to a square box
        boundary = TRANSPOSE([ [MIN(boundary[0,*]), MAX(boundary[0,*]), $
                                MAX(boundary[0,*]), MIN(boundary[0,*])], $
                               [MIN(boundary[1,*]), MIN(boundary[1,*]), $
                                MAX(boundary[1,*]), MAX(boundary[1,*])] ])
      ENDIF
        
      WIDGET_CONTROL, info.status, set_value="Generating mesh ..."
      
      fpsi = DBLARR(2, N_ELEMENTS((*(info.rz_grid)).fpol))
      fpsi[0,*] = (*(info.rz_grid)).simagx + (*(info.rz_grid)).npsigrid * ( (*(info.rz_grid)).sibdry - (*(info.rz_grid)).simagx )
      fpsi[1,*] = (*(info.rz_grid)).fpol

      mesh = create_grid((*(info.rz_grid)).psi, (*(info.rz_grid)).r, (*(info.rz_grid)).z, settings, $
                         boundary=boundary, strict=info.strict_bndry, $
                         /nrad_flexible, $
                         single_rad_grid=info.single_rad_grid, $
                         critical=(*(info.rz_grid)).critical, $
                         fast=info.fast, xpt_mul=xpt_mul, $
                         fpsi = fpsi, /simple, y_boundary_guards=y_boundary_guards)
      IF mesh.error EQ 0 THEN BEGIN
        PRINT, "Successfully generated mesh"
        WIDGET_CONTROL, info.status, set_value="Successfully generated mesh. All glory to the Hypnotoad!"
        oplot_mesh, *info.rz_grid, mesh
        
        info.flux_mesh_valid = 1
        info.flux_mesh = PTR_NEW(mesh)
        widget_control, event.top, set_UVALUE=info
      ENDIF ELSE BEGIN
        a = DIALOG_MESSAGE("Could not generate mesh", /error, dialog_parent=info.draw)
        WIDGET_CONTROL, info.status, set_value="  *** FAILED to generate mesh ***"
      ENDELSE
    END
    'mesh2': BEGIN
      ; Create a non-orthogonal mesh
      IF info.rz_grid_valid EQ 0 THEN BEGIN
        PRINT, "ERROR: No valid equilibrium data. Read from file first"
        a = DIALOG_MESSAGE("No valid equilibrium data. Read from file first", /error)
        RETURN
      ENDIF

      boundary = TRANSPOSE([[(*info.rz_grid).rlim], [(*info.rz_grid).zlim]])
      
      IF info.detail_set THEN BEGIN
        settings = {dummy:0}
      ENDIF ELSE BEGIN
        ; Get settings
        widget_control, info.nrad_field, get_value=nrad
        widget_control, info.npol_field, get_value=npol

        widget_control, info.psi_inner_field, get_value=psi_inner
        widget_control, info.psi_outer_field, get_value=psi_outer

        widget_control, info.rad_peak_field, get_value=rad_peak

        widget_control, info.parweight_field, get_value=parweight

        settings = {nrad:nrad, npol:npol, psi_inner:psi_inner, psi_outer:psi_outer, $
                    rad_peaking:rad_peak, parweight:parweight, $
                    nonorthogonal_weight_decay_power:info.nonorthogonal_weight_decay_power}
      ENDELSE
      
      widget_control, info.y_boundary_guards_field, get_value=y_boundary_guards

      ; Check if a simplified boundary should be used
      IF info.simple_bndry THEN BEGIN
        ; Simplify the boundary to a square box
        boundary = TRANSPOSE([ [MIN(boundary[0,*]), MAX(boundary[0,*]), $
                                MAX(boundary[0,*]), MIN(boundary[0,*])], $
                               [MIN(boundary[1,*]), MIN(boundary[1,*]), $
                                MAX(boundary[1,*]), MAX(boundary[1,*])] ])
      ENDIF
      
      WIDGET_CONTROL, info.status, set_value="Generating non-orthogonal mesh ..."
      
      mesh = create_nonorthogonal((*(info.rz_grid)).psi, (*(info.rz_grid)).r, (*(info.rz_grid)).z, settings, $
                         boundary=boundary, strict=info.strict_bndry, $
                         /nrad_flexible, $
                         single_rad_grid=info.single_rad_grid, $
                         critical=(*(info.rz_grid)).critical, xpt_only=info.xpt_only, /simple, $
                         y_boundary_guards=y_boundary_guards)
      IF mesh.error EQ 0 THEN BEGIN
        PRINT, "Successfully generated non-orthogonal mesh"
        WIDGET_CONTROL, info.status, set_value="Successfully generated mesh. All glory to the Hypnotoad!"
        oplot_mesh, *info.rz_grid, mesh
        
        info.flux_mesh_valid = 1
        info.flux_mesh = PTR_NEW(mesh)
        widget_control, event.top, set_UVALUE=info
      ENDIF ELSE BEGIN
        a = DIALOG_MESSAGE("Could not generate non-orthogonal mesh", /error, dialog_parent=info.draw)
        WIDGET_CONTROL, info.status, set_value="  *** FAILED to generate non-orthogonal mesh ***"
      ENDELSE
    END
    'process': BEGIN
      ; Process mesh to produce output
      PRINT, "Write output file"
      filename = DIALOG_PICKFILE(dialog_parent=event.top, file="bout.grd.nc", $
                                 /write, /overwrite_prompt, path=info.path, get_path=newpath)
      info.path=newpath
      
      IF STRLEN(filename) EQ 0 THEN BEGIN
        WIDGET_CONTROL, info.status, set_value="   *** Cancelled process mesh ***"
        BREAK
      ENDIF

      IF info.rz_grid_valid AND info.flux_mesh_valid THEN BEGIN
        ; Get settings
        settings = {calcp:info.calcp, calcbt:info.calcbt, $
                    calchthe:info.calchthe, calcjpar:info.calcjpar, $
                    orthogonal_coordinates_output:info.orthogonal_coordinates_output, $
                    y_boundary_guards:(*info.flux_mesh).y_boundary_guards}
        
        process_grid, *(info.rz_grid), *(info.flux_mesh), $
                      output=filename, poorquality=poorquality, /gui, parent=info.draw, $
                      curv=info.curv_ind, smoothpressure=info.smoothP, smoothhthe=info.smoothH, smoothCurv=info.smoothCurv, settings=settings
        
        IF poorquality THEN BEGIN
          r = DIALOG_MESSAGE("Poor quality equilibrium", dialog_parent=info.draw)
        ENDIF
        
        plot_rz_equil, *(info.rz_grid)
        oplot_mesh, *(info.rz_grid), *(info.flux_mesh)
      ENDIF ELSE BEGIN
        PRINT, "ERROR: Need to generate a mesh first"
        WIDGET_CONTROL, info.status, set_value="  *** Need to generate mesh first ***"
      ENDELSE
    END
    'print': BEGIN
      IF info.rz_grid_valid THEN BEGIN
        filename = DIALOG_PICKFILE(dialog_parent=event.top, file="bout.grd.ps", $
                                 /write, /overwrite_prompt, path=info.path, get_path=newpath)
        info.path=newpath
        
        IF STRLEN(filename) EQ 0 THEN BEGIN
          WIDGET_CONTROL, info.status, set_value="   *** Cancelled printing ***"
          BREAK
        ENDIF
        SET_PLOT, 'PS'
        DEVICE, file=filename
        plot_mesh, *(info.flux_mesh), xtitle="Major radius [m]", $
          ytitle="Height [m]", title="Generated: "+SYSTIME()

        ; Plot boundary
        IF in_struct(*info.rz_grid, "nlim") THEN BEGIN
          data = *info.rz_grid
          IF data.nlim GT 2 THEN BEGIN
            OPLOT, [REFORM(data.rlim), data.rlim[0]], [REFORM(data.zlim), data.zlim[0]], $
                   thick=2,color=2
          ENDIF
        ENDIF
        
        DEVICE, /close
        SET_PLOT, 'X'
        WIDGET_CONTROL, info.status, set_value="Plotted mesh to file "+filename
      ENDIF ELSE BEGIN
        WIDGET_CONTROL, info.status, set_value="  *** Need to generate mesh first ***"
      ENDELSE
    END
    'curv': BEGIN
      ; Combo box with curvature options
      info.curv_ind = event.index
      widget_control, event.top, set_UVALUE=info
    END
    'strict': BEGIN
      ; Checkbox with boundary strictness
      info.strict_bndry = event.select
      widget_control, event.top, set_UVALUE=info
    END
    'xpt_only': BEGIN
      ; Checkbox for X-point only non-orth option
      info.xpt_only = event.select
      widget_control, event.top, set_UVALUE=info
    END
    'simplebndry': BEGIN
      ; Simplify the boundary
      info.simple_bndry = event.select
      widget_control, event.top, set_UVALUE=info
    END
    'radgrid': BEGIN
      info.single_rad_grid = event.select
      widget_control, event.top, set_UVALUE=info
    END
    'smoothP': BEGIN
      info.smoothP = event.select
      widget_control, event.top, set_UVALUE=info
    END
    'smoothH': BEGIN
      info.smoothH = event.select
      widget_control, event.top, set_UVALUE=info
    END
    'smoothCurv' : BEGIN
      info.smoothCurv = event.select
      widget_control, event.top, set_UVALUE=info
    END
    'calcp' : BEGIN
      info.calcp = event.select
      widget_control, event.top, set_UVALUE=info
    END
    'calcbt' : BEGIN
      info.calcbt = event.select
      widget_control, event.top, set_UVALUE=info
    END
    'calchthe' : BEGIN
      info.calchthe = event.select
      widget_control, event.top, set_UVALUE=info
    END
    'calcjpar' : BEGIN
      info.calcjpar = event.select
      widget_control, event.top, set_UVALUE=info
    END
    'orthogonal_coordinates_output' : BEGIN
      info.orthogonal_coordinates_output = event.select
      widget_control, event.top, set_UVALUE=info
    END
    'fast': BEGIN
      info.fast = event.select
      widget_control, event.top, set_UVALUE=info
    END
    'draw': BEGIN
      IF info.flux_mesh_valid EQ 0 THEN RETURN
      
      pos = CONVERT_COORD(event.x, event.y, /device, /to_data)
      r = pos[0]
      z = pos[1]
      ; (r,z) position where clicked
      
      m = MIN( (REFORM((*info.flux_mesh).rxy - r))^2 + $
               (REFORM((*info.flux_mesh).zxy - r))^2 , ind)
      xi = FIX(ind / TOTAL((*info.flux_mesh).npol + (*info.flux_mesh).n_y_boundary_guards))
      yi = FIX(ind MOD TOTAL((*info.flux_mesh).npol + (*info.flux_mesh).n_y_boundary_guards))
      PRINT, xi, yi
    END
    'detail': BEGIN
      ; Control detailed settings. 
      IF info.rz_grid_valid EQ 0 THEN BEGIN
        WIDGET_CONTROL, info.status, set_value="   *** Cancelled printing ***"
        RETURN
      ENDIF
      
      IF info.flux_mesh_valid THEN BEGIN
        critical = (*info.flux_mesh).critical
      ENDIF ELSE BEGIN
        
        critical = (*(info.rz_grid)).critical
        IF (*(info.rz_grid)).nlim GT 2 THEN BEGIN
          ; Check that the critical points are inside the boundary
          bndryi = DBLARR(2, (*(info.rz_grid)).nlim)
          bndryi[0,*] = INTERPOL(FINDGEN((*(info.rz_grid)).nr), (*(info.rz_grid)).R, (*(info.rz_grid)).rlim)
          bndryi[1,*] = INTERPOL(FINDGEN((*(info.rz_grid)).nz), (*(info.rz_grid)).Z, (*(info.rz_grid)).zlim)
          critical = critical_bndry(critical, bndryi)
        ENDIF
        
        ; Restrict the psi range
        widget_control, info.psi_outer_field, get_value=psi_outer

        critical = restrict_psi_range(critical, psi_outer)
      ENDELSE
      
      n_xpoint = critical.n_xpoint

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Create a new window
      popup = WIDGET_BASE(title="Detailed settings", /COLUMN, $ ; mbar=mbar
                          EVENT_PRO = 'popup_event')
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Number of radial regions

      l = WIDGET_LABEL(popup, value="Number of points in radial direction")
      rad_base = WIDGET_BASE(popup, /ROW, EVENT_PRO='popup_event')
      
      IF info.flux_mesh_valid THEN BEGIN
        ; Use the result from create_grid
        
        nnrad = N_ELEMENTS((*info.flux_mesh).nrad)
        nrad_field = LONARR(nnrad)
        
        nrad_field[0] = CW_FIELD( rad_base,                            $
                                  title  = 'Core:',      $ 
                                  uvalue = 'nrad',                $ 
                                  /long,                          $ 
                                  value = (*info.flux_mesh).nrad[0], $
                                  xsize=8                         $
                                )
      
        IF nnrad GT 1 THEN BEGIN
          IF nnrad GT 2 THEN BEGIN
            ; Regions between separatrices
            
            FOR i=1, nnrad-2 DO BEGIN
              nrad_field[i] = CW_FIELD( rad_base,                            $
                                        title  = 'Inter-separatrix:',      $ 
                                        uvalue = 'nrad',                $ 
                                        /long,                          $ 
                                        value = (*info.flux_mesh).nrad[i], $
                                        xsize=8                         $
                                      )
            ENDFOR
          ENDIF
          
          ; SOL region
          nrad_field[nnrad-1] = CW_FIELD( rad_base,                            $
                                          title  = 'SOL:',      $ 
                                          uvalue = 'nrad',                $ 
                                          /long,                          $ 
                                          value = (*info.flux_mesh).nrad[nnrad-1], $
                                          xsize=8                         $
                                        )
        ENDIF
      ENDIF ELSE BEGIN
        nnrad = 1
        nrad_field = LONARR(nnrad)
        
        widget_control, info.nrad_field, get_value=nrad

        nrad_field[0] = CW_FIELD( rad_base,                       $
                                  title  = 'Total:',              $
                                  uvalue = 'nrad',                $ 
                                  /long,                          $ 
                                  value = nrad,                   $
                                  xsize=8                         $
                                )
      ENDELSE
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Inner psi
      
      l = WIDGET_LABEL(popup, value="Inner psi")
      l = WIDGET_LABEL(popup, value="(Clockwise from innermost x-point)")
      
      in_psi_base = WIDGET_BASE(popup, /ROW, EVENT_PRO='popup_event')
      in_psi_field = LONARR(n_xpoint+1)
      
      IF info.flux_mesh_valid THEN BEGIN
        psi_inner = (*info.flux_mesh).psi_inner
      ENDIF ELSE BEGIN
        widget_control, info.psi_inner_field, get_value=psi_in
        psi_inner = DBLARR(n_xpoint+1) + psi_in
      ENDELSE
      
      in_psi_field[0] = CW_FIELD( in_psi_base,                    $
                                  title  = 'Core: ',              $ 
                                  uvalue = 'in_psi',              $ 
                                  /double,                        $ 
                                  value = psi_inner[0],           $
                                  xsize=8                         $
                                )
      FOR i=1, n_xpoint DO BEGIN
        in_psi_field[i] = CW_FIELD( in_psi_base,                    $
                                    title  = 'PF '+STRTRIM(STRING(i),2)+': ', $ 
                                    uvalue = 'in_psi',              $ 
                                    /double,                        $ 
                                    value = psi_inner[i],           $
                                    xsize=8                         $
                                  )
      ENDFOR
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Outer psi
      
      l = WIDGET_LABEL(popup, value="Outer psi")
      out_psi_base = WIDGET_BASE(popup, /ROW, EVENT_PRO='popup_event')
      
      IF info.flux_mesh_valid THEN BEGIN
        psi_outer = (*info.flux_mesh).psi_outer
      ENDIF ELSE BEGIN
        widget_control, info.psi_outer_field, get_value=psi_out
        psi_outer = DBLARR(n_xpoint) + psi_out
      ENDELSE
      
      out_psi_field = LONARR(N_ELEMENTS(psi_outer))
      
      FOR i=0, N_ELEMENTS(psi_outer)-1 DO BEGIN
        out_psi_field[i] = CW_FIELD( out_psi_base,                    $
                                     title  = 'SOL '+STRTRIM(STRING(i),2)+': ', $ 
                                     uvalue = 'out_psi',              $ 
                                     /double,                         $ 
                                     value = psi_outer[i],            $
                                     xsize=8                          $
                                  )
      ENDFOR
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Poloidal points

      l = WIDGET_LABEL(popup, value="Number of points in poloidal direction")
      l = WIDGET_LABEL(popup, value="(Clockwise from innermost x-point)")
      pol_base = WIDGET_BASE(popup, /ROW, EVENT_PRO='popup_event')
      
      IF info.flux_mesh_valid THEN BEGIN
        nnpol = N_ELEMENTS((*info.flux_mesh).npol)
        
        npol_field = LONARR(nnpol)
      
        IF n_xpoint EQ 0 THEN BEGIN
          npol_field[0] = CW_FIELD( pol_base,                            $
                                    title  = 'Core: ',      $ 
                                    uvalue = 'npol',                $ 
                                    /long,                          $ 
                                    value = (*info.flux_mesh).npol[0], $
                                    xsize=8                         $
                                  )
        ENDIF ELSE BEGIN
          FOR i=0, nnpol-1 DO BEGIN
            IF i MOD 3 EQ 1 THEN title='Core: ' ELSE title  = 'Private Flux: '
            npol_field[i] = CW_FIELD( pol_base,                           $
                                      title  = title,                     $
                                      uvalue = 'npol',                    $
                                      /long,                              $
                                      value = (*info.flux_mesh).npol[i],  $
                                      xsize=8                             $
                                    )
          ENDFOR
        ENDELSE
      ENDIF ELSE BEGIN
        nnpol = 1
        npol_field = LONARR(nnpol)
        
        widget_control, info.npol_field, get_value=npol
        npol_field[0] = CW_FIELD( pol_base,                       $
                                  title  = 'Total: ',             $ 
                                  uvalue = 'npol',                $ 
                                  /long,                          $ 
                                  value = npol, $
                                  xsize=8                         $
                                )
      ENDELSE

      l = WIDGET_LABEL(popup, value='Strength of decay of non-orthogonality away from grid section boundary')
      l = WIDGET_LABEL(popup, value='(Larger exponent pushes the grid back toward orthogonality faster)')
      nonorthogonal_decay_base = WIDGET_BASE(popup, /ROW, EVENT_PRO='popup_event')
      nonorthogonal_weight_decay_power = $
          CW_FIELD( nonorthogonal_decay_base, $
                    title = 'Exponent: ', $
                    uvalue = 'nonorthogonal_weight_decay_power', $
                    /double, $
                    value = 2.7D, $
                    xsize = 8 $
                  )
      
      mesh_button = WIDGET_BUTTON(popup, VALUE='Generate mesh', $
                                  uvalue='mesh', tooltip="Generate a new mesh")
      
      mesh2_button = WIDGET_BUTTON(popup, VALUE='Generate non-orthogonal mesh', $
                                   uvalue='mesh2', tooltip="Generate a new non-orthogonal mesh")
      
      popup_info = {info:info, $ ; Store the main info too
                    nrad_field:nrad_field, $
                    in_psi_field:in_psi_field, $
                    out_psi_field:out_psi_field, $
                    npol_field:npol_field, $
                    nonorthogonal_weight_decay_power:nonorthogonal_weight_decay_power, $
                    top:event.top}

      WIDGET_CONTROL, popup, set_uvalue=popup_info 

      WIDGET_CONTROL, popup, /real
      XMANAGER, 'popup', popup, /just_reg
    END
    'save': BEGIN
      filename = DIALOG_PICKFILE(dialog_parent=event.top, file="hypnotoad.idl", $
                                 /write, /overwrite_prompt, path=info.path, get_path=newpath)
      info.path = newpath
      SAVE, info, file=filename
      WIDGET_CONTROL, info.status, set_value="Saved state to "+filename
    END
    'restore': BEGIN
      filename = DIALOG_PICKFILE(dialog_parent=event.top, file="hypnotoad.idl", /read, path=info.path, get_path=newpath)
      info.path=newpath
      IF STRLEN(filename) EQ 0 THEN BEGIN
        WIDGET_CONTROL, info.status, set_value="   *** Cancelled restore ***"
        RETURN
      ENDIF
      oldinfo = info
      
      CATCH, err
      IF err NE 0 THEN BEGIN
        WIDGET_CONTROL, info.status, set_value=!ERROR_STATE.MSG
      ENDIF ELSE BEGIN

        RESTORE, filename
        
        ; Copy the widget IDs, adding fields if needed for backwards compatability
        str_set, info, "nrad_field", oldinfo.nrad_field, /over
        str_set, info, "npol_field", oldinfo.npol_field, /over
        str_set, info, "draw", oldinfo.draw, /over
        str_set, info, "psi_inner_field", oldinfo.psi_inner_field, /over
        str_set, info, "psi_outer_field", oldinfo.psi_outer_field, /over
        str_set, info, "rad_peak_field", oldinfo.rad_peak_field, /over
        str_set, info, "xpt_dist_field", oldinfo.xpt_dist_field, /over
        str_set, info, "nonorthogonal_weight_decay_power", oldinfo.nonorthogonal_weight_decay_power, /over
        str_set, info, "y_boundary_guards_field", oldinfo.y_boundary_guards_field, /over

        str_set, info, "status", oldinfo.status, /over
        str_set, info, "leftbargeom", oldinfo.leftbargeom, /over

        ; Restore options bar settings
        str_set, info, "curv_select", oldinfo.curv_select, /over
        str_set, info, "curv_ind", oldinfo.curv_ind
        WIDGET_CONTROL, info.curv_select, set_combobox_select=info.curv_ind
        
        str_set, info, "strict_check", oldinfo.strict_check, /over
        str_set, info, "strict_bndry", oldinfo.strict_bndry
        Widget_Control, info.strict_check, Set_Button=info.strict_bndry
        
        str_set, info, "simple_check", oldinfo.simple_check, /over
        str_set, info, "simple_bndry", oldinfo.simple_bndry
        Widget_Control, info.simple_check, Set_Button=info.simple_bndry
        
        str_set, info, "smoothP_check", oldinfo.smoothP_check, /over
        str_set, info, "smoothP", oldinfo.smoothP
        Widget_Control, info.smoothP_check, Set_Button=info.smoothP
        
        str_set, info, "smoothH_check", oldinfo.smoothH_check, /over
        str_set, info, "smoothH", oldinfo.smoothH
        Widget_Control, info.smoothH_check, Set_Button=info.smoothH
        
        str_set, info, "smoothCurv_check", oldinfo.smoothCurv_check, /over
        str_set, info, "smoothCurv", oldinfo.smoothCurv
        Widget_Control, info.smoothCurv_check, Set_Button=info.smoothCurv

        str_set, info, "calcp_check", oldinfo.calcp_check, /over        
        str_set, info, "calcp", oldinfo.calcp
        Widget_Control, info.calcp_check, Set_Button=info.calcp
        
        str_set, info, "calcbt_check", oldinfo.calcbt_check, /over      
        str_set, info, "calcbt", oldinfo.calcbt
        Widget_Control, info.calcbt_check, Set_Button=info.calcbt
        
        str_set, info, "calchthe_check", oldinfo.calchthe_check, /over      
        str_set, info, "calchthe", oldinfo.calchthe
        Widget_Control, info.calchthe_check, Set_Button=info.calchthe
        
        str_set, info, "calcjpar_check", oldinfo.calcjpar_check, /over      
        str_set, info, "calcjpar", oldinfo.calcjpar
        Widget_Control, info.calcjpar_check, Set_Button=info.calcjpar

        str_set, info, "orthogonal_coordinates_output_check", oldinfo.orthogonal_coordinates_output_check, /over      
        str_set, info, "orthogonal_coordinates_output", oldinfo.orthogonal_coordinates_output
        Widget_Control, info.orthogonal_coordinates_output_check, Set_Button=info.orthogonal_coordinates_output

        str_set, info, "radgrid_check", oldinfo.radgrid_check, /over
        str_set, info, "single_rad_grid", oldinfo.single_rad_grid
        Widget_Control, info.radgrid_check, Set_Button=info.single_rad_grid
        
        str_set, info, "fast_check", oldinfo.fast_check, /over
        str_set, info, "fast", oldinfo.fast
        Widget_Control, info.fast_check, Set_Button=info.fast
        
        str_set, info, "path", oldinfo.path, /over
        
        IF info.rz_grid_valid THEN BEGIN
          plot_rz_equil, *info.rz_grid
        ENDIF
        
        IF info.flux_mesh_valid THEN BEGIN
          oplot_mesh, *info.rz_grid, *info.flux_mesh
        ENDIF
        
        widget_control, event.top, set_UVALUE=info
        WIDGET_CONTROL, info.status, set_value="Restored state from "+filename
      ENDELSE
    END
    ELSE: PRINT, "Unknown event", uvalue
  ENDCASE
END

PRO handle_resize, event
  
  IF WHERE(TAG_NAMES(event) EQ "X") EQ -1 THEN RETURN
  
  WIDGET_CONTROL, event.top, get_uvalue=info, /No_Copy
  
  statusgeom = WIDGET_INFO(info.status, /geom) 
  
  WIDGET_CONTROL, info.draw, $
                  Draw_XSize=(event.x - info.leftbargeom.xsize) > statusgeom.xsize, $
                  Draw_YSize=(event.y - statusgeom.ysize) > info.leftbargeom.ysize
  
  IF info.rz_grid_valid THEN BEGIN
    ; Plot the equilibrium
    plot_rz_equil, *info.rz_grid

     IF info.flux_mesh_valid THEN BEGIN
       ; Overplot the mesh
       mesh = *info.flux_mesh
       oplot_mesh, *info.rz_grid, *info.flux_mesh
     ENDIF
  ENDIF

  Widget_Control, event.top, Set_UValue=info, /No_Copy
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Main procedure

PRO hypnotoad
  
  ; Make IDL retain a backing store
  DEVICE, retain=2

  ; Create the main window
  base = WIDGET_BASE(title="Hypnotoad", /ROW, $ ; mbar=mbar
                     EVENT_PRO = 'event_handler', TLB_size_events=1)

  ; Put items in the menu
  ;input_menu = WIDGET_BUTTON(mbar, VALUE='Input', /MENU)
  ;input_bttn1=WIDGET_BUTTON(input_menu, VALUE='Open G-eqdsk (neqdsk)...',$
  ;                          UVALUE='aandg', EVENT_PRO = 'event_handler')
  ;input_bttn2=WIDGET_BUTTON(input_menu, VALUE='From IDAM...',$
  ;                          UVALUE='idam', EVENT_PRO = 'event_handler')
  ;input_bttn2=WIDGET_BUTTON(input_menu, VALUE='Test data...',$
  ;                          UVALUE='test', EVENT_PRO = 'event_handler')

  ; Create a bar down left side for buttons and settings
  bar = WIDGET_BASE(base, /COLUMN, EVENT_PRO = 'event_handler')
  
  ; Create a tabbed interface in the bar
  tab_base = WIDGET_TAB(base)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Mesh tab, for generating new meshes

  tab1 = WIDGET_BASE(tab_base, title="Mesh", /Column, EVENT_PRO = 'event_handler')
  
  read_button = WIDGET_BUTTON(tab1, VALUE='Read G-EQDSK', $
                              uvalue='aandg', tooltip="Read RZ equilibrium from EFIT")

  restore_rz_button = WIDGET_BUTTON(tab1, VALUE='Restore R-Z', $
                                    uvalue='restorerz', tooltip="Restore R-Z equilibrium")

  bndry_button = WIDGET_BUTTON(tab1, VALUE='Read boundary', $
                               uvalue='bndry', tooltip="Read boundary from g-eqdsk file")
  
  nrad_field = CW_FIELD( tab1,                            $
                         title  = 'Radial points:',      $ 
                         uvalue = 'nrad',                $ 
                         /long,                          $ 
                         value = 36,                     $
                         xsize=8                         $
                       )
  npol_field = CW_FIELD( tab1,                            $
                         title  = 'Poloidal points:',    $ 
                         uvalue = 'npol',                $ 
                         /long,                          $ 
                         value = 64,                     $
                         xsize=8                         $
                       )
  
  psi_inner_field = CW_FIELD( tab1,                            $
                              title  = 'Inner psi:',          $ 
                              uvalue = 'inner_psi',           $ 
                              /double,                        $ 
                              value = 0.9D,                    $
                              xsize=8                         $
                            )
  psi_outer_field = CW_FIELD( tab1,                            $
                              title  = 'Outer psi:',          $ 
                              uvalue = 'outer_psi',           $ 
                              /double,                        $ 
                              value = 1.1D,                    $
                              xsize=8                         $
                            )
  
  
  rad_peak_field = CW_FIELD( tab1,                            $
                             title  = 'Sep. spacing:',          $ 
                             uvalue = 'rad_peak',           $ 
                             /double,                        $ 
                             value = 1,                    $
                             xsize=8                         $
                           )
  
  parweight_field = CW_FIELD( tab1,                            $
                             title  = 'Par. vs pol:',          $ 
                             uvalue = 'parweight',           $ 
                             /double,                        $ 
                             value = 0.0D,                    $
                             xsize=8                         $
                           )

  xpt_dist_field = CW_FIELD( tab1,                            $
                             title  = 'Xpt dist x:',          $ 
                             uvalue = 'xpt_mul',           $ 
                             /double,                        $ 
                             value = 1,                    $
                             xsize=8                         $
                           )

  y_boundary_guards_field = CW_FIELD( tab1,                               $
                                      title = '# y-boundary guard cells', $
                                      uvalue = 'y_boundary_guards',       $
                                      /long,                              $
                                      value = 0,                          $
                                      xsize = 3                           $
                                    )

  l = WIDGET_LABEL(tab1, value = '(default 0 for backward compatibility,' + STRING(10B) $
                                 + 'recommended to set to number of' + STRING(10B)      $
                                 + 'y-guards in your simulation, e.g. 2)',              $
                                 /ALIGN_LEFT)
  
  detail_button = WIDGET_BUTTON(tab1, VALUE='Detailed settings', $
                                uvalue='detail', $
                                tooltip="Set quantities in each region")

  save_button = WIDGET_BUTTON(tab1, VALUE='Save state', $
                               uvalue='save', tooltip="Save current Hypnotoad state")
  restore_button = WIDGET_BUTTON(tab1, VALUE='Restore state', $
                               uvalue='restore', tooltip="Restore Hypnotoad state")
  
  ; Options
  checkboxbase = WIDGET_BASE(tab1, /COLUMN, EVENT_PRO = 'event_handler', /NonExclusive)
  strict_check = WIDGET_BUTTON(checkboxbase, VALUE="Strict boundaries", uvalue='strict', $
                               tooltip="Enforce boundaries strictly")
  Widget_Control, strict_check, Set_Button=0

  simple_check = WIDGET_BUTTON(checkboxbase, VALUE="Simplify boundary", uvalue='simplebndry', $
                               tooltip="Simplify the boundary to a square")
  Widget_Control, simple_check, Set_Button=0
  
  radgrid_check = WIDGET_BUTTON(checkboxbase, VALUE="Single radial grid", uvalue='radgrid', $
                             tooltip="Grid radially in one")
  Widget_Control, radgrid_check, Set_Button=1
  
  fast_check = WIDGET_BUTTON(checkboxbase, VALUE="Fast", uvalue='fast', tooltip="Uses faster but less acurate methods")
  Widget_Control, fast_check, set_button=0
  
  xptonly_check = WIDGET_BUTTON(checkboxbase, VALUE="Non-orth X-point only", uvalue='xpt_only', $
                               tooltip="Allows strikepoints to be orthogonal, but keeps X-point non-orthogonal")
  Widget_Control, xptonly_check, Set_Button=0
  
  gen_mesh_text = WIDGET_LABEL(tab1, VALUE='Generate Mesh:', frame=0)

  mesh_button = WIDGET_BUTTON(tab1, VALUE='Orthogonal mesh', $
                              uvalue='mesh', tooltip="Generate a new orthogonal mesh")

  mesh2_button = WIDGET_BUTTON(tab1, VALUE='Non-orthogonal mesh', $
                              uvalue='mesh2', tooltip="Generate a new non-orthogonal mesh")

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Profiles tab

  tab3 = WIDGET_BASE(tab_base, title="Profiles", /COLUMN, EVENT_PRO = 'event_handler')
  

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Output tab
  
  tab2 = WIDGET_BASE(tab_base, title="Output", /COLUMN, EVENT_PRO = 'event_handler')

  w = WIDGET_LABEL(tab2, value="Curvature method")
  curv_select = WIDGET_COMBOBOX(tab2, VALUE=["Toroidal, SVD method", $
                                            "Cylindrical+interpol", $
                                             "Curl(b/B)", $
                                            "Field-aligned coords"], $
                                EVENT_PRO = 'event_handler', $
                                UVALUE="curv")
  curv_index = 1 ; Default index
  WIDGET_CONTROL, curv_select, set_combobox_select=curv_index

  checkboxbase = WIDGET_BASE(tab2, /COLUMN, EVENT_PRO = 'event_handler', /NonExclusive)
  

  smoothP_check = WIDGET_BUTTON(checkboxbase, VALUE="Smooth pressure", uvalue='smoothP', $
                                tooltip="Smooth P profiles")
  
  Widget_Control, smoothP_check, Set_Button=0
  
  smoothH_default = 1
  smoothH_check = WIDGET_BUTTON(checkboxbase, VALUE="Smooth Hthe", uvalue='smoothH', $
                                tooltip="Smooth Hthe")

  Widget_Control, smoothH_check, Set_Button=smoothH_default

  smoothCurv_default = 1
  smoothCurv_check = WIDGET_BUTTON(checkboxbase, VALUE="Smooth Curvature", uvalue='smoothCurv', $
                                tooltip="Smooth Curvature")
  Widget_Control, smoothCurv_check, Set_Button=smoothCurv_default
  
  
  calcp_default = 0
  calcp_check = WIDGET_BUTTON(checkboxbase, VALUE="Recalc P", uvalue='calcp', $
                              tooltip="Recalculate Pressure")
  Widget_Control, calcp_check, Set_Button=calcp_default

  calcbt_default = 0
  calcbt_check = WIDGET_BUTTON(checkboxbase, VALUE="Recalc RBt", uvalue='calcbt', $
                              tooltip="Recalculate R*Bt")
  Widget_Control, calcbt_check, Set_Button=calcbt_default

  calchthe_default = 0
  calchthe_check = WIDGET_BUTTON(checkboxbase, VALUE="Recalc hthe", uvalue='calchthe', $
                              tooltip="Recalculate hthe")
  Widget_Control, calchthe_check, Set_Button=calchthe_default
  
  calcjpar_default = 0
  calcjpar_check = WIDGET_BUTTON(checkboxbase, VALUE="Recalc Jpar", uvalue='calcjpar', $
                              tooltip="Recalculate Jpar")
  Widget_Control, calcjpar_check, Set_Button=calcjpar_default

  orthogonal_coordinates_output_default = 0
  orthogonal_coordinates_output_check = WIDGET_BUTTON(checkboxbase, $
            VALUE="Output for orthogonal coords", uvalue='orthogonal_coordinates_output', $
            tooltip="Output metrics for simulations in orthogonal coordinates using ShiftedMetric (i.e. with zero integrated shear, I=0, when calculating metric terms).")
  Widget_Control, orthogonal_coordinates_output_check, Set_Button=orthogonal_coordinates_output_default

  process_button = WIDGET_BUTTON(tab2, VALUE='Output mesh', $
                                 uvalue='process', tooltip="Process mesh and output to file")

  print_button = WIDGET_BUTTON(tab2, VALUE='Plot to file', $
                               uvalue='print', tooltip="Produce a Postscript plot of the mesh")

  leftbargeom = WIDGET_INFO(bar, /Geometry)

  rightbar = WIDGET_BASE(base, /COLUMN, EVENT_PRO = 'event_handler')
  status_box = WIDGET_TEXT(rightbar)
  ; Create an area for drawing
  draw = WIDGET_DRAW(rightbar, xsize=400, ysize=600, /button_events, $
                     uvalue="draw", EVENT_PRO = 'event_handler')

  widget_control, status_box, set_value="Hypnotoad flux grid generator. Read equilibrium G-EQDSK file to begin"

  ; Get current working directory
  CD, current=path

  ; Create a structure for storing the state
  ; This is shared 

  info = { nrad_field:nrad_field, $ ; nrad input box
           npol_field:npol_field, $ ; npol input box
           rz_grid:(PTRARR(1))[0], $ ; Pointer to R-Z grid data
           rz_grid_valid:0, $  ; Flag to indicate if rz_mesh is valid
           flux_mesh:(PTRARR(1))[0], $ ; Pointer to flux surface aligned mesh
           flux_mesh_valid:0, $
           boundary:(PTRARR(1))[0], $ ; Pointer to boundary array [2,*]
           boundary_valid: 0, $
           settings:(PTRARR(1))[0], $ ; Settings structure
           draw:draw, $ ; Drawing widget
           detail_set:0, $ ; 1 if using detailed settings
           psi_inner_field:psi_inner_field, psi_outer_field:psi_outer_field, $
           rad_peak_field:rad_peak_field, $
           parweight_field:parweight_field, $
           xpt_dist_field:xpt_dist_field, $
           y_boundary_guards_field:y_boundary_guards_field, $
           status:status_box, $
           leftbargeom:leftbargeom, $
           $;;; Options tab 
           curv_select:curv_select, $
           curv_ind:curv_index, $ 
           strict_check:strict_check, $
           strict_bndry:0, $ ; 1 if boundaries should be strict
           simple_check:simple_check, $
           simple_bndry:0, $ ; Use simplified boundary?
           xptonly_check:xptonly_check, $ ; 
           xpt_only:0, $ ; x-point only non-orthogonal
           nonorthogonal_weight_decay_power:2.7D, $ ; how fast to decay towards orthogonal mesh
           radgrid_check:radgrid_check, $
           single_rad_grid:1, $
           smoothP_check:smoothP_check, $
           smoothP:0, $     ; Interpolate to make P smooth
           smoothH_check:smoothH_check, $
           smoothH:smoothH_default, $ 
           smoothCurv_check:smoothCurv_check, $
           smoothCurv:smoothCurv_default, $
           calcp_check:calcp_check, $
           calcp:calcp_default, $
           calcbt_check:calcbt_check, $
           calcbt:calcbt_default, $
           calchthe_check:calchthe_check, $
           calchthe:calchthe_default, $
           calcjpar_check:calcjpar_check, $
           calcjpar:calcjpar_default, $
           orthogonal_coordinates_output_check:orthogonal_coordinates_output_check, $
           orthogonal_coordinates_output:orthogonal_coordinates_output_default, $
           fast_check:fast_check, $
           fast:0, $
           $;;;
           path:path $
         } 

  ; Store this in the base UVALUE
  WIDGET_CONTROL, base, set_uvalue=info 

  ; Draw everything
  WIDGET_CONTROL, base, /real

  XMANAGER, 'hypnotoad', base, /no_block, event_handler='handle_resize'
END
