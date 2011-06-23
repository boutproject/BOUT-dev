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

PRO plot_rz_equil, data
  nlev = 100
  minf = MIN(data.psi)
  maxf = MAX(data.psi)
  levels = findgen(nlev)*(maxf-minf)/FLOAT(nlev-1) + minf
  
  safe_colors, /first
  CONTOUR, data.psi, data.r, data.z, levels=levels, /iso, color=1
  IF data.nlim GT 2 THEN OPLOT, [data.rlim, data.rlim[0]], $
                                [data.zlim, data.zlim[0]], $
                                color = 2, thick=2
  critical = data.critical
  IF data.nlim GT 2 THEN BEGIN
    OPLOT, [REFORM(data.rlim), data.rlim[0]], [REFORM(data.zlim), data.zlim[0]], $
           thick=2,color=2
    
    ; Check that the critical points are inside the boundary
    bndryi = FLTARR(2, data.nlim)
    bndryi[0,*] = INTERPOL(FINDGEN(data.nr), data.R, data.rlim)
    bndryi[1,*] = INTERPOL(FINDGEN(data.nz), data.Z, data.zlim)
    
    critical = critical_bndry(critical, bndryi)
  ENDIF
  ; Overplot the separatrices, O-points
  oplot_critical, data.psi, data.r, data.z, critical
END


PRO oplot_mesh, rz_mesh, flux_mesh
  ; Plot X-points and separatrices
  FOR i=0, flux_mesh.critical.n_xpoint-1 DO BEGIN
    ; plot the separatrix contour
    CONTOUR, rz_mesh.psi, rz_mesh.R, rz_mesh.Z, levels=[flux_mesh.critical.xpt_f[i]], c_colors=2, /overplot
    oplot, [INTERPOLATE(rz_mesh.R, flux_mesh.critical.xpt_ri[i])], [INTERPOLATE(rz_mesh.Z, flux_mesh.critical.xpt_zi[i])], psym=7, color=2
  ENDFOR

  ; Plot O-points
  FOR i=0, flux_mesh.critical.n_opoint-1 DO BEGIN
    oplot, [INTERPOLATE(rz_mesh.R, flux_mesh.critical.opt_ri[i])], [INTERPOLATE(rz_mesh.Z, flux_mesh.critical.opt_zi[i])], psym=7, color=3
  ENDFOR
  
  ypos = 0
  FOR i=0, N_ELEMENTS(flux_mesh.npol)-1 DO BEGIN
    plot_region, flux_mesh.rxy, flux_mesh.zxy, ypos, ypos+flux_mesh.npol[i]-1, color=i+1
    ypos = ypos + flux_mesh.npol[i]
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
  
  CASE uvalue OF
    'mesh': BEGIN
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
      psi_inner = FLTARR(ninpsi)
      FOR i=0, ninpsi-1 DO BEGIN
        widget_control, info.in_psi_field[i], get_value=inp
        psi_inner[i] = inp
      ENDFOR

      nnpol = N_ELEMENTS(info.npol_field)
      npol = LONARR(nnpol)
      FOR i=0, nnpol-1 DO BEGIN
        widget_control, info.npol_field[i], get_value=np
        npol[i] = np
      ENDFOR
      
      widget_control, base_info.psi_outer_field, get_value=psi_outer

      widget_control, base_info.rad_peak_field, get_value=rad_peak
      
      settings = {nrad:nrad, npol:npol, psi_inner:psi_inner, psi_outer:psi_outer}
      
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
      
      ; Create the mesh
      mesh = create_grid((*(base_info.rz_grid)).psi, (*(base_info.rz_grid)).r, (*(base_info.rz_grid)).z, $
                         settings, $
                         boundary=boundary, strict=base_info.strict_bndry, rad_peaking=rad_peak, $
                         single_rad_grid=base_info.single_rad_grid, $
                         critical=(*(base_info.rz_grid)).critical)
      
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
  
  ; Retrieve a copy of information stored in tlb
  widget_control, event.top, get_uvalue=info
  
  IF N_ELEMENTS(uvalue) EQ 0 THEN RETURN ; Undefined
  
  CASE uvalue OF
    'aandg': BEGIN
      PRINT, "Open G-eqdsk (neqdsk) file"
      filename = DIALOG_PICKFILE(dialog_parent=event.top, file="neqdsk", /read)
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
    'bndry': BEGIN
      IF info.rz_grid_valid NE 1 THEN BEGIN
        PRINT, "ERROR: Need to read an equilibrium first"
        WIDGET_CONTROL, info.status, set_value="   *** Need to read equilibrium first ***"
        BREAK
      ENDIF
      PRINT, "Read boundary from G-eqdsk (neqdsk) file"
      filename = DIALOG_PICKFILE(dialog_parent=event.top, file="neqdsk", /read)
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

        settings = {nrad:nrad, npol:npol, psi_inner:psi_inner, psi_outer:psi_outer}
      ENDELSE
      

      ; Check if a simplified boundary should be used
      IF info.simple_bndry THEN BEGIN
        ; Simplify the boundary to a square box
        boundary = TRANSPOSE([ [MIN(boundary[0,*]), MAX(boundary[0,*]), $
                                MAX(boundary[0,*]), MIN(boundary[0,*])], $
                               [MIN(boundary[1,*]), MIN(boundary[1,*]), $
                                MAX(boundary[1,*]), MAX(boundary[1,*])] ])
      ENDIF
        
      WIDGET_CONTROL, info.status, set_value="Generating mesh ..."
      
      mesh = create_grid((*(info.rz_grid)).psi, (*(info.rz_grid)).r, (*(info.rz_grid)).z, settings, $
                         boundary=boundary, strict=info.strict_bndry, rad_peaking=rad_peak, $
                         /nrad_flexible, $
                         single_rad_grid=info.single_rad_grid, $
                         critical=(*(info.rz_grid)).critical)
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

        settings = {nrad:nrad, npol:npol, psi_inner:psi_inner, psi_outer:psi_outer}
      ENDELSE
      

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
                         boundary=boundary, strict=info.strict_bndry, rad_peaking=rad_peak, $
                         /nrad_flexible, $
                         single_rad_grid=info.single_rad_grid, $
                         critical=(*(info.rz_grid)).critical)
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
                                 /write, /overwrite_prompt)
      
      IF STRLEN(filename) EQ 0 THEN BEGIN
        WIDGET_CONTROL, info.status, set_value="   *** Cancelled process mesh ***"
        BREAK
      ENDIF

      IF info.rz_grid_valid AND info.flux_mesh_valid THEN BEGIN
        
        process_grid, *(info.rz_grid), *(info.flux_mesh), $
                      output=filename, poorquality=poorquality, /gui, parent=info.draw, $
                      curv=info.curv_ind, smoothpressure=info.smoothP
        
        IF poorquality THEN BEGIN
          r = DIALOG_MESSAGE("Poor quality equilibrium", dialog_parent=info.draw)
        ENDIF
      ENDIF ELSE BEGIN
        PRINT, "ERROR: Need to generate a mesh first"
        WIDGET_CONTROL, info.status, set_value="  *** Need to generate mesh first ***"
      ENDELSE
    END
    'print': BEGIN
      IF info.rz_grid_valid THEN BEGIN
        filename = DIALOG_PICKFILE(dialog_parent=event.top, file="bout.grd.ps", $
                                 /write, /overwrite_prompt)
        
        IF STRLEN(filename) EQ 0 THEN BEGIN
          WIDGET_CONTROL, info.status, set_value="   *** Cancelled printing ***"
          BREAK
        ENDIF
        SET_PLOT, 'PS'
        DEVICE, file=filename
        plot_mesh, *(info.flux_mesh), xtitle="Major radius [m]", $
          ytitle="Height [m]", title="Generated: "+SYSTIME()
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
    'draw': BEGIN
      IF info.flux_mesh_valid EQ 0 THEN RETURN
      
      pos = CONVERT_COORD(event.x, event.y, /device, /to_data)
      r = pos[0]
      z = pos[1]
      ; (r,z) position where clicked
      
      m = MIN( (REFORM((*info.flux_mesh).rxy - r))^2 + $
               (REFORM((*info.flux_mesh).zxy - r))^2 , ind)
      xi = FIX(ind / TOTAL((*info.flux_mesh).npol))
      yi = FIX(ind MOD TOTAL((*info.flux_mesh).npol))
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
          bndryi = FLTARR(2, (*(info.rz_grid)).nlim)
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
        psi_inner = FLTARR(n_xpoint+1) + psi_in
      ENDELSE
      
      in_psi_field[0] = CW_FIELD( in_psi_base,                    $
                                  title  = 'Core: ',              $ 
                                  uvalue = 'in_psi',              $ 
                                  /float,                          $ 
                                  value = psi_inner[0],           $
                                  xsize=8                         $
                                )
      FOR i=1, n_xpoint DO BEGIN
        in_psi_field[i] = CW_FIELD( in_psi_base,                    $
                                    title  = 'PF '+STRTRIM(STRING(i),2)+': ', $ 
                                    uvalue = 'in_psi',              $ 
                                    /float,                          $ 
                                    value = psi_inner[i],           $
                                    xsize=8                         $
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
            npol_field[i] = CW_FIELD( pol_base,                          $
                                      title  = title,                    $ 
                                      uvalue = 'npol',                   $ 
                                      /long,                             $ 
                                      value = (*info.flux_mesh).npol[i], $
                                      xsize=8                            $
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
      
      mesh_button = WIDGET_BUTTON(popup, VALUE='Generate mesh', $
                                  uvalue='mesh', tooltip="Generate a new mesh")
      
      popup_info = {info:info, $ ; Store the main info too
                    nrad_field:nrad_field, $
                    in_psi_field:in_psi_field, $
                    npol_field:npol_field, $
                    top:event.top}

      WIDGET_CONTROL, popup, set_uvalue=popup_info 

      WIDGET_CONTROL, popup, /real
      XMANAGER, 'popup', popup, /just_reg
    END
    'save': BEGIN
      filename = DIALOG_PICKFILE(dialog_parent=event.top, file="hypnotoad.idl", $
                                 /write, /overwrite_prompt)
      SAVE, info, file=filename
      WIDGET_CONTROL, info.status, set_value="Saved state to "+filename
    END
    'restore': BEGIN
      filename = DIALOG_PICKFILE(dialog_parent=event.top, file="hypnotoad.idl", /read)
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
        
        ; Copy the widget IDs
        info.nrad_field = oldinfo.nrad_field
        info.npol_field = oldinfo.npol_field
        info.draw = oldinfo.draw
        info.psi_inner_field = oldinfo.psi_inner_field
        info.rad_peak_field = oldinfo.rad_peak_field
        info.status = oldinfo.status
        info.leftbargeom = oldinfo.leftbargeom
        info.strict_bndry = oldinfo.strict_bndry
        info.simple_bndry = oldinfo.simple_bndry
        info.smoothP = oldinfo.smoothP
        info.single_rad_grid = oldinfo.single_rad_grid
        
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
  
  read_button = WIDGET_BUTTON(bar, VALUE='Read G-EQDSK', $
                              uvalue='aandg', tooltip="Read RZ equilibrium from EFIT")

  bndry_button = WIDGET_BUTTON(bar, VALUE='Read boundary', $
                               uvalue='bndry', tooltip="Read boundary from g-eqdsk file")
  
  nrad_field = CW_FIELD( bar,                            $
                         title  = 'Radial points:',      $ 
                         uvalue = 'nrad',                $ 
                         /long,                          $ 
                         value = 36,                     $
                         xsize=8                         $
                       )
  npol_field = CW_FIELD( bar,                            $
                         title  = 'Poloidal points:',    $ 
                         uvalue = 'npol',                $ 
                         /long,                          $ 
                         value = 64,                     $
                         xsize=8                         $
                       )
  
  psi_inner_field = CW_FIELD( bar,                            $
                              title  = 'Inner psi:',          $ 
                              uvalue = 'inner_psi',           $ 
                              /floating,                      $ 
                              value = 0.9,                    $
                              xsize=8                         $
                            )
  psi_outer_field = CW_FIELD( bar,                            $
                              title  = 'Outer psi:',          $ 
                              uvalue = 'outer_psi',           $ 
                              /floating,                      $ 
                              value = 1.1,                    $
                              xsize=8                         $
                            )
  
  
  rad_peak_field = CW_FIELD( bar,                            $
                             title  = 'Sep. packing:',          $ 
                             uvalue = 'rad_peak',           $ 
                             /floating,                      $ 
                             value = 1,                    $
                             xsize=8                         $
                           )
  
  w = WIDGET_LABEL(bar, value="Curvature method")
  curv_select = WIDGET_COMBOBOX(bar, VALUE=["Toroidal, SVD method", $
                                            "Cylindrical+interpol", $
                                            "Field-aligned coords"], $
                                EVENT_PRO = 'event_handler', $
                                UVALUE="curv")
  curv_index = 1 ; Default index
  WIDGET_CONTROL, curv_select, set_combobox_select=curv_index

  checkboxbase = WIDGET_BASE(bar, /COLUMN, EVENT_PRO = 'event_handler', /NonExclusive)
  strict_check = WIDGET_BUTTON(checkboxbase, VALUE="Strict boundaries", uvalue='strict', $
                               tooltip="Enforce boundaries strictly")
  Widget_Control, strict_check, Set_Button=1

  simple_check = WIDGET_BUTTON(checkboxbase, VALUE="Simplify boundary", uvalue='simplebndry', $
                               tooltip="Simplify the boundary to a square")
  Widget_Control, simple_check, Set_Button=0
  
  radgrid_check = WIDGET_BUTTON(checkboxbase, VALUE="Single radial grid", uvalue='radgrid', $
                             tooltip="Grid radially in one")
  Widget_Control, radgrid_check, Set_Button=1

  smoothP_check = WIDGET_BUTTON(checkboxbase, VALUE="Smooth pressure", uvalue='smoothP', $
                                tooltip="Interpolate P to smooth derivative")
  
  Widget_Control, smoothP_check, Set_Button=1
  
  mesh_button = WIDGET_BUTTON(bar, VALUE='Generate mesh', $
                              uvalue='mesh', tooltip="Generate a new mesh")

  mesh2_button = WIDGET_BUTTON(bar, VALUE='Nonorthogonal mesh', $
                              uvalue='mesh2', tooltip="Generate a new nonorthogonal mesh")
  
  detail_button = WIDGET_BUTTON(bar, VALUE='Detailed settings', $
                                uvalue='detail', $
                                tooltip="Set quantities in each region")

  process_button = WIDGET_BUTTON(bar, VALUE='Output mesh', $
                                 uvalue='process', tooltip="Process mesh and output to file")

  print_button = WIDGET_BUTTON(bar, VALUE='Plot to file', $
                               uvalue='print', tooltip="Produce a Postscript plot of the mesh")

  save_button = WIDGET_BUTTON(bar, VALUE='Save state', $
                               uvalue='save', tooltip="Save current Hypnotoad state")
  restore_button = WIDGET_BUTTON(bar, VALUE='Restore state', $
                               uvalue='restore', tooltip="Restore Hypnotoad state")

  leftbargeom = WIDGET_INFO(bar, /Geometry)

  rightbar = WIDGET_BASE(base, /COLUMN, EVENT_PRO = 'event_handler')
  status_box = WIDGET_TEXT(rightbar)
  ; Create an area for drawing
  draw = WIDGET_DRAW(rightbar, xsize=400, ysize=600, /button_events, $
                     uvalue="draw", EVENT_PRO = 'event_handler')

  widget_control, status_box, set_value="Hypnotoad flux grid generator. Read equilibrium G-EQDSK file to begin"

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
           strict_bndry:1, $ ; 1 if boundaries should be strict
           simple_bndry:0, $ ; Use simplified boundary?
           smoothP:1, $     ; Interpolate to make P smooth
           single_rad_grid:1, $
           psi_inner_field:psi_inner_field, psi_outer_field:psi_outer_field, $
           rad_peak_field:rad_peak_field, $
           status:status_box, $
           leftbargeom:leftbargeom, $
           curv_ind:curv_index $
         } 

  ; Store this in the base UVALUE
  WIDGET_CONTROL, base, set_uvalue=info 

  ; Draw everything
  WIDGET_CONTROL, base, /real

  XMANAGER, 'hypnotoad', base, /no_block, event_handler='handle_resize'
END
