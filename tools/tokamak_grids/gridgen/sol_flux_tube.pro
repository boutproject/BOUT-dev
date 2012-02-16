;
; Generate a flux tube input file in the SOL of a tokamak equilibrium
;

PRO sol_flux_tube, gfile, psinorm, output=output, ny=ny
  
  IF NOT KEYWORD_SET(output) THEN output="fluxtube"+STR(psinorm)+".grd.nc"
  
  IF NOT KEYWORD_SET(ny) THEN ny = 128

  IF psinorm LE 1.0 THEN BEGIN
    PRINT, "Error: Normalised psi must be greater than 1"
    RETURN
  ENDIF
  
  ; Read the EFIT G-EQDSK file
  g = read_neqdsk(gfile)
  
  ; Check the return
  IF SIZE(g, /TYPE) NE 8 THEN BEGIN
     PRINT, "ERROR: Couldn't read G-EQDSK file '"+gfile+"'"
     RETURN
  ENDIF
  rzgrid = {nr:g.nx, nz:g.ny, $ ; Number of grid points
            r:REFORM(g.r[*,0]), z:REFORM(g.z[0,*]), $ ; R and Z as 1D arrays
            simagx:g.simagx, sibdry:g.sibdry, $       ; Range of psi
            psi:g.psi, $       ; Poloidal flux in Weber/rad on grid points
            pres:g.pres, $     ; Plasma pressure in nt/m^2 on uniform flux grid
            qpsi:g.qpsi, $     ; q values on uniform flux grid
            nlim:g.nlim, rlim:g.xlim, zlim:g.ylim} ; Wall boundary
  
  ; Plot psi
  nlev = 100
  minf = MIN(rzgrid.psi)
  maxf = MAX(rzgrid.psi)
  levels = findgen(nlev)*(maxf-minf)/FLOAT(nlev-1) + minf
  
  safe_colors, /first
  
  CONTOUR, rzgrid.psi, rzgrid.R, rzgrid.Z, $
           levels=levels, color=1, /iso, xstyl=1, ysty=1
  
  ; Find O- and X- points
  critical = analyse_equil(rzgrid.psi, rzgrid.R, rzgrid.Z)
  
  ; Overplot the separatrices, O-points
  oplot_critical, rzgrid.psi, rzgrid.R, rzgrid.Z, critical
  
  ; Use primary O-point and X-point to convert normalised psi into psi
  
  psi_axis = critical.opt_f[critical.primary_opt]
  psi_sep = critical.xpt_f[critical.inner_sep]
  
  ; See how this compares to grid values
  PRINT, "Psi Axis : ", psi_axis, g.simagx
  PRINT, "Psi Bndry: ", psi_sep, g.sibdry
  
  psi_axis = g.simagx
  psi_sep = g.sibdry

  psi = psi_axis + psinorm*(psi_sep - psi_axis)
  
  ;;;;;;;;;;;;;;; Calculate DCT ;;;;;;;;;;;;;;
  
  PRINT, "Calculating DCT..."
  dctpsi = DCT2D(rzgrid.psi)
  PRINT, "Finished DCT"
  
  ; Create a structure containing interpolation settings and data
  interp_data = {nx:g.nx, ny:g.ny, $
                 method:0, $
                 f: rzgrid.psi, $       ; Always include function
                 dct: dctpsi} ; Pass the DCT coefficients
  
  ;;;;;;;;;;;;;;; Get a contour line
  
  contour_lines, rzgrid.psi, findgen(rzgrid.nr), findgen(rzgrid.nz), $
                 levels=[psi], $
                 path_info=info, path_xy=xy
  
  IF N_ELEMENTS(info) GT 1 THEN BEGIN
     ; Find the surface closest to the outboard midplane
      
     ind = closest_line(info, xy, critical.opt_ri[critical.primary_opt], critical.opt_zi[critical.primary_opt])
     info = info[ind]
  ENDIF ELSE info = info[0]
  
  ; Extract r and z indices
  ri = REFORM(xy[0,info.offset:(info.offset+info.n-1)])
  zi = REFORM(xy[1,info.offset:(info.offset+info.n-1)])
  np = N_ELEMENTS(ri)

  ; Refine flux surface location using DCT
  FOR i=0, np-1 DO BEGIN
    follow_gradient, interp_data, rzgrid.R, rzgrid.Z, ri[i], zi[i], psi, $
                     rinew, zinew, status=status
    IF status THEN BEGIN
      PRINT, "ERROR refining surface location"
      RETURN
    ENDIF
    
    ri[i] = rinew
    zi[i] = zinew
  ENDFOR
  
  rpos = INTERPOLATE(rzgrid.R, ri)
  zpos = INTERPOLATE(rzgrid.Z, zi)

  ; Smooth positions
  rpos = SMOOTH(rpos, 4)
  zpos = SMOOTH(zpos, 4)
  
  IF g.nlim GT 2 THEN BEGIN
    ; Find intersections with boundary
    oplot, g.xlim, g.ylim, thick=2
    
    
    cpos = line_crossings(rpos, zpos, 0, $
                          g.xlim, g.ylim, 1, $
                          ncross=ncross, inds1=inds)
    
    print, "Number of crossings: ", ncross
    
    IF ncross NE 2 THEN BEGIN
      PRINT, "HELP! Don't know what to do..."
      STOP
    ENDIF
    
    rpos = rpos[inds[0]:inds[1]]
    zpos = zpos[inds[0]:inds[1]]
    
    ri = ri[inds[0]:inds[1]]
    zi = zi[inds[0]:inds[1]]
    
  ENDIF ELSE BEGIN
    PRINT, "WARNING: No boundary found"
  ENDELSE
    
  oplot, rpos, zpos, color=4, thick=2
  
  ;;;;;;;;; Construct constant arc lengths
  
  drdi = DERIV(rpos)
  dzdi = DERIV(zpos)
  dldi = SQRT(drdi^2 + dzdi^2) ; Length differential
  
  L = int_func(dldi, /simple) ; Integrate function
  
  lpos = max(L) * FINDGEN(ny)/FLOAT(ny-1)
  
  ; Find index
  inds = INTERPOL(findgen(N_ELEMENTS(L)), L, lpos)
  
  rpos = INTERPOLATE(rpos, inds)
  zpos = INTERPOLATE(zpos, inds)
  
  STOP
END
