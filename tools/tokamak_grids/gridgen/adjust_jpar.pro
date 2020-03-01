; Adjusts Jpar0 to get force balance in reduced MHD models
; i.e. the equilibrium terms in the vorticity equation: 
;
; B^2 Grad_par(Jpar0/B) + 2*b0xk dot Grad(P0) = 0
;
; First finds the Pfirsch-Schluter current, then matches
; the current with the input at the inboard midplane
; 
; Usage:
; 
; IDL> adjust_jpar, "some.grid.nc"
; 
; will calculate the new Jpar, then offer to write it to the file
; (overwriting the old value of Jpar)
;
; IDL> g = file_import("some.grid.nc")
; IDL> adjust_jpar, g, jpar=jpar
;
; Given a grid structure, just calculates the new Jpar
;
; B.Dudson, May 2011
;


FUNCTION grad_par, var, mesh
  dtheta = 2.D*!DPI / DOUBLE(TOTAL(mesh.npol))
  RETURN, (mesh.Bpxy / (mesh.Bxy * mesh.hthe)) * ddy(var, mesh)*dtheta / mesh.dy
END

PRO adjust_jpar, grid, smoothp=smoothp, jpar=jpar, noplot=noplot
  
  type = SIZE(grid, /type)
  IF type EQ 7 THEN BEGIN
    ; Input is a string. Read in the data
    data = file_import(grid)
  ENDIF ELSE IF type EQ 8 THEN BEGIN
    ; A structure, hopefully containing the grid data
    data = grid
  ENDIF ELSE BEGIN
    PRINT, "ERROR: Not sure what to do with this type of grid input"
    RETURN
  ENDELSE
  
  ; Find the inboard midplane. Use inboard since this is maximum B
  ; Matching here rather than outboard produces more realistic results
  ; (current doesn't reverse direction at edge)
  mid_ind = -1
  status = gen_surface_hypnotoad(mesh=data) ; Start generator
  REPEAT BEGIN
    yi = gen_surface_hypnotoad(last=last, xi=xi, period=period)
  
    IF period THEN BEGIN
      mr = MIN(data.rxy[xi, yi], mid_ind)
      mr = MAX(data.rxy[xi, yi], out_mid)
      BREAK
    ENDIF
  ENDREP UNTIL last
  IF mid_ind LT 0 THEN BEGIN
    PRINT, "ERROR: No closed flux surfaces?"
    RETURN
  ENDIF
  
  ; Calculate 2*b0xk dot Grad P
  
  kp = 2.D*data.bxcvx*DDX(data.psixy, data.pressure)
  
  ; Calculate B^2 Grad_par(Jpar0)
  
  gj = data.Bxy^2 * Grad_par(data.Jpar0/data.Bxy, data)
  
  ; Generate Jpar0 by integrating kp (Pfirsch-Schluter current)
  ; Grad_par = (Bp / (B*hthe))*d/dy
  
  gparj = -kp * data.hthe / (data.Bxy * data.Bpxy)
  ps = data.Bxy * int_y(gparj, data, /nosmooth, /simple) * data.dy
  
  
  ; In core region add divergence-free parallel current to match input at
  ; inboard midplane. Using inboard as if the outboard is matched then
  ; unphysical overshoots in jpar can result on the inboard side
  
  ; Need to make sure this bootstrap current is always in the same
  ; direction 
  dj = data.jpar0[*,mid_ind] - ps[*,mid_ind]
  m = MAX(ABS(dj), ind)
  s = SIGN(dj[ind])

  w = WHERE(dj * s LT 0.0D, count) ; find where contribution reverses
  IF count GT 0 THEN dj[w] = 0.0D ; just zero in this region

  jpar = ps
  status = gen_surface_hypnotoad(mesh=data) ; Start generator
  REPEAT BEGIN
    yi = gen_surface_hypnotoad(last=last, xi=xi, period=period)
    
    IF NOT period THEN BEGIN
      ; Due to multi-point differencing, dp/dx can be non-zero outside separatrix
      ps[xi,yi] = 0.0D
      jpar[xi,yi] = 0.0D
    ENDIF

    w = WHERE(yi EQ mid_ind, count)
    IF (count NE 0) AND period THEN BEGIN
      ; Crosses midplane
      
      dj_b = dj[xi] / data.Bxy[xi,mid_ind]
      jpar[xi,yi] = jpar[xi,yi] + dj_b * data.Bxy[xi,yi]
    ENDIF
  ENDREP UNTIL last
  
  IF NOT KEYWORD_SET(noplot) THEN BEGIN
    WINDOW, xsize=800, ysize=800
    !P.multi=[0,2,2,0,0]
    SURFACE, data.jpar0, tit="Input Jpar0", chars=2
    SURFACE, jpar, tit="New Jpar0", chars=2
    PLOT, data.jpar0[0,*], tit="jpar at x=0 Solid=input", yr=[MIN([data.jpar0[0,*],jpar[0,*]]), $
                                                               MAX([data.jpar0[0,*],jpar[0,*]])]
    OPLOT, jpar[0,*], psym=1
  
    ;x = data.ixseps1-1
    ;PLOT, data.jpar0[x,*], tit="Jpar at x="+STR(x)+" Solid=input", $
    ;  yr=[MIN([data.jpar0[x,*],jpar[x,*]]), $
    ;      MAX([data.jpar0[x,*],jpar[x,*]])]
    ;OPLOT, jpar[x,*], psym=1
    
    y = out_mid
    PLOT, data.jpar0[*,y], tit="Jpar at y="+STR(y)+" Solid=input", $
      yr=[MIN([data.jpar0[*,y],jpar[*,y]]), $
          MAX([data.jpar0[*,y],jpar[*,y]])]
    OPLOT, jpar[*,y], psym=1
    
    !P.multi=0
  ENDIF
  
  IF type EQ 7 THEN BEGIN
    ; Ask if user wants to write this new Jpar to file
    
    IF get_yesno("Write new Jpar0 to file?") THEN BEGIN
      
      f = file_open(grid, /write) ; Read/write mode
      
      status = file_write(f, "Jpar0", jpar)
      
      IF status THEN BEGIN
        PRINT, "ERROR writing Jpar0 to file '"+ grid+"'"
      ENDIF
      
      file_close, f
    ENDIF
  ENDIF
END
