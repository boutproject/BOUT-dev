;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Generator of continuous surfaces
;
; First call: 
;   status = gen_surface_hypnotoad(mesh=mesh)     - Initialisation
; Subsequent calls
;   yi = gen_surface_hypnotoad(period=period, last=last, xi=xi)
;
; period - Set to 1 if the surface is periodic, 0 otherwise
; last   - Set to 1 if this is the last surface
; xi     - The X index of this surface
;

FUNCTION range, first, last
  IF first LT last THEN BEGIN
    RETURN, first + INDGEN(last - first + 1)
  ENDIF ELSE BEGIN
    RETURN, last + REVERSE(INDGEN(first - last + 1))
  ENDELSE
END

FUNCTION gen_surface_hypnotoad, mesh=mesh, period=period, last=last, xi=xi
  COMMON gen_surf_com, m, ys, xind, nd, domain, visited
  IF KEYWORD_SET(mesh) THEN BEGIN
    ; Starting
    m = mesh
    xind = 0 ; Radial surface
    nd = N_ELEMENTS(mesh.npol) ; Number of domains
    domain = 0 ; The domain to start in
    
    ; Running total of npol to get starting y index
    ys = LONARR(nd)
    FOR i=1, nd-1 DO ys[i] = ys[i-1] + mesh.npol[i-1] + mesh.n_y_boundary_guards[i-1]
    
    ; visited marks which domains have been used
    visited = INTARR(nd)
    
    RETURN, 0
  ENDIF

  IF xind GE TOTAL(m.nrad) THEN BEGIN
    last = 1
    RETURN, 1 ; Error
  ENDIF
  
  ; Get the next surface
  ny = 0
  period = 0 ; Mark as non-periodic
  last = 0 ; mark as not the last
  xi = xind
  REPEAT BEGIN
    IF visited[domain] EQ 1 THEN BEGIN
      ; Already visited this domain
      period = 1 ; Means this domain is periodic
      BREAK
    ENDIF
    
    ; Get the range of indices for this domain
    yi = [range(ys[domain], ys[domain] + m.npol[domain] + m.n_y_boundary_guards[domain] - 1)]
    IF ny EQ 0 THEN yinds = yi ELSE yinds = [yinds, yi]
    ny = ny + m.npol[domain] + m.n_y_boundary_guards[domain]
    
    visited[domain] = 1 ; Mark domain as visited
    
    ; Find next domain
    IF xind LT m.yup_xsplit[domain] THEN BEGIN
      domain = m.yup_xin[domain]
    ENDIF ELSE BEGIN
      domain = m.yup_xout[domain]
    ENDELSE
  ENDREP UNTIL domain LT 0 ; Keep going until hit a boundary

  ; Find a domain which hasn't been visited
  w = WHERE(visited EQ 0, count)
  IF count NE 0 THEN BEGIN
    ; See if there are any regions with boundaries on lower side
    domain = -1
    FOR i=0, count-1 DO BEGIN
      IF xind LT m.ydown_xsplit[w[i]] THEN BEGIN
        d = m.ydown_xin[w[i]]
      ENDIF ELSE BEGIN
        d = m.ydown_xout[w[i]]
      END
      IF d LT 0 THEN BEGIN
        domain = w[i]
        BREAK
      ENDIF
    ENDFOR
    IF domain LT 0 THEN domain = w[0] ; Set the domain to the first one
    
  ENDIF ELSE BEGIN
    ; No domains left - increase x index (if possible)

    xind = xind + 1
    visited = INTARR(nd) ; Set all to zeros again
    domain = 0 ; Start again with the first domain
    IF xind EQ TOTAL(m.nrad) THEN last = 1 ; No more left
  ENDELSE
  
  IF ny EQ 0 THEN RETURN, 2 ; This shouldn't happen
  
  RETURN, yinds
END
