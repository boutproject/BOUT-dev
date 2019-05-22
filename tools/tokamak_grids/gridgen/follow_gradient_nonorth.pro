;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Follow the gradient from a given point to a target f
; 

; Calculate dR/df and dZ/df for use by LSODE
; Input: pos[0] = R, pos[1] = Z
; Output [0] = dR/df = -Bz/B^2 , [1] = dZ/df = Br/B^2
FUNCTION radial_differential_nonorth, fcur, pos
  COMMON rd_com_no, idata, lastgoodf, lastgoodpos, R, Z, ood, boundary, ri0, zi0, tol, vec_comm_up, weightc_up, vec_comm_down, weightc_down, bndry_periodic
  
  local_gradient, idata, pos[0], pos[1], status=status, dfdr=dfdr, dfdz=dfdz
  
  ood = 0

  IF status EQ 0 THEN BEGIN
    ; No error in localgradient.
    lastgoodf = fcur 
    lastgoodpos = pos

    ; If status NE 0 then an error occurred.
    ; Allow dfdz to cause an error so escape LSODE
  ENDIF ELSE ood = 1 ; Out Of Domain

  IF N_ELEMENTS(boundary) GT 1 THEN BEGIN
    ; Got a boundary - check for crossings
    
    cpos = line_crossings([ri0, pos[0]], [zi0, pos[1]], 0, $
                          boundary[0,*], boundary[1,*], bndry_periodic, ncross=ncross, inds2=inds2)
    IF (ncross MOD 2) EQ 1 THEN BEGIN ; Odd number of boundary crossings
      ; Check how far away the crossing is
      IF SQRT( (pos[0] - cpos[0,0])^2 + (pos[1] - cpos[1,0])^2 ) GT tol THEN BEGIN
        ; Need to trigger an error
        ;PRINT, "HIT BOUNDARY:", SQRT( (pos[0] - cpos[0,0])^2 + (pos[1] - cpos[1,0])^2 )
        status = a.bcd ; gibberish
      ENDIF
    ENDIF
  ENDIF
  
  dRdi = INTERPOLATE(DERIV(R), pos[0], /DOUBLE)
  dZdi = INTERPOLATE(DERIV(Z), pos[1], /DOUBLE)
  
  ; Check mismatch between fcur and f ? 
  Br = dfdz/dZdi
  Bz = -dfdr/dRdi
  B2 = Br^2 + Bz^2

  IF KEYWORD_SET(vec_comm_up) AND KEYWORD_SET(vec_comm_down) THEN BEGIN
  ;; V*mod(R)/mod(V) w + (1-w) vec(r)
     newvec = vec_comm_up * SQRT((-Bz/B2/dRdi)^2 + (Br/B2/dZdi)^2) * weightc_up $
              + vec_comm_down * SQRT((-Bz/B2/dRdi)^2 + (Br/B2/dZdi)^2) * weightc_down $
              + (1-weightc_up-weightc_down)*[-Bz/B2/dRdi, Br/B2/dZdi]
     RETURN, newvec
     ;;RETURN, [-Bz/B2/dRdi, Br/B2/dZdi]
     
  ENDIF ELSE IF KEYWORD_SET(vec_comm_up) THEN BEGIN
     newvec = vec_comm_up * SQRT((-Bz/B2/dRdi)^2 + (Br/B2/dZdi)^2) * weightc_up $
              + (1-weightc_up)*[-Bz/B2/dRdi, Br/B2/dZdi]
     RETURN, newvec
  ENDIF ELSE IF KEYWORD_SET(vec_comm_down) THEN BEGIN
     newvec = vec_comm_down * SQRT((-Bz/B2/dRdi)^2 + (Br/B2/dZdi)^2) * weightc_down $
              + (1-weightc_down)*[-Bz/B2/dRdi, Br/B2/dZdi]
     RETURN, newvec
  ENDIF ELSE BEGIN
  ;; stop
     RETURN, [-Bz/B2/dRdi, Br/B2/dZdi]
  ENDELSE
END

;
; interp_data (in)  Data for interpolation of 2D Psi
; R, Z      (in)    1D arrays of position vs. index
; ri0,zi0   (in)    Starting indices
; ftarget   (in)    The f to aim for
; ri,zi     (out)   Final position
; status    (out)   Non-zero if hits a boundary. 1 if out of range at
;                   the start, 2 if before reaching ftarget
; boundary  (in)    Optional 2D array [2,n] of points on the boundary
; fbndry    (out)   If hits boundary, gives final f value
; ibndry    (out)   If hits boundary, index where hit
;
PRO follow_gradient_nonorth, interp_data, R, Z, ri0, zi0, ftarget, ri, zi, status=status, $
                     boundary=boundary, fbndry=fbndry, ibndry=ibndry, $
                     vec_up=vec_up, weight_up=weight_up, $
                     vec_down=vec_down, weight_down=weight_down, $
                     bndry_noperiodic=bndry_noperiodic
  COMMON rd_com_no, idata, lastgoodf, lastgoodpos, Rpos, Zpos, ood, bndry, ri0c, zi0c, tol, vec_comm_up, weightc_up, vec_comm_down, weightc_down, bndry_periodic
  
  tol = 0.1D

  Rpos = R
  Zpos = Z
  
  ibndry = -1
  
  idata = interp_data

  IF KEYWORD_SET(boundary) THEN BEGIN
    bndry = boundary
    ri0c = ri0
    zi0c = zi0
  ENDIF ELSE bndry = 0
  ood = 0

  bndry_periodic = 1
  IF KEYWORD_SET(bndry_noperiodic) THEN bndry_periodic = 0

  IF NOT KEYWORD_SET(weight_up) THEN BEGIN
     weight_up = 0.D
     weightc_up = weight_up*1.0D
  ENDIF
  IF NOT KEYWORD_SET(weight_down) THEN BEGIN
     weight_down = 0.D
     weightc_down = weight_down*1.0D
  ENDIF

  IF KEYWORD_SET(vec_up) THEN BEGIN
     vec_comm_up = vec_up
     weightc_up = weight_up*1.0D
  ENDIF
  IF KEYWORD_SET(vec_down) THEN BEGIN
     vec_comm_down = vec_down
     weightc_down = weight_down*1.0D
  ENDIF

  IF SIZE(ftarget, /TYPE) EQ 0 THEN PRINT, ftarget

  ; Get starting f
  local_gradient, interp_data, ri0, zi0, status=status, f=f0
  IF status EQ 1 THEN BEGIN
    ri = ri0
    zi = zi0
    status = 1
    RETURN
  ENDIF

  fmax = ftarget ; Target (with maybe boundary in the way)

  CATCH, theError
  IF theError EQ 0 THEN BEGIN
    ; Call LSODE to follow gradient
    rzold = [ri0, zi0]
    rcount = 0
    REPEAT BEGIN
      rznew = LSODE(rzold,f0,ftarget - f0,'radial_differential_nonorth', lstat)
      IF lstat EQ -1 THEN BEGIN
        PRINT, "  -> Excessive work "+STR(f0)+" to "+STR(ftarget)+" Trying to continue..."
        lstat = 2 ; continue
        rcount = rcount + 1
        IF rcount GT 3 THEN BEGIN
          PRINT, "   -> Too many repeats. Giving Up."
          
          ri = lastgoodpos[0]
          zi = lastgoodpos[1]
          fmax = lastgoodf
          
          RETURN
        ENDIF
        ; Get f at this new location
        local_gradient, interp_data, rznew[0], rznew[1], status=status, f=f0
        IF status EQ 1 THEN BEGIN
          ri = ri0
          zi = zi0
          status = 1
          RETURN
        ENDIF
        rzold = rznew
  
        CONTINUE
      ENDIF ELSE BREAK
    ENDREP UNTIL 0
    
    IF lstat LT 0 THEN BEGIN
      PRINT, "Error in LSODE routine when following psi gradient."
      PRINT, "LSODE status: ", lstat
      ;STOP
    ENDIF
    CATCH, /CANCEL
    
    ri = rznew[0]
    zi = rznew[1]
    
  ENDIF ELSE BEGIN
    ; An error occurred in LSODE.
    ; lastgoodf contains the last known good f value
    ;PRINT, "CAUGHT ERROR "+!ERROR_STATE.MSG
    CATCH, /cancel
    ri = lastgoodpos[0]
    zi = lastgoodpos[1]
    fmax = lastgoodf
    IF ood THEN BEGIN
      ; Gone Out Of Domain
      status = 2
      fbndry = lastgoodf
      ;PRINT, "Out of domain at f = ", fbndry
      ; Repeat to verify that this does work
      rzold = [ri0, zi0]
      CATCH, theError
      fbndry = lastgoodf - 0.1D*(ftarget - f0)
      IF theError NE 0 THEN BEGIN
        PRINT, "   Error again at ", fbndry
      ENDIF
      rznew = LSODE(rzold,f0,fbndry - f0,'radial_differential')
      CATCH, /cancel
      
      RETURN
    ENDIF
    ; Otherwise just crossed a boundary
  ENDELSE
  CATCH, /cancel
  
  IF KEYWORD_SET(boundary) THEN BEGIN
    ; Check if the line crossed a boundary
    ;PRINT, "Checking boundary ", boundary[*,1:2], [ri0, ri], [zi0, zi]
    cpos = line_crossings([ri0, ri], [zi0, zi], 0, $
                          boundary[0,*], boundary[1,*], bndry_periodic, ncross=ncross, inds2=inds2)
    IF (ncross MOD 2) EQ 1 THEN BEGIN ; Odd number of boundary crossings
      
      IF SQRT( (ri - cpos[0,0])^2 + (zi - cpos[1,0])^2 ) GT 0.1D THEN BEGIN
        ;PRINT, "FINDING BOUNDARY", SQRT( (ri - cpos[0,0])^2 + (zi - cpos[1,0])^2 )
        ; Use divide-and-conquer to find crossing point
        
        tol = 1e-4 ; Make the boundary crossing stricter
        
        ibndry = inds2[0] ; Index in boundary where hit
        
        fcur = f0 ; Current known good position
        rzold = [ri0,zi0]
        rzcur = rzold
        REPEAT BEGIN
          fbndry = (fcur + fmax) / 2
          ; Try to go half-way to fmax
          CATCH, theError
          IF theError NE 0 THEN BEGIN
            ; Crossed boundary. Change fmax
            CATCH, /cancel
            fmax = fbndry
            ibndry = inds2[0] ; refined boundary index
          ENDIF ELSE BEGIN
            rznew = LSODE(rzold,f0,fbndry-f0,'radial_differential')
            
            ; Didn't cross. Make this the new current location
            fcur = fbndry
            rzcur = rznew
            
            CATCH, /cancel
          ENDELSE
          
        ENDREP UNTIL ABS(fmax - fcur) LT 0.01D*ABS(ftarget - f0)
        ri = rzcur[0]
        zi = rzcur[1]
        fbndry = fcur
        
        ;PRINT, "Hit boundary", ri, zi, " f =", fbndry
        status = 2
        RETURN
      ENDIF
    ENDIF
  ENDIF

  status = 0
END
