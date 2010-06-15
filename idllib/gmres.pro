;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; General Minimal Residual (GMRES)
;
; Inverts a linear operator, given an initial guess x0
; and RHS vector b
; 
; The operator should take a vector x and return operator(x) of the
; same length
;
; optional arguments:
;
; restart   Number of iterations between restarts
; max_iter  Maximum number of iterations
; tol       Tolerance
; stats     Returns statistics (errors etc)
; show      if set, plots a graph of error against iteration
; output    if set, regularly outputs the current error
;
;
; Ben Dudson, University of York, March 2008
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION dot, a, b
  RETURN, TRANSPOSE(a) # b
END

PRO update, x, it, H, s, v
  y = s
  
  ; backsolve
  FOR i=it, 0, -1 DO BEGIN
      y[i] = y[i] / H[i,i]
      FOR j=i-1, 0, -1 DO BEGIN
          y[j] = y[j] - H[j,i]*y[i]
      ENDFOR
  ENDFOR

  FOR i=0, it DO BEGIN
      x = x + REFORM(v[i,*]) * y[i]
  ENDFOR
END

FUNCTION gmres, x0, operator, b, restart=restart, max_iter=max_iter, tol=tol, $
                stats=stats, show=show, output=output
  
  x = x0 ; make a copy

  n = N_ELEMENTS(x) ; size of the problem

  IF NOT KEYWORD_SET(restart) THEN restart = n
  IF NOT KEYWORD_SET(max_iter) THEN max_iter = n
  IF NOT KEYWORD_SET(tol) THEN tol = 1.0e-4

  ; allocate memory

  s = FLTARR(restart+1)
  cs = s
  sn = s
  v = FLTARR(restart+1, n)
  H = FLTARR(restart+1, restart+1)

  normb = NORM(b)

  IF ABS(normb) LT tol THEN normb = 1.0

  ; r = b - Ax

  r = b - CALL_FUNCTION(operator, x)

  beta = NORM(r)

  resid = beta / normb
  IF resid LE tol THEN BEGIN
      stats = {iter:0, resid:resid}
      RETURN, x
  ENDIF


  outint = FIX(n/50)
  IF KEYWORD_SET(show) THEN BEGIN
      iarr = [0]
      rarr = [resid]
  ENDIF

  it = 1

  WHILE it LE max_iter DO BEGIN
      v[0,*] = r / beta
      
      s[0] = beta

      itt = 0
      WHILE (itt LT restart) AND (it LE max_iter) DO BEGIN

          ; w = A*v[itt]
          w = CALL_FUNCTION(operator, REFORM(v[itt,*]))

          FOR n=0, itt DO BEGIN
              H[n, itt] = dot(w, REFORM(v[n,*]))
              ; w = w - H[n][itt] * v[n]
              w = w - H[n,itt] * v[n,*]
          ENDFOR
          ; H[itt+1][itt] = |w|
          H[itt+1,itt] = NORM(w)

          ; v[itt+1] = w / |w|
          v[itt+1,*] = w / H[itt+1,itt]

          ; apply plane rotations
          
          FOR n=0, itt-1 DO BEGIN
              tmp = H[n, itt]
              H[n,itt]   = cs[n] * tmp        + sn[n] * H[n+1,itt]
              H[n+1,itt] = cs[n] * H[n+1,itt] - sn[n] * tmp
          ENDFOR
          
          ; generate plane rotation
          
          IF H[itt+1,itt] EQ 0.0 THEN BEGIN
              cs[itt] = 1.0
              sn[itt] = 0.0;
          ENDIF ELSE IF ABS(H[itt+1,itt]) GT ABS(H[itt,itt]) THEN BEGIN
              tmp = H[itt,itt] / H[itt+1,itt]
              sn[itt] = 1.0 / SQRT(1.0 + tmp^2)
              cs[itt] = tmp * sn[itt]
          ENDIF ELSE BEGIN
              tmp = H[itt+1,itt] / H[itt,itt]
              cs[itt] = 1.0 / SQRT(1.0 + tmp^2)
              sn[itt] = tmp * cs[itt]
          ENDELSE

          tmp = H[itt,itt]
          H[itt,itt]   = cs[itt] * tmp  + sn[itt] * H[itt+1,itt]
          H[itt+1,itt] = cs[itt] * H[itt+1,itt] - sn[itt] * tmp

          s[itt+1] = -sn[itt] * s[itt]
          s[itt] = cs[itt] * s[itt]
       
          resid = ABS(s[itt+1] / normb)

          IF it MOD outint EQ 0 THEN BEGIN
              IF KEYWORD_SET(output) THEN PRINT, it, itt, resid
              IF KEYWORD_SET(show) THEN BEGIN
                  iarr = [iarr, it]
                  rarr = [rarr, resid]
                  !P.multi=[0,0,1,0,0]
                  plot, iarr, rarr
              ENDIF
          ENDIF

          IF resid LT tol THEN BEGIN
              update, x, itt, H, s, v

              stats = {iter:it, resid:resid}
              RETURN, x
          ENDIF
          
          it = it + 1
          itt = itt + 1
      ENDWHILE

      update, x, itt-1, H, s, v

      ; r = b - Ax
      r = b - CALL_FUNCTION(operator, x)

      beta = NORM(r)
      
      resid = beta / normb
      
      ;PRINT, "=>", it, resid

      IF resid LE tol THEN BEGIN
          stats = {iter:it, resid:resid}
          RETURN, x
      ENDIF
  ENDWHILE
  
  PRINT, "WARNING: Maximum iterations reached"
  stats = {iter:it, resid:resid}
  RETURN, x
END

FUNCTION gmfunc, x
  COMMON gmres_test, M
  RETURN, M # x
END


PRO test_gmres
  COMMON gmres_test, M
  
  n = 10

  ; generate a random matrix

  M = RANDOMN(2378461, n,n)

  x = RANDOMN(2813, n)

  b = M # x

  x0 = FLTARR(n)

  x1 = GMRES(x0, "gmfunc", b)

  PRINT, MAX(x - x1)

  STOP
END
