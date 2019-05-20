; IDL routine to read a 'g' file neqdsk
;
; Example:
;   g = read_neqdsk("neqdsk")
; 
; See comments at end for meaning of structure elements
;
; Format of G-EQDSK file is specified here:
;   https://fusion.gat.com/THEORY/efit/g_eqdsk.html
; 
; Ben Dudson, University of York, Feb 2010

FUNCTION read_neqdsk, file
  CATCH, theError
  IF theError THEN BEGIN
    CATCH, /cancel
    PRINT, "ERROR Occurred whilst reading G-EQDSK file"
    RETURN, 0
  ENDIF
  
  OPENR, fid, file, /GET_LUN, error=errid
  IF errid THEN BEGIN
    PRINT, "Error whilst reading '"+file+"' :" + !ERROR_STATE.MSG
    RETURN, 0
  ENDIF
  
  line=' '
  READF, fid, line  ; Read the first line
  s = STRSPLIT(line, ' ', /extract)  ; Strings separated by spaces
  ns = N_ELEMENTS(s)
  IF ns LT 3 THEN BEGIN
    PRINT, "ERROR: Expecting at least 3 numbers on first line"
    RETURN, 0
  ENDIF
  
  idum = 0
  nxefit = 0
  nyefit = 0
  READS, s[ns-3], idum
  READS, s[ns-2], nxefit
  READS, s[ns-1], nyefit
  
  PRINT, "   nxefit = ", nxefit, " nyefit = ", nyefit

  ; Set up the generator
  status = next_double(fid=fid)

  xdim   = next_double()
  zdim   = next_double()
  rcentr = next_double()
  rgrid1 = next_double()
  zmid   = next_double()
  
  rmagx  = next_double()
  zmagx  = next_double()
  simagx = next_double()
  sibdry = next_double()
  bcentr = next_double()
  
  cpasma = next_double()
  simagx = next_double()
  xdum   = next_double()
  rmagx  = next_double()
  xdum   = next_double()
  
  zmagx  = next_double()
  xdum   = next_double()
  sibdry = next_double()
  xdum   = next_double()
  xdum   = next_double()
  
  ; Read arrays
  fpol   = read_1d(nxefit)
  pres   = read_1d(nxefit)
  workk1 = read_1d(nxefit)
  workk2 = read_1d(nxefit)
  fold   = read_2d(nxefit, nyefit)
  qpsi   = read_1d(nxefit)
  
  nbdry  = LONG(next_double())
  nlim   = LONG(next_double())
  
  PRINT, nbdry, nlim
  
  IF nbdry GT 0 THEN BEGIN
    rbdry = DBLARR(nbdry)
    zbdry = DBLARR(nbdry)
    FOR i=0, nbdry-1 DO BEGIN
      rbdry[i] = next_double()
      zbdry[i] = next_double()
    ENDFOR
  ENDIF ELSE BEGIN
    rbdry = [0]
    zbdry = [0]
  ENDELSE
  
  IF nlim GT 0 THEN BEGIN
    xlim = DBLARR(nlim)
    ylim = DBLARR(nlim)
    FOR i=0, nlim-1 DO BEGIN
      xlim[i] = next_double()
      ylim[i] = next_double()
    ENDFOR
  ENDIF ELSE BEGIN
    xlim = [0]
    ylim = [0]
  ENDELSE
  
  FREE_LUN, fid

  ; Reconstruct the (R,Z) mesh
  r=DBLARR(nxefit, nyefit)
  z=DBLARR(nxefit, nyefit)
  
  FOR i=0,nxefit-1 DO BEGIN
    FOR j=0,nyefit-1 DO BEGIN
      r[i,j] = rgrid1 + xdim*i/(nxefit-1)
      z[i,j] = (zmid-0.5D*zdim) + zdim*j/(nyefit-1)
    ENDFOR
  ENDFOR

  ; Put data into structure
  
  result = {nx:nxefit, ny:nyefit, $ ; Number of horizontal and vertical points
            r:r, z:z, $             ; Location of the grid-points
            xdim:xdim, zdim:zdim, $ ; Size of the domain in meters
            rcentr:rcentr, bcentr:bcentr, $ ; Reference vacuum toroidal field (m, T)
            rgrid1:rgrid1, $ ; R of left side of domain
            zmid:zmid, $     ; Z at the middle of the domain
            rmagx:rmagx, zmagx:zmagx, $ ; Location of magnetic axis
            simagx:simagx, $ ; Poloidal flux at the axis (Weber / rad)
            sibdry:sibdry, $ ; Poloidal flux at plasma boundary (Weber / rad)
            cpasma:cpasma, $
            psi:fold, $   ; Poloidal flux in Weber/rad on grid points
            fpol:fpol, $  ; Poloidal current function on uniform flux grid
            pres:pres, $  ; Plasma pressure in nt/m^2 on uniform flux grid
            qpsi:qpsi, $  ; q values on uniform flux grid
            nbdry:nbdry, rbdry:rbdry, zbdry:zbdry, $  ; Plasma boundary
            nlim:nlim, xlim:xlim, ylim:ylim} ; Wall boundary
  CATCH, /cancel
  RETURN, result
END
