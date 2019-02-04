;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; BOUT / BOUT++ input grid file generator                         ;
;                                                                 ;
; takes a PDB file generated from a GATO or ELITE input file      ;
; (gato2pdb or elite2pdb), and generates output in various forms, ;
; either for input to UEDGE, or directly to BOUT/BOUT++           ;
;                                                                 ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; interpolate a function on periodic domain. x is an array of indices
; (floating point)
FUNCTION loop_interp, f, x
  nx = N_ELEMENTS(x)

  ff = FFT(f)
  result = DBLARR(nx)
  
  FOR i=0, nx-1 DO BEGIN
      result[i] = fft_interp(ff, x[i])
  ENDFOR

  RETURN, result
END

function pdiff_rz, r, z, f
;
; Calculate partial derivatives df/dr and df/dz
; for function f given on a set of 5 points (r,z)
; using singular-value decomposition
;
; Inputs: arrays r[5],z[5],f[5]
; Output: structure {r:df/dr, z:df/dz}
;-------------------------------------------------

   n = N_ELEMENTS(r)

   A=TRANSPOSE([[fltarr(n)+1],[r-r(0)],[z-z(0)]])

   SVDC, A,W,U,V

   res=SVSOL(U,W,V,f)

   pdiff={r:res[1],z:res[2],phi:0.0}

return, pdiff
end

; returns structure index if exists, otherwise -1
FUNCTION var_present, st, name
  l = TAG_NAMES(st)
  n = N_ELEMENTS(l)

  FOR i=0, n-1 DO BEGIN
      IF STRCMP(name, l[i], /fold_case) THEN RETURN, i
  ENDFOR
  
  RETURN, -1
END

FUNCTION apply_transform, var, transform
  s = SIZE(var, /dim)
  s2 = SIZE(transform, /dim)
  
  IF N_ELEMENTS(s) NE 2 THEN BEGIN
      PRINT, "Error in apply_transform: var must be 2D"
      STOP
  ENDIF
  IF N_ELEMENTS(s2) NE 2 THEN BEGIN
      PRINT, "Error in apply_transform: trasform must be 2D"
      STOP
  ENDIF

  nx = s[0]
  IF s2[0] NE nx THEN BEGIN
      PRINT, "Error in apply_transform: number of radial points not preserved"
      STOP
  ENDIF

  ny = s2[1]

  newvar = DBLARR(nx, ny)

  FOR i=0, nx-1 DO BEGIN
      newvar[i,*]  = loop_interp(var[i,*], transform[i,*])
  ENDFOR

  RETURN, newvar
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; UEDGE style output
; Creates a set of grid files in the same format as UEDGE
;
; This set of routines specific to UEDGE output
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; takes a variable f[nx, ny] and produces cells like uedge f[ny2, nx, 5]
; transform gives the indices of the centre of the cells
FUNCTION generate_uedge, f, transform
   s = SIZE(f, /dimensions)
   nx = s[0]
   ny = s[1]

   s = SIZE(transform, /dimensions)
   IF nx NE s[0] THEN BEGIN
       PRINT, "Error: X size of f and transform do not match"
       STOP
   ENDIF
   ny2 = s[1]

   ; perform poloidal FFT
   FF = COMPLEXARR(nx, ny)
   FOR i=0, nx-1 DO BEGIN
       FF[i,*] = FFT(f[i,*])
   ENDFOR

   result = DBLARR(ny2+2, nx, 5)

   ; Go through bulk of  points (not guard cells)

   FOR i=1, nx-2 DO BEGIN
       FOR j=0, ny2-1 DO BEGIN
           ; centres
           result[j+1, i, 0] = fft_interp(FF[i,*], transform[i,j])

           ;IF result[j+1,i,0] LT 2.3 THEN BEGIN
           ;    PRINT, "Here1"
           ;    STOP
           ;ENDIF

           ; index 1: psi minus, theta minus
           ; get index half-way between this point and the last
           IF j EQ 0 THEN BEGIN
               tm1 = (transform[i-1,0] + transform[i-1, ny2-1] - ny)/2.0
               tm2 = (transform[i,0] + transform[i, ny2-1] - ny)/2.0
           ENDIF ELSE BEGIN
               tm1 = (transform[i-1,j] + transform[i-1,j-1])/2.0
               tm2 = (transform[i,j] + transform[i,j-1])/2.0
           ENDELSE
           result[j+1,i,1] = (fft_interp(FF[i-1,*], tm1) + fft_interp(FF[i,*], tm2))/2.0

           ; this point belongs to 3 other cells
           IF j EQ 0 THEN BEGIN
               jm = ny2-1
           ENDIF ELSE BEGIN
               jm = j - 1
           ENDELSE
           result[jm+1, i  , 2] = result[j+1,i,1]
           result[j+1 , i-1, 3] = result[j+1,i,1]
           result[jm+1, i-1, 4] = result[j+1,i,1]
       ENDFOR       
   ENDFOR

   ; Radial guard cells. These 1/2 size of normal cells

   FOR j=0, ny2-1 DO BEGIN
       ; inner guard cells
       
       result[j+1, 0, 0] = 0.25*result[j+1,1,0] + 0.75*fft_interp(FF[0,*], transform[0,j])
       
       ;IF result[j+1,0,0] LT 2.3 THEN BEGIN
       ;    PRINT, "Here1"
       ;    STOP
       ;ENDIF

       IF j EQ 0 THEN BEGIN
           tm = (transform[0,0] + transform[0, ny2-1] - ny)/2.0
       ENDIF ELSE BEGIN
           tm = (transform[0,j] + transform[0,j-1])/2.0
       ENDELSE
       result[j+1, 0, 1] = fft_interp(FF[0,*], tm)
       
       IF j EQ 0 THEN BEGIN
           jm = ny2-1
       ENDIF ELSE BEGIN
           jm = j - 1
       ENDELSE
       ; only need to set one other corner
       result[jm+1, 0, 2] = result[j+1,0,1]
       
       
       ; outer guard cells
       
       i = nx-1
       result[j+1,i,0] = 0.25*result[j+1,i-1,0] + 0.75*fft_interp(FF[i,*], transform[i,j])
       
       ;IF result[j+1,i,0] LT 2.3 THEN BEGIN
       ;    PRINT, "Here1"
       ;    STOP
       ;ENDIF

       IF j EQ 0 THEN BEGIN
           tm1 = (transform[i-1,0] + transform[i-1, ny2-1] - ny)/2.0
           tm2 = (transform[i,0] + transform[i, ny2-1] - ny)/2.0
       ENDIF ELSE BEGIN
           tm1 = (transform[i-1,j] + transform[i-1,j-1])/2.0
           tm2 = (transform[i,j] + transform[i,j-1])/2.0
       ENDELSE
       result[j+1,i,3] = fft_interp(FF[i,*], tm2);
       result[j+1,i,1] = (fft_interp(FF[i-1,*], tm1) + result[j+1,i,3])/2.0
       
       ; set other points
       result[jm+1,i, 2] = result[j+1,i,1]
       result[jm+1,i, 4] = result[j+1,i,3]

       result[j+1,i-1,3] = result[j+1,i,1]
       result[jm+1,i-1,4] = result[j+1,i,1]
   ENDFOR

   RETURN, result
END

PRO add_poloidal_guards, f
   ; Add poloidal guard cells. These have zero poloidal size

   s = SIZE(f, /dimensions)
   ny = s[0]
   nx = s[1]

   FOR i=0, nx-1 DO BEGIN
       f[0,i,1] = f[1,i,1]
       f[0,i,2] = f[1,i,1]
       f[0,i,3] = f[1,i,3]
       f[0,i,4] = f[1,i,3]
       f[0,i,0] = (f[0,i,1] + f[0,i,3])/2.0
       
       f[ny-1,i,1] = f[ny-2,i,2]
       f[ny-1,i,2] = f[ny-2,i,2]
       f[ny-1,i,3] = f[ny-2,i,4]
       f[ny-1,i,4] = f[ny-2,i,4]
       f[ny-1,i,0] = (f[ny-2,i,2] + f[ny-2,i,4])/2.0
   ENDFOR

END

; re-arranges f[y,x,5] so that the y index is moved to zero
; need to keep the guard cells where they are
FUNCTION rearrange, f, index
  s = SIZE(f, /dimensions)
  ny = s[0]
  nx = s[1]
  
  IF index EQ 0 THEN RETURN, f ; nothing to do
  
  result = DBLARR(ny, nx, 5)
  FOR i=0, 4 DO BEGIN
      FOR j=1, ny-2 DO BEGIN
          jold = 1 + ((j + index - 1) MOD (ny-2))
          result[j, *, i] = f[jold, *, i]
      ENDFOR
  ENDFOR

  RETURN, result
END

; main UEDGE output file generator
; produce 3 PDB output files like UEDGE
PRO save_uedge, a, min_x, max_x, edge, transform
   PRINT, "Generating cells for uedge-style output (may take some time)"

   nx = max_x - min_x + 1 ; size of the output array

   s = SIZE(transform, /dim)
   IF s[0] NE nx THEN BEGIN
       PRINT, "Error: X size of transform array not consistent with x range"
       RETURN
   ENDIF
   ny = s[1]  ; number of output poloidal points
   
   psixy = DBLARR(nx, a.ny)
   FOR i=0, nx-1 DO BEGIN
       psixy[i,*] = a.psi[min_x + i]
   ENDFOR

   ; quantities on cell centres and corners
   rm   = generate_uedge(a.Rxy[min_x:max_x, *], transform)
   PRINT, "R"
   zm   = generate_uedge(a.Zxy[min_x:max_x, *], transform)
   PRINT, "Z"
   psi  = generate_uedge(psixy, transform)
   PRINT, "Psi"
   br   = generate_uedge(a.BRxy[min_x:max_x, *], transform)
   PRINT, "Br"
   bz   = generate_uedge(a.BZxy[min_x:max_x, *], transform)
   PRINT, "Bz"
   bpol = generate_uedge(a.BPxy[min_x:max_x, *], transform)
   PRINT, "Bpol"
   bphi = generate_uedge(a.BTxy[min_x:max_x, *], transform)
   PRINT, "Btor"
   b    = generate_uedge(a.Bxy[min_x:max_x, *], transform)
   PRINT, "B"
   jparue = generate_uedge(a.Jpar[min_x:max_x, *], transform)
   PRINT, "Jpar"

   ; find the inboard midplane index
   
   m = MIN(rm[*, nx-1, 0], ind)

   IF ind NE 0 THEN BEGIN
       PRINT, "Re-arranging grid to have origin on inner midplane, currently index:",ind 

       rm   = rearrange(rm, ind)
       zm   = rearrange(zm, ind)
       psi  = rearrange(psi, ind)
       br   = rearrange(br, ind)
       bz   = rearrange(bz, ind)
       bpol = rearrange(bpol, ind)
       bphi = rearrange(bphi, ind)
       b    = rearrange(b, ind)
       jparue = rearrange(jparue, ind)
   ENDIF

   ; now add poloidal guard cells
   add_poloidal_guards, rm
   add_poloidal_guards, zm
   add_poloidal_guards, psi
   add_poloidal_guards, br
   add_poloidal_guards, bz
   add_poloidal_guards, bpol
   add_poloidal_guards, bphi
   add_poloidal_guards, b
   
   ; quantities on cell centres only
   ni  = fltarr(ny+2, nx)
   te  = ni
   ti  = ni
   up  = ni
   phi = ni
   dx  = ni
   jpar0 = REFORM(jparue[*,*,0]) ; just take central values
      
   FOR j=0, ny+1 DO BEGIN
       ;-calculate radial diameter of cell
       rbot=0.5*(rm[j,*,1]+rm[j,*,2])
       zbot=0.5*(zm[j,*,1]+zm[j,*,2])
       rtop=0.5*(rm[j,*,3]+rm[j,*,4])
       ztop=0.5*(zm[j,*,3]+zm[j,*,4])
       dx[j,*]=sqrt((rtop-rbot)^2+(ztop-zbot)^2)
   ENDFOR
      
   PRINT, "========== SET PLASMA PARAMETERS ========="

   pressure = a.mu0p[min_x:max_x] ;/ (4.0*!PI*1.0e-7)  ; get pressure

   REPEAT BEGIN
       Temp = get_float("Constant temperature (eV): ")
       density = pressure / (2.0*Temp * 1.602e-19) ; /2 to get ion (or electron) density
       
       PRINT, "Density range: ", MIN(density), MAX(density)
       done = get_yesno("Is this ok? (no = change temperature)")
   ENDREP UNTIL done EQ 1
   
   min_density = get_float("Minimum density:")
   
   te[*,*] = Temp
   ti[*,*] = Temp
   FOR j=0, ny+1 DO BEGIN
       ni[j,*] = density + min_density - MIN(density)
   ENDFOR
   
   ; output is MKS like UEDGE. GATO file should also be MKS

   PRINT, "========== WRITING OUTPUT FILES =========="

   d1={$
        rm_com:rm,$
        zm_com:zm,$
        psi_com:psi,$
        br_com:br,$
        bz_com:bz,$
        bpol_com:bpol,$
        bphi_com:bphi,$
        b_com:b,$
        nx_com:ny,$
        ny_com:nx-2,$
        ixpt1_com:0,$
        ixpt2_com:ny,$
        nxpt_com:0,$
        ixlb_com:0,$
        ixrb_com:ny,$
        iysptrx1_com:nx-2,$ 
        iysptrx2_com:nx-2,$
        ix_lim_com:0,$
        runidg_grd:a.file}
   
   d2={$
        gy_com:1./dx,$ 
        xcs_com:0.0,$           ;-never used
        yyc_com:rm[*,0],$ 
        rmagx_flx:0.0,$         ;-never used       
        rseps_flx:0.0,$         ;-never used 
        zmagx_flx:0.0,$         ;-never used
        zseps_flx:0.0,$         ;-never used
        psi0max_flx:0.0,$       ;-never used
        psi0min1_flx:0.0,$      ;-never used
        psi0min2_flx:0.0,$      ;-never used
        sibdryg_grd:a.psi[edge],$
        simagxg_grd:a.psi[0]}  

   d3={$
        ni___1__value:ni,$
        up___1__value:up,$
        te____ev_value:te,$
        ti____ev_value:ti,$
        jpar____value:jpar0, $
        phi____value:phi}   
   
   file1='gridue.pdb' 
   file2='uedgegrd.pdb'
   file3='uedgeout.pdb'
      
   ; need to ensure files don't already exist (for some reason)
   SPAWN, "rm -f "+file1
   SPAWN, "rm -f "+file2
   SPAWN, "rm -f "+file3
   
   PRINT, 'Writing files...'
   PRINT, file1  
   pd_export, file1, d1

   PRINT, file2  
   pd_export, file2, d2
   
   PRINT, file3        
   pd_export, file3, d3

   PRINT, "Finished writing UEDGE-style output"
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                Main PDB2BOUT routine                         ;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


PRO pdb2bout, file, uedge=uedge, output=output, input=input, _extra=_extra

  IF NOT KEYWORD_SET(output) THEN output="bout.grd.pdb"

  IF KEYWORD_SET(output) THEN output=EXPAND_PATH(output)
  file = EXPAND_PATH(file)

  MU = 4.0*!PI*1.0e-7

  safe_colors, /first

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Load input file
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  PRINT, "Loading input file..."

  a = file_import(file)

  ; add the name of the file to the structure
  a = CREATE_STRUCT(a, "file", file)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Check the input file
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ; list of needed variables
  needed = ["nx", "ny", "psi", "mu0p", "mu0pprime", "f", "ffprime", "qsafe", "Rxy", "Zxy"]
  ; number of dimensions
  dims =   [  0 ,  0  ,   1  ,   1   ,      1     ,  1 ,    1     ,    1   ,   2  ,   2  ]

  FOR i=0, N_ELEMENTS(needed)-1 DO BEGIN
      ; check the variable is present

      sind = var_present(a, needed[i])
      IF sind EQ -1 THEN BEGIN
          ; structure does not contain variable
          PRINT, "Error: Input file must contain "+ needed[i]
          RETURN
      ENDIF

      ; check the number of dimensions
      nd = SIZE(a.(sind), /n_dim)
      
      IF nd NE dims[i] THEN BEGIN
          PRINT, "Error: Input variable "+needed[i]+" must have "+STRTRIM(STRING(dims[i]),2)+" dimensions"
          RETURN
      ENDIF

      ; check the sizes of the dimensions
      s = SIZE(a.(sind), /dim)
      IF nd EQ 1 THEN BEGIN
          IF s[0] NE a.nx THEN BEGIN
              PRINT, "Error: Input variable "+needed[i]+" must be of length nx"
              RETURN
          ENDIF
      ENDIF ELSE IF nd EQ 2 THEN BEGIN
          IF s[0] NE a.nx THEN BEGIN
              PRINT, "Error: First index of input "+needed[i]+" must be length nx"
              RETURN
          ENDIF
          IF s[1] NE a.ny THEN BEGIN
              PRINT, "Error: Second index of input "+needed[i]+" must be length ny"
              RETURN
          ENDIF
      ENDIF
  ENDFOR

  IF ( var_present(a, "Brxy") NE -1 ) OR (var_present(a, "Bzxy") NE -1) THEN BEGIN
      ; at least one of Brxy and Bzxy are given in the input file
      ; hence both of them should be
      
      IF ( var_present(a, "Brxy") EQ -1 ) OR (var_present(a, "Bzxy") EQ -1) THEN BEGIN
          PRINT, "Error: If one of Brxy or Bzxy are specified then the other must be"
          RETURN
      ENDIF
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Test if the last poloidal point is the
  ; same as the first
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  dl = (a.rxy[*,a.ny-1] - a.rxy[*,0])^2 + (a.zxy[*,a.ny-1] - a.zxy[*,0])
  IF MAX(dl) LT 1.0e-4 THEN BEGIN
      PRINT, "**Last poloidal point duplicates the first. Removing..."
      
      tn = TAG_NAMES(a)
      n = N_ELEMENTS(tn)
      
      j = 0
      FOR i=0, n-1 DO BEGIN
          nd = SIZE(a.(i), /n_dim)

          IF nd EQ 2 THEN BEGIN
              ; 2D array - remove the last point
              
              IF j EQ 0 THEN BEGIN
                  b = CREATE_STRUCT(tn[i], a.(i)[*,0:a.ny-2])
                  j = 1
              ENDIF ELSE  b = CREATE_STRUCT(b, tn[i], a.(i)[*,0:a.ny-2])
          ENDIF ELSE BEGIN
              ; just add it to the new structure
              IF j EQ 0 THEN BEGIN
                  b = CREATE_STRUCT(tn[i], a.(i))
                  j = 1
              ENDIF ELSE BEGIN
                  b = CREATE_STRUCT(b, tn[i], a.(i))
              ENDELSE
          ENDELSE
      ENDFOR

      b.ny = a.ny - 1
      
      a = b
  ENDIF

  ; get size of the array
  nx = a.nx
  ny = a.ny

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Calculate Br, Bz and Bpol as a check
  ; using local differentials of psi
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  PRINT, "Calculating poloidal field using SVD..."

  Brxy = DBLARR(nx, ny)
  Bzxy = DBLARR(nx, ny)
  
  FOR i=1, nx-2 DO BEGIN
      FOR j=0, ny-1 DO BEGIN
          IF j EQ 0 THEN BEGIN
              r = [a.rxy[i,j], a.rxy[i-1,j], a.rxy[i+1,j], a.rxy[i,ny-1], a.rxy[i,j+1]]
              z = [a.zxy[i,j], a.zxy[i-1,j], a.zxy[i+1,j], a.zxy[i,ny-1], a.zxy[i,j+1]]
              p = [a.psi[i], a.psi[i-1], a.psi[i+1], a.psi[i], a.psi[i]]
          ENDIF ELSE IF j EQ ny-1 THEN BEGIN
              r = [a.rxy[i,j], a.rxy[i-1,j], a.rxy[i+1,j], a.rxy[i,j-1], a.rxy[i,0]]
              z = [a.zxy[i,j], a.zxy[i-1,j], a.zxy[i+1,j], a.zxy[i,j-1], a.zxy[i,0]]
              p = [a.psi[i], a.psi[i-1], a.psi[i+1], a.psi[i], a.psi[i]]
          ENDIF ELSE BEGIN
              r = [a.rxy[i,j], a.rxy[i-1,j], a.rxy[i+1,j], a.rxy[i,j-1], a.rxy[i,j+1]]
              z = [a.zxy[i,j], a.zxy[i-1,j], a.zxy[i+1,j], a.zxy[i,j-1], a.zxy[i,j+1]]
              p = [a.psi[i], a.psi[i-1], a.psi[i+1], a.psi[i], a.psi[i]]
          ENDELSE
              
          diff = pdiff_rz(r, z, p)
          
          Brxy[i,j] = -1.0*diff.z / a.Rxy[i,j]
          Bzxy[i,j] = diff.r / a.Rxy[i,j]
      ENDFOR
  ENDFOR

  FOR j=0, ny-1 DO BEGIN
      jp = (j + 1) MOD ny
      jm = (j - 1 + ny) MOD ny
      
      ; first point
      i = 0

      r = [a.rxy[i,j], a.rxy[i,jm], a.rxy[i,jp], a.rxy[i+1,j], a.rxy[i+1,jm], a.rxy[i+1,jp]]
      z = [a.zxy[i,j], a.zxy[i,jm], a.zxy[i,jp], a.zxy[i+1,j], a.zxy[i+1,jm], a.zxy[i+1,jp]]
      p = [a.psi[i], a.psi[i], a.psi[i], a.psi[i+1], a.psi[i+1], a.psi[i+1]]

      diff = pdiff_rz(r, z, p)
          
      Brxy[i,j] = -1.0*diff.z / a.Rxy[i,j]
      Bzxy[i,j] = diff.r / a.Rxy[i,j]

      ; last point
      i = nx-1
      
      r = [a.rxy[i,j], a.rxy[i,jm], a.rxy[i,jp], a.rxy[i-1,j], a.rxy[i-1,jm], a.rxy[i-1,jp]]
      z = [a.zxy[i,j], a.zxy[i,jm], a.zxy[i,jp], a.zxy[i-1,j], a.zxy[i-1,jm], a.zxy[i-1,jp]]
      p = [a.psi[i], a.psi[i], a.psi[i], a.psi[i-1], a.psi[i-1], a.psi[i-1]]

      diff = pdiff_rz(r, z, p)
          
      Brxy[i,j] = -1.0*diff.z / a.Rxy[i,j]
      Bzxy[i,j] = diff.r / a.Rxy[i,j]
  ENDFOR

  Bpxy = SQRT(Brxy^2 + Bzxy^2)
  
  IF var_present(a, "Bpxy") NE -1 THEN BEGIN
      ; Bp supplied in the input file - keep, but compare methods
      PRINT, "Difference in Bpxy computed two different ways:", MAX(ABS(Bpxy[1:*,*] - a.Bpxy[1:*,*]))
      PRINT, "Maximum percentage difference: ", 100.0*MAX(ABS( (Bpxy[1:*,*] - a.Bpxy[1:*,*])/a.Bpxy[1:*,*] ))

      err = FLTARR(nx)
      FOR i=0, nx-1 DO BEGIN
          err[i] = 100.0*MAX(ABS( (Bpxy[i,*] - a.Bpxy[i,*])/a.Bpxy[i,*] ))
      ENDFOR
      !P.multi=[0,0,1,0,0]
      plot, a.psi, err, xtitle="Psi", ytitle="Maximum % Error in Bp", /ylog, $
        title="Error in SVD Bp calculation", color=1
      
      IF var_present(a, "Brxy") EQ -1 THEN BEGIN
          ; have Bp but not Br and Bz
          
          a = CREATE_STRUCT(a, "Brxy", Brxy, "Bzxy", Bzxy)
      ENDIF ELSE BEGIN
          ; already have Br and Bz
          
          PRINT, "Using given Bp, Br and Bz values"
      ENDELSE
      
      WAIT, 2

  ENDIF ELSE BEGIN
      ; Bp not supplied

      IF var_present(a, "Brxy") NE -1 THEN BEGIN
          ; calculate Bp from Br and Bz
          
          Bp2 = SQRT(a.Brxy^2 + a.Bzxy^2)
          
          PRINT, "Difference in Bpxy computed two different ways:", MAX(ABS(Bpxy - Bp2))
          PRINT, "Maximum percentage difference: ", 100.0*MAX(ABS( (Bpxy - Bp2)/Bp2 ))

          PRINT, "Using Bp from given Br and Bz"
          
          a = CREATE_STRUCT(a, "Bpxy", Bp2)
      ENDIF ELSE BEGIN
          ; don't have Brxy and Bzxy either
          PRINT, "Using Bp, Br and Bz from SVD calculation"

          a = CREATE_STRUCT(a, "Bpxy", Bpxy, "Brxy", Brxy, "Bzxy", Bzxy)
      ENDELSE
      
  ENDELSE

  ; overwrite given Bp, Br and Bz values
  ;a.Bpxy = Bpxy
  ;a.Brxy = Brxy
  ;a.Bzxy = Bzxy

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; structure 'a' now contains all poloidal B components
  ; Check if toroidal field is included
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  PRINT, "Calculating toroidal field from input f"

  Btxy = DBLARR(nx, ny)
  FOR i=0, ny-1 DO BEGIN
      Btxy[*,i] = a.f / a.Rxy[*,i]
  ENDFOR
  
  IF var_present(a, "Btxy") EQ -1 THEN BEGIN
      
      PRINT, "Using Bt from f / R"
      a = CREATE_STRUCT(a, "Btxy", Btxy)
  ENDIF ELSE BEGIN
      ; check values of Btxy
      
      PRINT, "Difference in Bt computed two different ways:", MAX(ABS(Btxy - a.Btxy))
      PRINT, "Maximum percentage difference: ", 100.0*MAX(ABS( (Btxy - a.Btxy)/a.Btxy ))

      PRINT, "Using Bt from input file"
  ENDELSE

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Get p and pprime from mu0p and mu0pprime
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; peak pressure should be in the range of 1e3 -> 1e4
  ; mu0 ~ 1e-6 -> mu0p should be ~ 1e-3 -> 1e-2
   
  IF MAX(a.mu0p) GT 1.0 THEN BEGIN
      PRINT, "***Maximum mu0p is ", MAX(a.mu0p)

      if keyword_set(INPUT) then ispressure=input.ispressure else $
        ispressure=get_yesno("Is this pressure (not mu0*pressure)?")

      IF ispressure EQ 0 THEN BEGIN
          a.mu0p = a.mu0p / MU
          a.mu0pprime = a.mu0pprime / MU
      ENDIF
      
  ENDIF ELSE BEGIN
      ; probably is mu0*pressure
      ; convert to pressure
      a.mu0p = a.mu0p / MU
      a.mu0pprime = a.mu0pprime / MU
  ENDELSE

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Get edge of the plasma
  ; Either where p = 0, or edge of grid
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  w = WHERE((a.mu0p - MIN(a.mu0p)) LT 0.1, count)
  
  IF count LE 1 THEN BEGIN
      PRINT, "Grid entirely inside plasma"
      edge = a.nx-1;
  ENDIF ELSE BEGIN
      PRINT, "Grid contains vacuum region"
      edge = w[0]-1;
  ENDELSE

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Get magnitude of B and parallel current
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; get field magnitude
  Bxy = SQRT(a.BPxy^2 + a.BTxy^2)
  
  ; calculate parallel current
  Jpar = DBLARR(nx, ny)
  FOR j=0, ny-1 DO BEGIN
      Jpar[*,j] = Bxy[*,j]*a.ffprime/(a.f*MU) + a.Rxy[*,j]*a.mu0pprime
  ENDFOR

  ; add to the a structure

  a = CREATE_STRUCT(a, "Bxy", Bxy, "Jpar", Jpar)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Plot the grid, pressure and safety factor
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  !P.MULTI=[0,2,0,0,0]

  plot, a.rxy, a.zxy, psym=3, /iso, color=1
  oplot, a.rxy[edge, *], a.zxy[edge, *], color=2, thick=1.5

  ; get psi normalised
  psi_min = a.psi[0]
  psi_bndry = a.psi[edge]

  psin = (a.psi - psi_min) / (psi_bndry - psi_min);

  PRINT, "PSI normalised range: 0 to "+STRTRIM(STRING(max(psin)),2)
  PRINT, "Number of radial grid points  : "+STRTRIM(STRING(a.nx),2)
  PRINT, "Number of poloidal grid points: "+STRTRIM(STRING(a.ny),2)
  
  !P.multi=[3,2,2,0,0]
  plot, psin, a.qsafe, title="Safety Factor", color=1
  oplot, [1.0,1.0], [-1.0, 1000.0], lines=2, color=2
  !P.multi=[1,2,2,0,0]
  plot, psin, a.mu0p, title="Pressure [Pa]", color=1
  oplot, [1.0,1.0], [-1.0, 1.0e6], lines=2, color=2

  PRINT, "Edge Q: "+STRTRIM(STRING(a.qsafe[edge]),2)
  PRINT, "Range of parallel current density [A/m^2]: ", min(jpar), max(jpar)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Get Normalisations consistent with ELITE
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  bmag = 0.5*(MIN(a.btxy[edge, *]) + MAX(a.btxy[edge, *]))
  rmag = 0.5*(MIN(a.rxy[edge, *]) + MAX(a.rxy[edge, *]))

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Get (approximate) region for grid from user
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  done = 0
  REPEAT BEGIN
      
      if keyword_set(INPUT) then begin
          min_psin = input.min_psin
          max_psin = input.max_psin
      endif else begin
          min_psin = get_float("Inner Psi boundary:")
          max_psin = get_float("Outer Psi boundary:")
      endelse


      PRINT, "Psi range ", min_psin, max_psin

      w = WHERE((psin GE min_psin) AND (psin LE max_psin),count)

      IF count LT 2 THEN BEGIN
          PRINT, "Region contains insufficient radial points"
      ENDIF ELSE BEGIN
          min_x = w[0]
          max_x = w[count-1]
          
          !P.multi=[0,2,0,0,0]
          
          plot, a.rxy, a.zxy, psym=3, /iso, color=1
          oplot, a.rxy[edge, *], a.zxy[edge, *], color=2, thick=1.5

          oplot, a.rxy[min_x,*], a.zxy[min_x, *], color=3, thick=1.5
          oplot, a.rxy[max_x,*], a.zxy[max_x, *], color=3, thick=1.5

          !P.multi=[3,2,2,0,0]
          plot, psin, a.qsafe, title="Safety Factor", color=1
          oplot, [1.0,1.0], [-1.0, 1000.0], lines=2, color=2
          oplot, [min_psin, min_psin], [-1.0, 1000.0], color=3
          oplot, [max_psin, max_psin], [-1.0, 1000.0], color=3
          !P.multi=[1,2,2,0,0]
          plot, psin, a.mu0p, title="Pressure [Pa]", color=1
          oplot, [1.0,1.0], [-1.0, 1.0e6], lines=2, color=2
          oplot, [min_psin, min_psin], [-1.0, 1.0e6], color=3
          oplot, [max_psin, max_psin], [-1.0, 1.0e6], color=3
          
          PRINT, "Number of radial grid points: ", count
          w2 = WHERE(w LE edge, count2)
          PRINT, "Of which inside the plasma: ", count2
          IF count2 LT 2 THEN BEGIN
              PRINT, "Insufficient points inside plasma"
          ENDIF ELSE BEGIN

              IF keyword_set(INPUT) then done=input.range_ok else $
                done = get_yesno("Is this range ok?")

          ENDELSE
      ENDELSE
  ENDREP UNTIL done EQ 1

  nx = count

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; GET PLASMA PROFILES
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  PRINT, "====== SETTING PLASMA PROFILES ======="
  
  gotprof = 0
  
  IF (var_present(a, "Te") NE -1) OR (var_present(a, "Ti") NE -1) OR (var_present(a, "Ni") NE -1) THEN BEGIN
      ; got some plasma parameters
      
      PRINT, "Some plasma parameters given in input file"

      IF keyword_set(INPUT) then use_given=input.use_given else $
        use_given=get_yesno("Use given parameters?")

      IF use_given EQ 1 THEN BEGIN
      
          IF (var_present(a, "Te") NE -1) OR (var_present(a, "Ti") NE -1) THEN BEGIN
              ; got one or both of Te or Ti
              
              IF var_present(a, "Te") EQ -1 THEN BEGIN
                  ; don't have Te
                  
                  PRINT, "Setting Te = Ti"
                  a = CREATE_STRUCT(a, "Te", a.Ti)
                  
              ENDIF ELSE IF var_present(a, "Ti") EQ -1 THEN BEGIN
                  PRINT, "Setting Ti = Te"
                  
                  a = CREATE_STRUCT(a, "Ti", a.Te)
              ENDIF ELSE BEGIN
                  PRINT, "Got both Ti and Te"
              ENDELSE
              
              ; set density to be consistent with pressure
              density = a.mu0p / ( (a.Te + a.Ti)*1.602e-19*1.0e20 )
              
              PRINT, "Setting density to be consistent with pressure"
              PRINT, "Maximum density (10^20 m^-3):", max(density)
              
              IF var_present(a, "Ni") NE -1 THEN BEGIN
                  PRINT, "   Density given in input: ", max(a.ni)
                  PRINT, "   Maximum difference: ", MAX(ABS(a.ni - density))
                  
                  a.ni = density
              ENDIF ELSE a = CREATE_STRUCT(a, "Ni", density)
              
          ENDIF ELSE BEGIN
              ; don't have any temperature data
              ; have got density - set Ti = Te to be consistent
                  
              PRINT, "Got density, setting temperature (Ti = Te)"
              
              a = CREATE_STRUCT(a, "Te", (a.mu0p / (2.*a.ni*1.602e-19*1.0e20)))
              a = CREATE_STRUCT(a, "Ti", a.Te)
              
              PRINT, "Maximum temperature (eV):", max(a.Te)
          ENDELSE
          gotprof = 1
      ENDIF
  ENDIF
  
  IF gotprof EQ 0 THEN BEGIN
      ; ignoring input plasma profiles (or none given)
             
      PRINT, "Generating plasma profiles:"
          
      PRINT, "  1. Flat temperature profile"
      PRINT, "  2. Flat density profile"
      PRINT, "  3. Te proportional to density"
      REPEAT BEGIN

          if keyword_set(INPUT) then opt=input.profile_option else $
            opt = get_integer("Profile option:")

      ENDREP UNTIL (opt GE 1) AND (opt LE 3)
      
      IF opt EQ 1 THEN BEGIN
                                ; flat temperature profile
          
          PRINT, "Setting flat temperature profile"
          REPEAT BEGIN

              if keyword_set(INPUT) then Te_x=input.te_x else $
                Te_x = get_float("Temperature (eV):")
                            
                                ; get density
              Ni = a.mu0p / (2.*Te_x* 1.602e-19*1.0e20)
              
              PRINT, "Maximum density (10^20 m^-3):", MAX(Ni)

              if keyword_set(INPUT) then done=input.max_dens_ok else $
                done = get_yesno("Is this ok?")

          ENDREP UNTIL done EQ 1
          
          Te = FLTARR(a.nx)+Te_x
          Ti = Te
          
      ENDIF ELSE IF opt EQ 2 THEN BEGIN
          PRINT, "Setting flat density profile"
          
          REPEAT BEGIN
              ni_x = get_float("Density [10^20 m^-3]:")
              
                                ; get temperature
              Te = a.mu0p / (2.*ni_x* 1.602e-19*1.0e20)
              
              PRINT, "Maximum temperature (eV):", MAX(te)
          ENDREP UNTIL get_yesno("Is this ok?") EQ 1
          
          Ti = Te
          Ni = FLTARR(a.nx) + ni_x
      ENDIF ELSE BEGIN
          PRINT, "Setting te proportional to density"
          
          REPEAT BEGIN
              te_x = get_float("Maximum temperature [eV]:")
              
              ni_x = max(a.mu0p) / (2.*Te_x* 1.602e-19*1.0e20)
              
              PRINT, "Maximum density [10^20 m^-3]:", ni_x

              shape = sqrt(a.mu0p / max(a.mu0p))
              Te = te_x * shape
              Ni = ni_x * shape
          ENDREP UNTIL get_yesno("Is this ok?") EQ 1
          Ti = Te
      ENDELSE

      ; put into the structure
      IF var_present(a, "Te") EQ -1 THEN a = CREATE_STRUCT(a, "Te", Te) ELSE a.Te = Te
      IF var_present(a, "Ti") EQ -1 THEN a = CREATE_STRUCT(a, "Ti", Ti) ELSE a.Ti = Ti
      IF var_present(a, "Ni") EQ -1 THEN a = CREATE_STRUCT(a, "Ni", Ni) ELSE a.Ni = Ni
  ENDIF
  
  !P.multi=[0,2,0,0,0]

  plot, psin, a.ni, title="Density [10^20 m^-3]", xtitle="Normalised Psi", color=1
  plot, psin, a.te, title="Te and Ti [eV]", xtitle="Normalised Psi", $
    yrange=[0.0, MAX([a.te, a.ti])], color=1
  ;cursor, x, y, /down

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Extract data from structure
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; 1D arrays
  pressure = a.mu0p[min_x:max_x]
  dpdpsi = a.mu0pprime[min_x:max_x]

  ; 2D arrays
  Rxy = a.Rxy[min_x:max_x, *]
  Zxy = a.Zxy[min_x:max_x, *]
  Jpar = Jpar[min_x:max_x, *]
  Bxy = a.Bxy[min_x:max_x, *]
  Bpxy = a.Bpxy[min_x:max_x, *]
  Btxy = a.Btxy[min_x:max_x, *]

  qsafe = a.qsafe[min_x:max_x]

  Te = DBLARR(nx, a.ny)
  Ti = Te
  Ni = Te

  psixy = Te
  
  FOR i=0, nx-1 DO BEGIN
      Te[i,*] = a.Te[i+min_x]
      Ti[i,*] = a.Ti[i+min_x]
      Ni[i,*] = a.Ni[i+min_x]
      
      psixy[i,*] = a.psi[i+min_x]
  ENDFOR

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Increase radial resolution (optional)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if keyword_set(INPUT) then increase_radial=input.increase_radial else $
    increase_radial=get_yesno("Increase radial resolution?")

  IF increase_radial THEN BEGIN

      if keyword_set(INPUT) then nx2=input.nx2 else $
        nx2 = get_integer("Number of radial points:")

      IF KEYWORD_SET(smooth) THEN BEGIN
        PRINT, "Adding "+STRTRIM(STRING(2*width),2)+" X points for smoothing"
        nx2 = nx2 + 2*width
      ENDIF
          
      Rxy = xinterp(Rxy, nx2)
      Zxy = xinterp(Zxy, nx2)
      Bxy = xinterp(Bxy, nx2)
      Bpxy = xinterp(Bpxy, nx2)
      Btxy = xinterp(Btxy, nx2)
      Jpar = xinterp(Jpar, nx2)

      qsafe = xinterp(qsafe, nx2)

      Te = xinterp(Te, nx2) > 0.0
      Ti = xinterp(Ti, nx2) > 0.0
      Ni = xinterp(Ni, nx2) > 0.0

      psixy = xinterp(psixy,nx2)
      
      dpdpsi = xinterp(dpdpsi, nx2)
      pressure = xinterp(pressure, nx2)
      
      nx = nx2

      plot, Ni[*,0], title="Density [10^20 m^-3]", xtitle="Index", color=1
      plot, Te[*,0], title="Te [eV]", xtitle="Index", color=1
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; GENERATE ORTHOGONAL COORDINATES
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  PRINT, "======== GENERATING ORTHOGONAL COORDINATES FOR BOUT ========= "

  REPEAT BEGIN

      if keyword_set(INPUT) then ny=input.ny else $
        ny = get_integer("Number of poloidal grid points:")
      
      REPEAT BEGIN

          if keyword_set(INPUT) then index=input.index else $
            index = get_integer("Enter x index of equal hthe [0, "+$
                                STRTRIM(STRING(nx-1),2)+"] :")

      ENDREP UNTIL (index GE 0) AND (index LT nx)
      
      ; create orthogonal coordinate system

      transform = gen_orthog2(Rxy, Zxy, index, n=ny)
      transform = transform * DOUBLE(a.ny) / (2.0*!PI)
      
      ; generate new Rxy and Zxy

      PRINT, "Interpolating Rxy"
      Rxy2 = apply_transform(Rxy, transform)
      PRINT, "Interpolating Zxy"
      Zxy2 = apply_transform(Zxy, transform)

      ; plot result

      !P.multi=[0,2,0,0,0]
          
      plot, a.rxy, a.zxy, psym=3, /iso, title="Original grid", color=1
      oplot, a.rxy[edge, *], a.zxy[edge, *], color=2, thick=1.5
      
      plot, Rxy2, Zxy2, psym=3, /iso, title="Poloidal grid", color=1
      oplot, a.rxy[edge, *], a.zxy[edge, *], color=2, thick=1.5
     
      if keyword_set(INPUT) then done=input.grid_ok else $
        done = get_yesno("Is this ok?")
  ENDREP UNTIL done EQ 1

  Rxy = Rxy2
  Zxy = Zxy2

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Output UEDGE style PDB files (optional)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  IF KEYWORD_SET(uedge) THEN BEGIN
      save_uedge, a, min_x, max_x, edge, transform
      RETURN
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Output directly to BOUT++
  ; requires calculating all the metrics, curvature
  ; terms etc. directly
  ; NOTE: Uses shifted coordinates (NO I IN METRIC)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF NOT KEYWORD_SET(output) THEN output="reduced.grd.pdb"

  ; interpolate all values onto new grid points
  PRINT, "Interpolating values onto new grid"
  
  Jpar = apply_transform(Jpar, transform)
  PRINT, "Jpar"
  Bxy  = apply_transform(Bxy,  transform)
  PRINT, "Bxy"
  Bpxy = apply_transform(Bpxy, transform)
  PRINT, "Bpxy"
  Btxy = apply_transform(Btxy, transform)
  PRINT, "Btxy"

  ; no poloidal dependence - don't interpolate
  Te = Te[*,0:ny-1]
  Ti = Ti[*,0:ny-1]
  Ni = Ni[*,0:ny-1]
  psixy = psixy[*, 0:ny-1]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Output to BOUT/BOUT++
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  data = {nx:nx, ny:ny, $
          Rxy:Rxy, Zxy:Zxy, $
          Bpxy:Bpxy, Btxy:Btxy, $
          psixy:psixy, $
          psi_axis:psi_min, psi_bndry:psi_bndry, $
          jpar:jpar, $
          te:te, ti:ti, ni:ni, $
          dpdpsi:dpdpsi, qsafe:qsafe, $
          rmag:rmag, bmag:bmag}
  
  bout_output, data, output=output, fix=index, input=input,_extra=_extra
  
END
