;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Converts a set of UEDGE output files into BOUT++ input grids    ;
;                                                                 ;
; NOTE: Doesn't work with X-points yet                            ;
;
; Example use
; -----------
;
; d = read_uedge(path=".")
; 
; OPENW, f, "test.inp" , /get
; PRINTF, f, d.nx
; PRINTF, f, d.ny
; FOR i=0, d.nx-1 DO BEGIN & FOR j=0,d.ny-1 DO BEGIN & PRINTF, f, d.Rxy[i,j] & ENDFOR & ENDFOR
; CLOSE, f
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION ue_convert, u_var

  s = SIZE(u_var, /dim)
  ny = s[0]
  nx = s[1]
  
  ; remove guard cells
  
  var = DBLARR(nx-2, ny-2)

  IF N_ELEMENTS(s) EQ 3 THEN BEGIN
      FOR x=1, nx-2 DO BEGIN
          FOR y=1, ny-2 DO BEGIN
              var[x-1,y-1] = u_var[y,x,0]
          ENDFOR
      ENDFOR
  ENDIF ELSE BEGIN
      FOR x=1, nx-2 DO BEGIN
          FOR y=1, ny-2 DO BEGIN
              var[x-1,y-1] = u_var[y,x]
          ENDFOR
      ENDFOR
  ENDELSE

  RETURN, var
END

FUNCTION read_uedge, path=path, reform=reform, _extra=_extra
  IF NOT KEYWORD_SET(path) THEN path="."

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Read input files
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  file1=EXPAND_PATH(path+'/gridue.nc')
  file2=EXPAND_PATH(path+'/uedgegrd.nc')
  file3=EXPAND_PATH(path+'/uedgeout.nc')

  PRINT, "READING "+file1
  d1 = file_import(file1)
  PRINT, "READING "+file2
  d2 = file_import(file2)
  PRINT, "READING "+file3
  d3 = file_import(file3)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Reformat input from BASIS
  ;
  ; NOTE: Should be possible to detect when
  ; this is needed automatically
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF KEYWORD_SET(reform) THEN BEGIN
      d1 = reformat_struc(d1)
      d2 = reformat_struc(d2)
      d3 = reformat_struc(d3)
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Extract needed information
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  Rxy = ue_convert(d1.rm_com)
  Zxy = ue_convert(d1.zm_com)
  
  s = SIZE(Rxy, /dim)
  nx = s[0]
  ny = s[1]

  Bpxy = ue_convert(d1.bpol_com)
  Btxy = ue_convert(d1.bphi_com)
  psixy = ue_convert(d1.psi_com)

  ;; Create structure with the necessary values
  data = {nx:nx, ny:ny, $
          Rxy:Rxy, Zxy:Zxy, $
          Bpxy:Bpxy, Btxy:Btxy, $
          psixy:psixy, $
          q:d1.qarr[1:(nx-2)]}

  ;; Check which variables are in d3. These are optional

  list = TAG_NAMES(d3)

  IF in_list(list, "JPAR____VALUE")  THEN data = CREATE_STRUCT(data, "Jpar", ue_convert(d3.jpar____value))
  IF in_list(list, "TE____EV_VALUE") THEN data = CREATE_STRUCT(data, "Te",   ue_convert(d3.te____ev_value))
  IF in_list(list, "TI____EV_VALUE") THEN data = CREATE_STRUCT(data, "Ti",   ue_convert(d3.ti____ev_value))
  IF in_list(list, "NI___1__VALUE")  THEN BEGIN
      data = CREATE_STRUCT(data, "Ni",   ue_convert(d3.ni___1__value))

      ;; NOTE: BOUT_OUTPUT expects density in units of 10^20 m^-3
      data.Ni = data.Ni / 1.0e20
  ENDIF

  ;STOP

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Call BOUT++ grid generator
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
  RETURN, data
END
