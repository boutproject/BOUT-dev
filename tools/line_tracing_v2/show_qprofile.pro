PRO Show_qProfile, path=path, g_file=g_file, xmin=xmin, xmax=xmax,             $
       nvalue=nvalue, color=color, position=position, debug=debug, reduce=reduce

  COMMON griddata, g, deltaZtor, Ntor

  IF NOT KEYWORD_SET(path) THEN path='./data'
  IF NOT KEYWORD_SET(g_file) THEN g_file='cbm18_8_y064_x516_090309.nc'
  g=file_import(path+'/'+g_file)
  IF NOT KEYWORD_SET(nvalue) THEN nvalue=5.0

  IF NOT KEYWORD_SET(color) THEN color=1

  IF NOT KEYWORD_SET(xmin) THEN xmin=MIN(g.psixy)
  IF NOT KEYWORD_SET(xmax) THEN xmax=MAX(g.psixy)

  IF KEYWORD_SET(reduce) THEN BEGIN
    xmax=MAX(g.psixy)-0.05*(MAX(g.psixy)-MIN(g.psixy))
    xmin=MIN(g.psixy)+0.05*(MAX(g.psixy)-MIN(g.psixy))
  ENDIF

  ;array of npts evenly spaced points between xmin and xmax
  npts=100
  xarr=xmin+(xmax-xmin)*FINDGEN(npts)/(npts-1)
  qarr=FLTARR(npts)

  FOR i=0,npts-1 DO BEGIN
    qarr[i]= SafetyFactor(xarr[i])
  ENDFOR

  xmin=MIN(xarr)
  xmax=MAX(xarr)
  qmin=1. ;min(qarr)
  qmax=6. ;max(qarr)

  plot,qarr,xarr,ytit="x [weber]", xtit="q", chars=2,yr=[xmin,xmax], $
       xr=[qmin,qmax],/xst,/yst
;;-show some rational q values
  ;array of 50 m/n rational surfaces
  qrat=(FLOAT(nvalue)+FINDGEN(500))/FLOAT(nvalue)
  nrat=N_ELEMENTS(qrat)
  FOR i=0,nrat-1 DO BEGIN
    qval=qrat[i]
    xval=INTERPOL(xarr,qarr,qval)
    oplot, [qmin,qmax], [0,0]+xval,lin=2, color=color
;    oplot, [0,0]+qval, [xmin,xmax],lin=3, color=color
  ENDFOR

  IF KEYWORD_SET(DEBUG) THEN STOP
END
