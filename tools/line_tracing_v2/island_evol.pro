;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;              ISLAND_EVOL
; 
; Draws Poincare plots at subsequent time steps
; 
; Keywords:
;  path       location of data (DEFAULT='./data/')
;  tstart     starting time index 
;  tstop      stopping time index
;  method     what method the Poincare plots were produced with
;  delay      time to wait between plots
;  pressure   if set, draws Poincare plots on top of pressure contours
;  ps         if set, writes output to postscript files
;
;  Written by Joshua Sauppe, University of Wisconsin-Madison
;
PRO Island_Evol,path=path,tstart=tstart,tstop=tstop,method=method,delay=delay, $
                pressure=pressure,ps=ps
  IF NOT KEYWORD_SET(path) THEN path='./data/'
  IF NOT KEYWORD_SET(tstart) THEN tstart=50
  IF NOT KEYWORD_SET(tstop) THEN tstop=75
  IF NOT KEYWORD_SET(delay) THEN delay=0.5
  IF NOT KEYWORD_SET(method) THEN method='tricube.rk4.100'

  loadct,39
  device,decomposed=0
  tek_color
  IF KEYWORD_SET(ps) THEN set_plot,'ps'

  FOR i=tstart,tstop DO BEGIN
    ts=STRTRIM(STRING(i),2)
    s_file='poincare.'+ts+'.'+method+'.idl.dat'
    s_res=file_info(path+'/'+s_file)
    IF s_res.exists THEN BEGIN
      ;check if pressure contours are requested
      IF KEYWORD_SET(pressure) THEN BEGIN
        p_file='p_plot.'+ts+'.idl.dat'
        p_res=file_info(path+'/'+p_file) 
        IF KEYWORD_SET(ps) THEN device,/col,bits=8,file='p_con.'+ts+'.'        $
         +method+'.ps'
        ;if pressure file already exists, no need to reproduce it
        IF p_res.exists THEN BEGIN
          pressure_plot,path=path,s_file=s_file,p_file=p_file 
        ;if pressure file doesn't exist, produce it
        ENDIF ELSE BEGIN
          pressure_plot,path=path,s_file=s_file
        ENDELSE          
        IF KEYWORD_SET(ps) THEN device,/close
      ;otherwise, just make Poincare plots
      ENDIF ELSE BEGIN
        IF KEYWORD_SET(ps) THEN device,/col,bits=8,file='poincare.'+ts+'.'     $
         +method+'.ps'
        poincare,path=path,s_file=s_file
        IF KEYWORD_SET(ps) THEN device,/close
      ENDELSE
      IF NOT KEYWORD_SET(ps) THEN WAIT, delay
    ENDIF
  ENDFOR
  set_plot,'x'
END
