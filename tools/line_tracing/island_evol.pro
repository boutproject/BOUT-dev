pro Island_Evol, tstart=tstart, tstop=tstop, path=path, delay=delay
  IF NOT KEYWORD_SET(tstart) THEN tstart=50
  IF NOT KEYWORD_SET(tstop) THEN tstop=75
  IF NOT KEYWORD_SET(path) THEN path='./data/'
  IF NOT KEYWORD_SET(delay) THEN delay=0.5

  FOR i=tstart,tstop DO BEGIN
    s_file='puncture_plot.'+STRTRIM(STRING(i),2)+'.idl.dat'
    p_file='p_plot.'+STRTRIM(STRING(i),2)+'.idl.dat'
    ;print,p_file
    result=file_info(path+'puncture_plot.'+STRTRIM(STRING(i),2)+'.idl.dat')
    IF result.exists THEN pressure_plot,path=path,s_file=s_file,p_file=p_file 
    WAIT, delay
  ENDFOR
END 
