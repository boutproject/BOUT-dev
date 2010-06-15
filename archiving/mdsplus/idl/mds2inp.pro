;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Reads settings from an MDSplus shot, outputting a BOUT.inp file
; Currently just prints the settings to stdout
;
; KEYWORDS
;
; server     Name (and optionally the port) of MDS server
;            Default 'localhost'
; tree       MDS tree. Default 'bout_mds'
; shot       Shot number. By default reads the latest
;
; Run server: mdsip -p 8000 -m -h $MDS_ROOT/etc/mdsip.hosts
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO print_settings, path, nremove, prefix, fd
  
  setting_nid=MDSVALUE('GETNCI("'+path+'","NID_NUMBER")',status=status)
  IF NOT (status AND 1) THEN RETURN

  n = N_ELEMENTS(setting_nid)
  ;PRINT, "Number of settings:", n

  nvars = 0  ; keep track of the number of variables
  FOR i=0, n-1 DO BEGIN
    str = mdsvalue('GETNCI($,"MINPATH")',setting_nid[i])
    ;PRINT, str
    
    ; look for members
    pos = STRPOS(str, ':')

    IF pos NE -1 THEN BEGIN
      ; this string is a leaf
      
      len = STRLEN(str)

      name = STRMID(str, pos+1, len - pos-1)
      path = STRMID(str, nremove, pos - nremove)
      
      ; replace '.' with '/'
      STRREPLACE, path, '.', '/'
      
      ;PRINT, "-> '" + path + "' - '"+name+"'"
      
      IF nvars EQ 0 THEN BEGIN
        paths = [path]
        names = [name]
        nids = [setting_nid[i]]
      ENDIF ELSE BEGIN
        paths = [paths, path]
        names = [names, name]
        nids = [nids, setting_nid[i]]
      ENDELSE
      nvars = nvars + 1
      
    ENDIF
    
  ENDFOR
  
  ; get the settings without section

  IF STRLEN(prefix) GT 0 THEN PRINTF, fd, '['+prefix+']'

  w = WHERE(paths EQ '', count)
  IF count GT 0 THEN BEGIN
    FOR i=0, count-1 DO BEGIN
      ; get the value
      p = MDSVALUE('GETNCI($, "PATH")', nids[i])
      
      v = MDSVALUE(p, status=status, /quiet)
      
      IF status AND 1 THEN BEGIN
        ; success
        PRINTF, fd, names[i]+" = " + STRING(v)
      ENDIF ELSE BEGIN
        ; failed - no data
        
        PRINTF, fd, '# '+names[i]+ ' = '
      ENDELSE
    ENDFOR

    w = WHERE(paths NE '', count)
    IF count GT 0 THEN BEGIN
      paths = paths[w]
      names = names[w]
      nids = nids[w]
    ENDIF
    nvars = count
  ENDIF

  WHILE nvars GT 0 DO BEGIN
    IF STRLEN(prefix) GT 0 THEN PRINTF, fd, '['+prefix+'/'+paths[0]+']' ELSE PRINTF, fd, '['+paths[0]+']'
    
    w = WHERE(paths EQ paths[0], count)

    FOR i=0, count-1 DO BEGIN
      ; get the value
      p = MDSVALUE('GETNCI($, "PATH")', nids[w[i]])
      v = MDSVALUE(p, status=status, /quiet)
      
      IF status AND 1 THEN BEGIN
        ; success
        PRINTF, fd, names[i]+" = " + STRING(v)
      ENDIF ELSE BEGIN
        ; failed - no data
        
        PRINTF, fd, '# '+names[i]+ ' = '
      ENDELSE
    ENDFOR

    w = WHERE(paths NE paths[0], count)
    IF count GT 0 THEN BEGIN
      paths = paths[w]
      names = names[w]
      nids = nids[w]
    ENDIF
    nvars = count
  ENDWHILE
END

PRO mds2inp, server=server, tree=tree, shot=shot
  
  IF NOT KEYWORD_SET(server) THEN server = "localhost"
  IF NOT KEYWORD_SET(tree) THEN tree = "bout_mds"

  mdsconnect, server, status = status

  IF NOT (status AND 1) THEN BEGIN
    PRINT, "Could not connect"
    RETURN
  ENDIF

  IF NOT KEYWORD_SET(shot) THEN BEGIN
    ; get the latest shot
    shot = MDSVALUE('current_shot("'+tree+'")')
    PRINT, "Latest shot = ", shot
  ENDIF

  mdsopen, tree, shot, status=status
  IF NOT (status AND 1) THEN BEGIN
    PRINT, "Could not open tree"
    mdsdisconnect
    RETURN
  ENDIF

  ; get the description and model

  desc = MDSVALUE(':DESC', status=status)
  IF NOT (status AND 1) THEN desc = 'No description'
  model = MDSVALUE(':MODEL', status=status)
  IF NOT (status AND 1) THEN model = 'NO MODEL'
  
  ; NOTE: This method could be replaced with calls for CHILDREN_NIDS
  ; but that doesn't seem to work (says no data for node)

  fd = -1
  
  PRINTF, fd, "# BOUT++ input file. Generated from MDSplus data"
  PRINTF, fd, "# Date: " + SYSTIME()
  PRINTF, fd, "# Shot: " + STRTRIM(STRING(shot),2)
  PRINTF, fd, "# Description: " + desc
  PRINTF, fd, "# Model: " + model
  PRINTF, fd, ""
  
  settings_path = "\\TOP.SETTINGS***"
  print_settings, settings_path, 10, '', fd

  print_settings, "\\TOP."+model+".SETTINGS***", 10+STRLEN(model), model, fd
  

  mdsclose
  mdsdisconnect
  
END


