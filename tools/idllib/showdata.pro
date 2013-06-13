;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;              SHOWDATA
; 
; Display animations of 1D,2D and 3D data.
; 
; Arguments:
;   data     2D,3D or 4D array. Final index assumed to be time
; 
; Defaults:
;  2D data   Animate a line plot
;  3D data   Animate a surface plot
;  4D data   Animate a poloidal cross-section (tokamaks only)
;
; Keywords:
;  /contour       For 3D input, show color contour plot
;  yr=[min,max]   Y range
;  /color         Sets color
;  profile=array  Background profile
;                 data is 3D -> profile is 1D (X)
;                 data is 4D -> profile is 2D (X,Y)
;  chars=size     character size
;  az=angle       Rotate surface plots
;  delay=time     Time delay between plots (default 0.2 seconds)
;  /noscale       By default, all plots are on the same scale.
;                 This changes the scale for each plot's range
;
;  uedge=struct   For 4D data, an imported BOUT++ grid file is needed
;  period=integer What fraction of the torus is simulated (e.g. 10)
;  
;  /addsym        For 2D data (1D plots), add symbols to mark data points
;  
;  /bw            Make contour plots greyscale
;  
;  bmp="string"    Output each plot to a numbered BMP file
;
;  output="string" Output each plot to a numbered PS file
;
;  numoff=integer  Add this to the numbering of output files
;
;  mpeg = "string"  Generate a truecolor MPEG movie with each plot as a separate frame
;
; B.Dudson, University of York
;
; Modified:
;
;   2012/03/12  I. Joseph: add mpeg write capability
;

PRO showdata, data, contour=contour, yr=yr, color=color, $
              profile=profile, chars=chars, $
              az=az, delay=delay, _extra=_extra, $
              noscale=noscale, uedge=uedge, period=period, $
              addsym=addsym, bw=bw, bmp=bmp, output=output, $
              numoff=numoff, nobuffer=nobuffer, mpeg=mpeg, $
              loop=loop

  on_error,2  ; If an error occurs, return to caller

  data = REFORM(data)
  IF NOT KEYWORD_SET(chars) THEN chars = 2
  IF NOT KEYWORD_SET(delay) THEN delay = 0.2
  IF NOT KEYWORD_SET(az) THEN az = -45
  IF NOT KEYWORD_SET(numoff) THEN numoff = 0 ELSE PRINT, "Offsetting numbers by "+STRTRIM(STRING(FIX(numoff)),2)
  numoff = FIX(numoff) ; make sure it's an integer

  IF KEYWORD_SET(bmp) OR KEYWORD_SET(mpeg) THEN nobuffer = 1 ; No buffering

  s = SIZE(data, /dimensions)

  ndims = N_ELEMENTS(s)

  wcur = !D.WINDOW ; Get current window ID
  xsize = !D.X_size ; Size of the current window
  ysize = !D.Y_size 

  IF NOT KEYWORD_SET(nobuffer) THEN BEGIN
    ; Create a pixmap window for double-buffering
    WINDOW, xsize=xsize, ysize=ysize, /PIXMAP, /FREE
    pixid = !D.WINDOW  ; Get window ID
    WSET, pixid
  ENDIF

  IF KEYWORD_SET(mpeg) THEN BEGIN
	 mpegfile = mpeg+".mpeg"
	 mpegid = MPEG_OPEN([xsize,ysize],file=mpegfile)
         PRINT, "Writing to: " +mpegfile
  ENDIF

  IF ndims EQ 2 THEN BEGIN
    IF NOT KEYWORD_SET(yr) THEN yr = [MIN(data),MAX(data)]

    nt = s[1]
    FOR t=0, nt-1 DO BEGIN
      IF KEYWORD_SET(noscale) THEN BEGIN
        plot, data[*,t], title="time = "+strtrim(string(t),2), chars=chars, _extra=_extra
        
        IF KEYWORD_SET(addsym) THEN oplot, data[*,t], psym=addsym
      ENDIF ELSE BEGIN
        plot, data[*,t], yr=yr, chars=chars, title="time = "+strtrim(string(t),2), _extra=_extra
        IF KEYWORD_SET(addsym) THEN oplot, data[*,t], psym=addsym
      ENDELSE
      
      IF NOT KEYWORD_SET(nobuffer) THEN BEGIN
        ; Copy buffer
        WSET, wcur
        DEVICE, copy=[0,0,xsize,ysize,0,0,pixid]
        WSET, pixid
      ENDIF
      
      IF delay LT 0.0 THEN BEGIN
        cursor, x, y, /down
      ENDIF ELSE WAIT, delay
    ENDFOR
    
  ENDIF ELSE IF ndims EQ 3 THEN BEGIN
      nx = s[0]
      nz = s[1]
      nt = s[2]

      val = data
      
      IF KEYWORD_SET(profile) THEN BEGIN
          FOR i=0, nx-1 DO BEGIN
              FOR j=0, nz-1 DO BEGIN
                  FOR t=0,nt-1 DO BEGIN
                      val[i,j,t] = val[i,j,t] + profile[i]
                  ENDFOR
              ENDFOR
          ENDFOR
      ENDIF

      IF NOT KEYWORD_SET(yr) THEN yr = [MIN(val),MAX(val)]

      IF KEYWORD_SET(contour) THEN BEGIN
          IF KEYWORD_SET(bw) THEN BEGIN
             LOADCT, 0
          ENDIF ELSE BEGIN
             loadct, 39

          ENDELSE
          device, decomposed=0
                                ;safe_colors, /first      

          REPEAT BEGIN
          FOR i=0, nt-1 DO BEGIN
              contour, reform(val[*,*,i]), chars=chars, zr=yr, zstyle=1, $
                /fill, nlev=50, color=color, $
                title="time = "+strtrim(string(i),2), _extra=_extra
              
              IF NOT KEYWORD_SET(nobuffer) THEN BEGIN
                ; Copy buffer
                WSET, wcur
                DEVICE, copy=[0,0,xsize,ysize,0,0,pixid]
                WSET, pixid
              ENDIF
              
              IF delay LT 0.0 THEN BEGIN
                  cursor, x, y, /down
              ENDIF ELSE WAIT, delay

              IF KEYWORD_SET(bmp) THEN BEGIN
                  IF i EQ 0 THEN BEGIN
                      PRINT, "Click window to begin"
                      cursor, x, y, /down
                  ENDIF
                  file = bmp + STRTRIM(STRING(i+numoff, FORMAT='(I04)'),2)+".bmp"
                  PRINT, "Writing file: " +file
                  
                  WRITE_BMP, file, TVRD(TRUE=1)
              ENDIF
              IF KEYWORD_SET(mpeg) THEN BEGIN
                  PRINT, "  frame "+strtrim(string(i+1),1)+"/"+strtrim(string(nt),1)
	          frame = tvrd(true=1,/order)
                  MPEG_PUT, mpegid, image=frame, frame=i
              ENDIF

            ENDFOR
        ENDREP UNTIL NOT KEYWORD_SET(loop)
      ENDIF ELSE BEGIN
          PRINT, "chars=", chars
          FOR i=0, nt-1 DO BEGIN
              surface, reform(val[*,*,i]), chars=chars, zr=yr, zstyle=1, $
                color=color, az=az, title="time = "+strtrim(string(i),2), _extra=_extra
              IF NOT KEYWORD_SET(nobuffer) THEN BEGIN
                ; Copy buffer
                WSET, wcur
                DEVICE, copy=[0,0,xsize,ysize,0,0,pixid]
                WSET, pixid
              ENDIF
              
              IF delay LT 0.0 THEN BEGIN
                  cursor, x, y, /down
              ENDIF ELSE WAIT, delay
          ENDFOR
      ENDELSE
      
  ENDIF ELSE IF ndims EQ 4 THEN BEGIN
      IF NOT KEYWORD_SET(uedge) THEN BEGIN
          PRINT, "For 4D, need to give uedge data"
          RETURN
      ENDIF

      nx = s[0]
      ny = s[1]
      nz = s[2]
      nt = s[3]

      loadct, 39
      device, decomposed=0
       
      IF NOT KEYWORD_SET(noscale) THEN yr=[MIN(data), MAX(data)]
          
      FOR i=0, nt-1 DO BEGIN
          ;d2d = get_pol_slice(reform(data[*,*,*,i]), uedge, n=period)
          ;s = SIZE(d2d, /dimensions)
          ;nx = s[0]
          ;ny = s[1]
          ;IF KEYWORD_SET(profile) THEN BEGIN
          ;    FOR x=0, nx-1 DO BEGIN
          ;        FOR y=0, ny-1 DO BEGIN
          ;            d2d[x,y] = d2d[x,y] + profile[x,y]
          ;        ENDFOR
          ;    ENDFOR
          ;ENDIF
          ;plot2d, d2d, uedge, yr=yr, title="Time = "+STRTRIM(STRING(i),2)

          IF KEYWORD_SET(output) THEN  BEGIN
              psfile = output + STRTRIM(STRING(i, FORMAT='(I04)'),2)+".ps"
              PRINT, "Setting psfile = "+psfile
          ENDIF

          plotpolslice, REFORM(data[*,*,*,i]), uedge, $
            period=period, profile=profile, output=psfile, _extra=_extra
          
          IF NOT KEYWORD_SET(nobuffer) THEN BEGIN
            ; Copy buffer
            WSET, wcur
            DEVICE, copy=[0,0,xsize,ysize,0,0,pixid]
            WSET, pixid
          ENDIF
          
          IF delay LT 0.0 THEN BEGIN
              cursor, x, y, /down
          ENDIF ELSE WAIT, delay

          IF KEYWORD_SET(bmp) THEN BEGIN
              IF i EQ 0 THEN BEGIN
                  PRINT, "Click window to begin"
                  cursor, x, y, /down
              ENDIF
              file = bmp + STRTRIM(STRING(i+numoff, FORMAT='(I04)'),2)+".bmp"
              PRINT, "Writing file: " +file
              
              WRITE_BMP, file, TVRD(TRUE=1)
          ENDIF

          IF KEYWORD_SET(mpeg) THEN BEGIN
;              IF i EQ 0 THEN BEGIN
;                  PRINT, "Click window to begin"
;                  cursor, x, y, /down
;              ENDIF		
              PRINT, "  frame "+strtrim(string(i+1),1)+"/"+strtrim(string(nt),1)
              frame = tvrd(true=1,/order)
              MPEG_PUT, mpegid, image=frame, frame=i
          ENDIF

      ENDFOR
  ENDIF ELSE BEGIN
      PRINT, "Data must be either 2, 3 or 4 dimensional"
  ENDELSE

  IF KEYWORD_SET(mpeg) THEN BEGIN
	PRINT, "Saving MPEG file: "+mpegfile
	MPEG_SAVE,mpegid
	MPEG_CLOSE,mpegid
  ENDIF

  IF NOT KEYWORD_SET(nobuffer) THEN BEGIN
    ; Switch back to original window
    WSET, wcur
    
    ; Delete the pixmap window
    WDELETE, pixid
  ENDIF

END
