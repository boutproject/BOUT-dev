pro safe_colors, first=first,ci=col

;Set the propper device

IF !D.NAME EQ 'X' THEN device, decompose=0

IF !D.TABLE_SIZE GT 11 THEN BEGIN
;Get the current color table vectors
  tvlct,red,green,blue,/GET

;Getting the length of the color vectors

  elem = N_ELEMENTS(red)

;Setting the first color Index to white
  red(0)= 255
  green(0)= 255
  blue(0) = 255

  IF KEYWORD_SET(first) THEN BEGIN
;Setting the first 10 colors to specific indices
    red[1:10] = 255*[0,1,0,0,0,1,1,0.25,0.5,0.75]
    green[1:10] = 255*[0,0,1,0,1,0,1,0.25,0.5,0.75]
    blue[1:10] = 255*[0,0,0,1,1,1,0,0.25,0.5,0.75]
    !p.color=1
    if (!D.N_COLORS LE 256) AND (!D.NAME EQ 'X') THEN BEGIN
        red[1]=255
        blue[1]=255
        green[1]=255
    ENDIF
  ENDIF ELSE BEGIN
;Setting the last 7 color indices to specific colors
    red(elem-7:*) = 255*[0,1,0,0,0,1,1]
    green(elem-7:*) = 255*[0,0,1,0,1,0,1]
    blue(elem-7:*) =  255*[0,0,0,1,1,1,0]
    !p.color=elem-7
    if (!D.N_COLORS LE 256) AND (!D.NAME EQ 'X') THEN BEGIN
        red[elem-7]=255
        blue[elem-7]=255
        green[elem-7]=255
    ENDIF
  ENDELSE
  IF KEYWORD_SET(ci) THEN col=!p.color
;Writing new colortable
  tvlct, red, green, blue
ENDIF ELSE BEGIN
  print,'Not enough colors available!'
ENDELSE

return
end
