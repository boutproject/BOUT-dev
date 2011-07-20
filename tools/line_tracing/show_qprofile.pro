pro Show_qProfile,path=path,g_file=g_file,xmin=xmin,xmax=xmax,xarr,qarr,color=color,$
     position=position,debug=debug
;
;
;
  COMMON griddata, g

  IF NOT KEYWORD_SET(path) THEN path='./data'
  IF NOT KEYWORD_SET(g_file) THEN g_file='cbm18_8_y064_x516_090309.nc'
  g=file_import(path+'/'+g_file)

if not keyword_set(color) then color=1

n=100
if not keyword_set(xmin) then xmin=MIN(g.psixy)
if not keyword_set(xmax) then xmax=MAX(g.psixy)


xarr=xmin+(xmax-xmin)*findgen(n)/(n-1)
qarr=fltarr(n)

for i=0,n-1 do begin
    qarr[i]= SafetyFactor(xarr[i])
endfor

xmin=min(xarr)
xmax=max(xarr)
qmin=1. ;min(qarr)
qmax=6. ;max(qarr)

;plot, qarr, xarr, xtit='Safety factor q', ytit='x [weber]', chars=2
;plot, qarr, xarr, xtit="Surface average q", ytit="x [weber]",/xst,/yst, chars=2,xr=[qmin,qmax]
;plot, POS=[0.15, 0.075, 0.3, 0.975],qarr, xarr, xtit="q", ytit="x [weber]",/xst,/yst, chars=2,charthick=3,xr=[qmin,qmax], color=color, thick=3
plot,qarr,xarr,ytit="x [weber]", xtit="q", chars=2,yr=[xmin,xmax],xr=[qmin,qmax],/xst,/yst

;;-show some rational q values

qrat=(18+findgen(20))/15.

qrat=(5+findgen(20))/5.

print,qrat

nrat=n_elements(qrat)


for i=0,nrat-1 do begin

    qval=qrat[i]
    xval=INTERPOL(xarr,qarr,qval)

    oplot, [qmin,qmax], [0,0]+xval,lin=2, color=color
    oplot, [0,0]+qval, [xmin,xmax],lin=3, color=color

endfor



if keyword_set(DEBUG) then STOP
end



