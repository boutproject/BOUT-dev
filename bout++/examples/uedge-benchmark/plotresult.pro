; Produce plots of the result
; Keywords:
;    /orig       Over-plot original result as symbols
;    /output     Write plots to a PS file result.ps
;    /nodata     Reproduce original plot (no new data)
;    /nouedge    Don't overplot UEDGE result

PRO plotresult, orig=orig, output=output, nodata=nodata, nouedge=nouedge

path = "data/"

; set the time axis range and units
tmax = 4.0
tfactor = 1000.
xtitle="Time [ms]"

IF NOT KEYWORD_SET(nodata) THEN BEGIN

;; read time-base

tarr = collect(path=path, var="t_array")
wci = collect(path=path, var="wci")

time = tarr*tfactor / wci

nt = N_ELEMENTS(tarr)

;; collect data

ni = collect(path=path, var="Ni", y=40)
vi = collect(path=path, var="Vi", y=40)
ti = collect(path=path, var="Ti", y=40)
te = collect(path=path, var="Te", y=40)

;; read normalisations

ni_x = collect(path=path, var="Ni_x") 
rho_s = collect(path=path, var="rho_s") 
Vi_x = rho_s * wci
Te_x = collect(path=path, var="Te_x") 

;; read grid file

u = file_import("uedge.grd_Up_Ni_Tei_2d.nc")

;; Normalise

nit = reform(ni[*,0,0,*])*ni_x
FOR t=0, nt-1 DO nit[*,t] = nit[*,t] + u.ni0[*,40]*1.e14

vit = REFORM(vi) * vi_x

tit = reform(ti[*,0,0,*])*te_x
FOR t=0, nt-1 DO tit[*,t] = tit[*,t] + u.Ti0[*,40]

tet = reform(te[*,0,0,*])*te_x
FOR t=0, nt-1 DO tet[*,t] = tet[*,t] + u.Te0[*,40]

ENDIF ELSE BEGIN
    ; Make original plot (no new data)

    RESTORE, "result_080917.idl"

    time = time_orig
    nit = ni_orig
    vit = vi_orig
    tit = ti_orig
    tet = te_orig
    orig = 0
ENDELSE

;; restore the original data (if needed)

IF KEYWORD_SET(orig) THEN BEGIN
    RESTORE, "result_080917.idl"
ENDIF

IF NOT KEYWORD_SET(nouedge) THEN BEGIN
    RESTORE, "ue_bmk.idl"

    time_ue = time_ue * tfactor

    ni_ue = ni_ue * 1.e-6 ; convert to cm^-3
ENDIF

;; Produce plots

IF KEYWORD_SET(output) THEN BEGIN
    SET_PLOT, 'PS'
    DEVICE, file="result.ps", /color, /landscape
ENDIF

safe_colors, /first

!P.MULTI=[0,2,2,0,0]

;; NI

plot, time, nit[40,*], color=1, xr=[0,tmax], xstyle=1, xtitle=xtitle, title="Ni"
oplot, time, nit[20,*], color=2
oplot, time, nit[10,*], color=4

IF KEYWORD_SET(orig) THEN BEGIN
    oplot, time_orig, ni_orig[40,*], color=1, lines=1
    oplot, time_orig, ni_orig[20,*], color=2, lines=1
    oplot, time_orig, ni_orig[10,*], color=4, lines=1
ENDIF

IF NOT KEYWORD_SET(nouedge) THEN BEGIN
    oplot, time_ue, ni_ue[*,10], color=3, thick=1.5
    oplot, time_ue, ni_ue[*,20], color=3, thick=1.5
    oplot, time_ue, ni_ue[*,40], color=3, thick=1.5
ENDIF

;; VI

plot, time, vit[40,*], color=1, xr=[0,tmax], xstyle=1, xtitle=xtitle, title="Vi"
oplot, time, vit[20,*], color=2
oplot, time, vit[10,*], color=4

IF KEYWORD_SET(orig) THEN BEGIN
    oplot, time_orig, vi_orig[40,*], color=1, lines=1
    oplot, time_orig, vi_orig[20,*], color=2, lines=1
    oplot, time_orig, vi_orig[10,*], color=4, lines=1
ENDIF

IF NOT KEYWORD_SET(nouedge) THEN BEGIN
    oplot, time_ue, vi_ue[*,10], color=3, thick=1.5
    oplot, time_ue, vi_ue[*,20], color=3, thick=1.5
    oplot, time_ue, vi_ue[*,40], color=3, thick=1.5
ENDIF

;; TI

plot, time, tit[40,*], color=1, xr=[0,tmax], xstyle=1, xtitle=xtitle, title="Ti", yr=[0,MAX(tit)]
oplot, time, tit[20,*], color=2
oplot, time, tit[10,*], color=4

IF KEYWORD_SET(orig) THEN BEGIN
    oplot, time_orig, ti_orig[40,*], color=1, lines=1
    oplot, time_orig, ti_orig[20,*], color=2, lines=1
    oplot, time_orig, ti_orig[10,*], color=4, lines=1
ENDIF

IF NOT KEYWORD_SET(nouedge) THEN BEGIN
    oplot, time_ue, ti_ue[*,10], color=3, thick=1.5
    oplot, time_ue, ti_ue[*,20], color=3, thick=1.5
    oplot, time_ue, ti_ue[*,40], color=3, thick=1.5
ENDIF

;; TE

plot, time, tet[40,*], color=1, xr=[0,tmax], xstyle=1, xtitle=xtitle, title="Te", yr=[0,MAX(tet)]
oplot, time, tet[20,*], color=2
oplot, time, tet[10,*], color=4

IF KEYWORD_SET(orig) THEN BEGIN
    oplot, time_orig, te_orig[40,*], color=1, lines=1
    oplot, time_orig, te_orig[20,*], color=2, lines=1
    oplot, time_orig, te_orig[10,*], color=4, lines=1
ENDIF

IF NOT KEYWORD_SET(nouedge) THEN BEGIN
    oplot, time_ue, te_ue[*,10], color=3, thick=1.5
    oplot, time_ue, te_ue[*,20], color=3, thick=1.5
    oplot, time_ue, te_ue[*,40], color=3, thick=1.5
ENDIF

IF KEYWORD_SET(output) THEN BEGIN
    DEVICE, /close
    SET_PLOT, 'X'
ENDIF

!P.MULTI=0

STOP
END

