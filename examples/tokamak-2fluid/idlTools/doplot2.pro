;===============================================;

Get_bpp_data, du=dubpp, t=tbpp, a=armsbpp

loadct,39
SHOW_NPHI, arms=armsbpp, du=dubpp,/fill,/al, it=99


plot, tbpp, armsbpp.NI_XYZT[10,15,*],/yl, xtitle="time, [s]", ytitle="RMS<Ni>"
mark, x=x1, y=y1 & mark, x=x2, y=y2 & gamma=alog(y2/y1)/(x2-x1) & print, "gamma=", gamma
