.run pdb2idl.pro
.run moment_xyzt.pro

d = collect(path="data_1", var="Ni", x=5, y=30)
moment_xyzt, d, rms=ni1

d = collect(path="data_10", var="Ni", x=5, y=30)
moment_xyzt, d, rms=ni10
tt = collect(path="data_10", var="t_array")
wci = collect(path="data_10", var="wci")
tt = tt[1:*] / wci

ni1 = REFORM(ni1[0,0,1:*])
ni10 = REFORM(ni10[0,0,1:*])

;-compare with original test results (grendel, 31-jan-2007)
 RESTORE, 'orig_test.idl.dat'
 error1=max(abs((ni1orig-ni1)/ni1orig)) + max(abs((ni10orig-ni10)/ni10orig))
 print, "Deviation from original test result is" , error1*1e2, " %"

set_plot, 'PS'
device, file='interchange_inst_test.ps'
safe_colors, /first

xtit='t, s' & ytit='rms <Ni>' & tit='Interchange instability test'
plot, tt, ni1, /yl, psym=4, xtit=xtit, ytit=ytit, tit=tit, chars=1.5,col=1
 oplot, tt, 1e-4*exp(tt*2.2e5),col=1
 oplot, tt, ni10, psym=4,col=1
 oplot, tt, 1e-4*exp(tt*6.3e4),col=1, lin=2

 oplot, tt, ni1orig,psym=7,col=1
 oplot, tt, ni10orig, psym=7,col=1

 xyouts, 7e-5, 0.15, "R0=10 m",col=1
 xyouts, 4e-5, 1e2, "R0=1 m",col=1
 
; calculate growth rates

n1g0 = MEAN(DERIV(tt, ALOG(ni1orig)))
n1g  = MEAN(DERIV(tt, ALOG(ni1)))

n1d = 100. * ABS(n1g0 - n1g) / n1g0

n10g0 = MEAN(DERIV(tt, ALOG(ni10orig)))
n10g  = MEAN(DERIV(tt, ALOG(ni10)))

n10d  = 100. * ABS(n10g0 - n10g) / n10g0

PRINT, ""
PRINT, "Growth rates:"
PRINT, "============="
PRINT, "R0 = 1m Growth = " + STRTRIM(STRING(n1g),2) + $
  " Analytic = 2.2e5 Orig = "+STRTRIM(STRING(n1g0),2)
PRINT, "    Deviation from original = " + STRTRIM(STRING(n1d),2)+ "%"
PRINT, ""
PRINT, "R0 = 10m Growth = " + STRTRIM(STRING(n10g),2) + $
  " Analytic = 6.3e4 Orig = "+STRTRIM(STRING(n10g0),2)
PRINT, "    Deviation from original = " + STRTRIM(STRING(n10d),2)+ "%"
PRINT, ""
PRINT, "DONE. Graph written to interchange_inst_test.ps"
PRINT, ""
device, /close
set_plot, 'X'

exit, status=status
