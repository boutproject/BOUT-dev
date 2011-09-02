PRO pressure_plot, path=path, g_file, s_file=s_file, p_file=p_file
; Plot data from the savefile
  IF NOT KEYWORD_SET(path) THEN path='./data/'
  IF NOT KEYWORD_SET(g_file) THEN g_file='cbm18_8_y064_x516_090309.nc'
  IF NOT KEYWORD_SET(s_file) THEN s_file='puncture_plot.idl.dat'

  ;restore the poincare file
  restore, path+'/'+s_file

  IF NOT KEYWORD_SET(p_file) THEN BEGIN
    p_file=path+'/p_plot.'+STRTRIM(STRING(tindex),2)+'.idl.dat'
    grid=file_import(path+'/'+g_file)
    x_psi=grid.psixy[*,32]
    z_psi=DINDGEN(64)/(64-1.0D)*(zmax-zmin)+zmin
    p0=collect(var='P0',path=path,quiet=2)
    ;p0=p0(x) is in (x,y) but need it in terms of (x,z)
    ;no need to do anything but only for MZ=65!!!
    p0=transpose(p0)
    p=collect(var='P',path=path,tind=[tindex,tindex],quiet=2)
    p=transpose(reform(p[*,32,*]))
    pt=p0+p
    save, pt, x_psi, z_psi, tindex, time, f=p_file
  ENDIF ELSE restore,path+'/'+p_file ;restore the pressure file

;  loadct,20
  loadct,39
  contour,pt,z_psi,x_psi,/xst,/yst,nlevels=60,/fill,title='time = '+time,      $
          xtit='z [rad]',ytit='x [weber]',chars=2
  nPts=N_ELEMENTS(allPts.x)
;  loadct,0
  tek_color
  FOR i=0l,nPts-1 DO plots, allPts.z[i], allPts.x[i], col=allpts.c[i], psym=3

END
