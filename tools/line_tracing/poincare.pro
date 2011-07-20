PRO poincare, path=path, g_file=g_file, s_file=s_file, tindex=tindex, xmin=xmin, xmax=xmax, $
              nstart=nstart, nmap=nmap, nsubcycle=nsubcycle, debug=debug, quiet=quiet
;
COMMON griddata, g, deltaZtor, Ntor
COMMON BDATA, bd  ;-bd={apar:apar, nx:nx, ny:ny, nz:nz}

IF NOT KEYWORD_SET(path) THEN path='./data/'
IF NOT KEYWORD_SET(g_file) THEN g_file='cbm18_8_y064_x516_090309.nc'
IF KEYWORD_SET(quiet) THEN quiet=2 ELSE quiet=0

IF NOT KEYWORD_SET(s_file) THEN BEGIN
  IF NOT KEYWORD_SET(tindex) THEN tindex=0
;BEGIN LOAD GRID
  time=collect(var='t_array',path=path,quiet=2)
  tstring=STRING(time[tindex])
  g=file_import(path+'/'+g_file)
  Zmax=collect(var='ZMAX',path=path,quiet=2)
  Zmin=collect(var='ZMIN',path=path,quiet=2) 
  deltaZtor=(Zmax-Zmin)*2*!DPI ;in radians
  ;load bfield and set Ntor
  apar=collect(var='Psi',path=path,tind=[tindex,tindex],quiet=2)
  sz=SIZE(apar)
  nx=sz[1]
  ny=sz[2]
  nz=sz[3]
  apar=apar*g.Rmag ; Convert to SI units [m]
  FOR k=0, nz-1 DO apar[*,*,k]=apar[*,*,k]*g.Bxy ; Convert to [Tm]
  ; Convert apar to closed periodic form
  apar_cp=fltarr(nx,ny,nz+1)
  FOR ix=0,nx-1 DO BEGIN
    FOR jy=0,ny-1 DO BEGIN
      apar_cp[ix,jy,*]=[REFORM(apar[ix,jy,*]),apar[ix,jy,0]]
    ENDFOR
  ENDFOR
  Ntor=nz
  bd={apar:apar_cp, nx:nx, ny:ny, nz:nz}
;END GRID LOAD

  s_file=path+'/puncture_plot.'+STRTRIM(STRING(tindex),2)+'.idl.dat'

  IF NOT KEYWORD_SET(nstart) THEN nstart=10
  IF NOT KEYWORD_SET(nmap) THEN nmap=50
  IF NOT KEYWORD_SET(nsubcycle) THEN nsubcycle=1
  IF NOT KEYWORD_SET(xmin) THEN xmin=MIN(g.psixy)
  IF NOT KEYWORD_SET(xmax) THEN xmax=MAX(g.psixy)
  ;call poincare_main 
  Poincare_Main, tindex=tindex, time=tstring, nstart=nstart, nmap=nmap, nsubcycle=nsubcycle, $
                 xmin=xmin, xmax=xmax, savefile=s_file, debug=debug
ENDIF ELSE s_file=path+'/'+s_file

;PLOTTING BEGINS HERE
;restore the data from the savefile (Poincare_Main produces no output)
IF quiet NE 2 THEN BEGIN
  restore,s_file
  plot, [zmin,zmax], [xminplot,xmaxplot], title="time = "+time ,xtit='z [rad]', $
      ytit='x [weber]', /nod, chars=2,/xst,/yst,color=color
  nPts=n_elements(allPts.x)
  for i=0l,nPts-1 do begin
    plots, allPts.z[i], allPts.x[i], col=allPts.c[i], psym=3
  endfor
; Draw a second, logarithmic axis on the right-hand side of the plot.
  AXIS, YAXIS=1, YRANGE=[4.25,4.85], /SAVE, ystyle=1, ytit='R[m]',color=color, chars=2
ENDIF

END
