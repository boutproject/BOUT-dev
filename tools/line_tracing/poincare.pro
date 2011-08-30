;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;              POINCARE
; 
; Produce Poincare plot on poloidal slice at y=pi (outer midplane)
; 
; Keywords:
;  path       location of data (DEFAULT='./data/')
;  g_file     grid file for BOUT++ run (DEFAULT='cbm18_8_y064_x516_090309.nc')
;  s_file     save file for Poincare plots, if set data is loaded and plotted
;  tindex     time index to load data from (DEFAULT=0)
;  xmin       minimum Psi for plots (DEFAULT=MIN(g.psixy)+0.05*delta_x
;  xmax       maximum Psi for plots (DEFAULT=MAX(g.psixy)-0.05*delta_x
;  nstart     number of initial conditions in x (DEFAULT=10)
;  nmap       number of toroidal transits for field line (DEFAULT=50)
;  nsubcycle  number of initial conditions in z (DEFAULT=1) 
;  /debug     if set, enables STOP statements throughout for debugging
;  /quiet     if set, collects data without printing to screen
;  /model     if set, use the function in apar_model.pro on grid points
;             the default behavior is to load data from the BOUT dump files
;  method     the interpolation scheme to use (DEFAULT='trilin')
;             'trilin'    - trilinear interpolation in x,y,z
;             'bilin_fft' - bilinear interpolation in x,y and Fourier in z
;             'bilin_cube'- bilinear interpolation in x,y and cubic in z 
;             'tricube'   - tricubic interpolation in x,y,z
;             'bicube_fft'- bicubic interpolation in x,y and Fourier in z
;  mc_method  interpolation scheme to use for metric coeffs (DEFAULT=bilinear)
;             'close'     - original scheme, just uses lower-left index value
;             'bilinear'  - bilinear interpolation of metric coefficients
;             'bicubic'   - bicubic interpolation of metric coefficients
;
;  Modifications by Joshua Sauppe, University of Wisconsin-Madison
;  Based on code initially written by Maxim Umansky, LLNL
;
PRO poincare, path=path, g_file=g_file, s_file=s_file, xmin=xmin, xmax=xmax,   $
              tindex=tindex, nstart=nstart, nmap=nmap, nsubcycle=nsubcycle,    $
              model=model, method=method, mc_method=mc_method, debug=debug,    $
              quiet=quiet
;
COMMON griddata, g, deltaZtor, Ntor
COMMON BDATA, bd                  ; bd={apar:apar, nx:nx, ny:ny, nz:nz}
COMMON flags, flag, mc_flag       ; stores keywords for interpolation schemes
COMMON metric_data, mc_data       ; stores coefficients of grid interpolation
COMMON bicube_data, ad            ; ad={TF:TF,coeffs:FLTARR(16) interp coeffs}

IF NOT KEYWORD_SET(path) THEN path='./data/'
IF NOT KEYWORD_SET(g_file) THEN g_file='cbm18_8_y064_x516_090309.nc'
IF KEYWORD_SET(quiet) THEN quiet=2 ELSE quiet=0

IF NOT KEYWORD_SET(s_file) THEN BEGIN
  IF NOT KEYWORD_SET(tindex) THEN tindex=0
;;select the interpolation method to use
  IF NOT KEYWORD_SET(method) THEN method='trilin'
  flag=method
  IF NOT KEYWORD_SET(mc_method) THEN mc_method='bilinear'
  mc_flag=mc_method

;;collect some basic info about the grid
  time=collect(var='t_array',path=path,quiet=2)
  tstring=STRING(time[tindex])
  g=file_import(path+'/'+g_file)
  Zmax=collect(var='ZMAX',path=path,quiet=2)
  Zmin=collect(var='ZMIN',path=path,quiet=2) 
  deltaZtor=(Zmax-Zmin)*2*!DPI ;length of z direction in radians

;;load either the perturbation data, or use model function evaluated on grid
  IF NOT KEYWORD_SET(model) THEN BEGIN
  ;load bfield and set Ntor
    IF quiet NE 2 THEN PRINT,'Loading the perturbation data...'
    apar=collect(var='Psi',path=path,tind=[tindex,tindex],quiet=2)
    sz=SIZE(apar)
    nx=sz[1]
    ny=sz[2]
    nz=sz[3]
    Ntor=nz
    ;Apar=Psi*B so multiply by B_0
    apar=apar*g.Rmag ; Convert to SI units [m]
    FOR k=0, nz-1 DO apar[*,*,k]=apar[*,*,k]*g.Bxy ; Convert to [Tm]
  ENDIF ELSE BEGIN
  ;use the model apar evaluated at the grid points
    IF quiet NE 2 THEN PRINT,'Using the model function...'
    nx=g.nx
    ny=g.ny
    nz=64         ;use 64 toroidal grid points
    Ntor=nz
    apar=FLTARR(g.nx,g.ny,nz)
    FOR ix=0,g.nx-1 DO BEGIN
      FOR iy=0,g.ny-1 DO BEGIN
        FOR iz=0,nz-1 DO BEGIN
          dApar=apar_model(g.psixy[ix,iy],iy*2*!DPI/FLOAT(g.ny),               $
                           iz*deltaZtor/FLOAT(nz))
          apar[ix,iy,iz]=dApar.f
        ENDFOR
      ENDFOR
    ENDFOR
  ENDELSE

;;Transform the perturbation into closed periodic form
  CASE flag OF
    'trilin'      :   cp_trilin,apar=apar
    'bilin_cube'  :   cp_bilin_cube,apar=apar
    'bilin_fft'   :   cp_bilin_fft,apar=apar
    'tricube'     :   cp_tricube,apar=apar
    'bicube_fft'  :   BEGIN
                        cp_bicube_fft,apar=apar
                        ad_pt={TF:0,coeffs:FLTARR(16)}
                        ;NOTE the coeff array indices differ from apar indices!
                        ad=REPLICATE(ad_pt,g.nx,g.ny,nz)
                      END
    'model'       :   cp_trilin,apar=apar
     else         :   BEGIN
                        cp_trilin,apar=apar
                        flag='trilin'
                      END
  ENDCASE

;;Initialize the storage arrays for metric coefficients
  ;;   mc_data[ix_int,iy_int] = 
  ;;             {    TF:TF,       True/False flag if computed
  ;;                 rxy:rxy,    \
  ;;                 bxy:bxy,     \    FLTARR(coeffs)
  ;;                bpxy:bpxy,     \   coeffs= 0 if close
  ;;                btxy:btxy,      >  coeffs= 4 if bilinear
  ;;               sinty:sinty,    /   coeffs=16 if bicubic
  ;;                hthe:hthe,    / 
  ;;               bxcvx:bxcvx,  / 
  ;;               bxcvz:bxcvz    }
  CASE mc_flag OF
    'close'    :  mc_pt={TF:0}
    'bilinear' :  mc_pt={TF:0,rxy:FLTARR(4),bxy:FLTARR(4),bpxy:FLTARR(4),      $
                         btxy:FLTARR(4),sinty:FLTARR(4),hthe:FLTARR(4),        $
                         bxcvx:FLTARR(4),bxcvz:FLTARR(4)                 }
    'bicubic'  :  mc_pt={TF:0,rxy:FLTARR(16),bxy:FLTARR(16),bpxy:FLTARR(16),   $
                         btxy:FLTARR(16),sinty:FLTARR(16),hthe:FLTARR(16),     $
                         bxcvx:FLTARR(16),bxcvz:FLTARR(16)                 }
     else      :  BEGIN
                    mc_pt={TF:0,rxy:FLTARR(4),bxy:FLTARR(4),bpxy:FLTARR(4),    $
                           btxy:FLTARR(4),sinty:FLTARR(4),hthe:FLTARR(4),      $
                           bxcvx:FLTARR(4),bxcvz:FLTARR(4)                 }
                    mc_flag='bilinear'
                  END
  ENDCASE
  mc_data=REPLICATE(mc_pt,g.nx,g.ny)

  bd={apar:apar, nx:nx, ny:ny, nz:nz}
;END GRID LOAD

  s_file=path+'/poincare.'+STRTRIM(STRING(tindex),2)+'.'+flag+'.idl.dat'

  IF NOT KEYWORD_SET(nstart) THEN nstart=10
  IF NOT KEYWORD_SET(nmap) THEN nmap=50
  IF NOT KEYWORD_SET(nsubcycle) THEN nsubcycle=1
  ;if starting points are too close to boundary, field lines may exit domain
  delta_x=MAX(g.psixy)-MIN(g.psixy)
  IF NOT KEYWORD_SET(xmin) THEN xmin=MIN(g.psixy)+0.05*delta_x
  IF NOT KEYWORD_SET(xmax) THEN xmax=MAX(g.psixy)-0.05*delta_x
  ;call poincare_main 
  Poincare_Main, tindex=tindex, time=tstring, nstart=nstart, nmap=nmap,        $
                 nsubcycle=nsubcycle, xmin=xmin, xmax=xmax, savefile=s_file,   $
                 debug=debug
ENDIF ELSE s_file=path+'/'+s_file

;PLOTTING BEGINS HERE
;restore the data from the savefile (Poincare_Main produces no output)
IF quiet NE 2 THEN BEGIN
  RESTORE,s_file
  PLOT, [zmin,zmax], [xminplot,xmaxplot], title="time = "+time ,xtit='z [rad]',$
      ytit='x [weber]', /nod, chars=2,/xst,/yst,color=color
  nPts=N_ELEMENTS(allPts.x)
  FOR i=0l,nPts-1 DO BEGIN
;    PLOTS, allPts.z[i], allPts.x[i], col=allPts.c[i], psym=3
    PLOTS, allPts.z[i],allPts.x[i],col=allPts.c[i],psym=6,symsize=0.05
  ENDFOR
; Draw a second, logarithmic axis on the right-hand side of the plot.
;  AXIS, YAXIS=1, YRANGE=[4.25,4.85], /SAVE, ystyle=1, ytit='R[m]',color=color
ENDIF

END
