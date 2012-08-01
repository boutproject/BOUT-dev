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
;  method     the interpolation scheme to use (DEFAULT='tricube')
;             'tricube'    - tricubic interpolation in x,y,z
;             'bicube_fft' - bicubic interpolation in x,y Fourier in z
;             'model'      - use analytic derivatives 
;
;  Modifications by Joshua Sauppe, University of Wisconsin-Madison
;  Based on code initially written by Maxim Umansky, LLNL
;
PRO poincare, path=path, g_file=g_file, s_file=s_file, xmin=xmin, xmax=xmax,   $
              tindex=tindex, nstart=nstart, nmap=nmap, nsubcycle=nsubcycle,    $
              model=model, method=method, int_method=int_method,               $
              int_steps=int_steps, quiet=quiet, debug=debug

  COMMON griddata, g, deltaZtor, Ntor
  COMMON BDATA, bd                   ; bd={apar:apar, nx:nx, ny:ny, nz:nz}
  COMMON flags, flag                 ; stores interpolation method to use
  COMMON metric_data, mc_data        ; stores grid data bicubic coefficients 
  COMMON bicube_data, ad             ; ad={TF:(0 or 1),coeffs:FLTARR(16)}
  COMMON tricube_data, A_int, A_inv  ; A_int={TF:(0 or 1),coeffs:FLTARR(64)}
  COMMON int_data, int_meth, nsteps  ; data for numerical field-line integrator

  IF NOT KEYWORD_SET(method) THEN method='tricube'
  CASE method OF
    'tricube'    : flag='tricube'
    'bicube_fft' : flag='bicube_fft'
    'model'      : flag='model'
     else        : flag='tricube'
  ENDCASE

  IF NOT KEYWORD_SET(int_method) THEN int_method='rk4'
  IF NOT KEYWORD_SET(int_steps) THEN int_steps=100
  CASE int_method OF
    'rk4'   :  int_meth='rk4'
    'cn'    :  int_meth='cn'
     else   :  int_meth='rk4'
  ENDCASE
  nsteps=int_steps


  IF NOT KEYWORD_SET(path) THEN path='./data/'
  IF NOT KEYWORD_SET(g_file) THEN g_file='cbm18_8_y064_x516_090309.nc'
  IF KEYWORD_SET(quiet) THEN quiet=2 ELSE quiet=0


  IF NOT KEYWORD_SET(s_file) THEN BEGIN
;-------------------BEGIN DATA LOAD---------------------------------------------
    IF NOT KEYWORD_SET(tindex) THEN tindex=0
    ;;collect some basic info
    time=collect(var='t_array',path=path,quiet=2)
    tstring=STRING(time[tindex])
    g=file_import(path+'/'+g_file)
    Zmax=collect(var='ZMAX',path=path,quiet=2)
    Zmin=collect(var='ZMIN',path=path,quiet=2) 
    deltaZtor=(Zmax-Zmin)*2*!DPI ;length of z direction in radians

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
      nz=128        ;use 64 toroidal grid points by default
      Ntor=nz
      apar=FLTARR(g.nx,g.ny,nz)
      FOR ix=0,g.nx-1 DO BEGIN
        FOR iy=0,g.ny-1 DO BEGIN
          FOR iz=0,nz-1 DO BEGIN
            dApar=apar_model(g.psixy[ix,iy],iy*2*!DPI/FLOAT(g.ny),             $
                             iz*deltaZtor/FLOAT(nz))
            apar[ix,iy,iz]=dApar.f
          ENDFOR
        ENDFOR
      ENDFOR
    ENDELSE
;---------------------END DATA LOAD---------------------------------------------
    ;;Initialize the storage arrays for metric coefficients
    ;;   mc_data[ix_int,iy_int] = 
    ;;             {    TF:TF,       True/False flag if computed
    ;;                 rxy:rxy,    \
    ;;                 bxy:bxy,     \
    ;;                bpxy:bpxy,     \
    ;;                btxy:btxy,      \  FLTARR(coeffs)
    ;;               sinty:sinty,      > coeffs=0 if close
    ;;                hthe:hthe,      /  coeffs=4 if bilinear 
    ;;               bxcvx:bxcvx,    /   coeffs=16 if bicubic 
    ;;               bxcvx:bxcvx,   /
    ;;               bxcvz:bxcvz,  /
    ;;               jpar0:jpar0  /  }
    mc_pt={TF:0,rxy:FLTARR(16),bxy:FLTARR(16),bpxy:FLTARR(16),btxy:FLTARR(16), $
           sinty:FLTARR(16),hthe:FLTARR(16),bxcvx:FLTARR(16),bxcvy:FLTARR(16), $
           bxcvz:FLTARR(16),jpar0:FLTARR(16)  }
    mc_data=REPLICATE(mc_pt,g.nx,g.ny)


    ;;Initialize storage arrays for (bi/tri)cubic interpolation coefficients
    ;;   A_int[ix_int,iy_int] = 
    ;;             {    TF:TF,       True/False flag if computed
    ;;              coeffs:coeffs    64-element vector of coefficients  }
    CASE flag OF
      'tricube'    : BEGIN
                       A_int_pt={TF:0,coeffs:FLTARR(64)}
                       A_int=REPLICATE(A_int_pt,nx,ny,nz)
                       ;;Calculate the inversion matrix A_inv 
                       A_inv=make_tricube_mat()
                     END
      'bicube_fft' : BEGIN
                       cp_bicube_fft,apar=apar
                       ad_pt={TF:0,coeffs:FLTARR(16)}
                       ;;Note the coeff array indices differ from apar indices!
                       ad=REPLICATE(ad_pt,nx,ny,nz)
                     END
      'model'      : PRINT,'Using the model derivatives'
       else        : STOP,'Error: interpolation flag not specified correctly!'
    ENDCASE

    bd={apar:apar, nx:nx, ny:ny, nz:nz}

;-------------------POINCARE MAIN CALL------------------------------------------
    s_file=path+'/poincare.'+STRTRIM(STRING(tindex),2)+'.'+flag+'.'+int_meth   $
              +'.'+STRTRIM(STRING(nsteps),2)+'.idl.dat'
    IF NOT KEYWORD_SET(nstart) THEN nstart=10
    IF NOT KEYWORD_SET(nmap) THEN nmap=50
    IF NOT KEYWORD_SET(nsubcycle) THEN nsubcycle=1
    ;if starting points are too close to boundary, field lines may exit domain
    delta_x=MAX(g.psixy)-MIN(g.psixy)
    IF NOT KEYWORD_SET(xmin) THEN xmin=MIN(g.psixy)+0.05*delta_x
    IF NOT KEYWORD_SET(xmax) THEN xmax=MAX(g.psixy)-0.05*delta_x
    ;call poincare_main 
    Poincare_Main, tindex=tindex, time=tstring, nstart=nstart, nmap=nmap,      $
                   nsubcycle=nsubcycle, xmin=xmin, xmax=xmax, savefile=s_file, $
                   debug=debug
  ENDIF ELSE s_file=path+'/'+s_file

;-------------------PLOTTING BEGINS HERE----------------------------------------
  ;restore the data from the savefile (Poincare_Main produces no output)
  IF quiet NE 2 THEN BEGIN
    RESTORE,s_file
  ;a quick fix to plot over all computational domain, not just domain where
  ;we started the poincare points
;  xminplot=-0.482359
;  xmaxplot= 0.257340
    loadct,39
    device,decomposed=0
    tek_color
    PLOT,[zmin,zmax],[xminplot,xmaxplot], title="time = "+time ,xtit='z [rad]',$
        ytit='x [weber]', /nod, chars=2,/xst,/yst,color=color
    nPts=N_ELEMENTS(allPts.x)
    FOR i=0l,nPts-1 DO BEGIN
;      PLOTS, allPts.z[i], allPts.x[i], col=allPts.c[i], psym=3
      PLOTS, allPts.z[i],allPts.x[i],col=allPts.c[i],psym=6,symsize=0.05
    ENDFOR
;   Draw a second, logarithmic axis on the right-hand side of the plot.
;    AXIS,YAXIS=1,YRANGE=[4.25,4.85],/SAVE,ystyle=1,ytit='R[m]',color=color
  ENDIF
;--------------------PLOTTING ENDS HERE-----------------------------------------

END
