; Generate the grid file for a straight linear device like LAPD. 
;
;
; IDL> set_mesh_cyl, /export, Nr=54, Nz=64, rMin=0.15, rMax=0.45, $
;                    ni0=2.5e18, te0=5., Bz0=0.08, ni_profile_type=1, $
;                    phi_profile_type=0, phi0V=0.0, /NOPLOTS
;
; -> Outputs "gridue.nc", "uedgegrd.nc" and "uedgeout.nc"
;
; IDL> !path=!path+":../tokamak_grids/all/"
; IDL> uedge2bout
;
; Add vacuum region? N
; Equilibrium correction option: 0 (no correction)
; Use new hthe? N
;
; -> Output "uedge.grd.nc" input grid for BOUT++
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; for the Helimak
; 
;set_mesh_cyl,/export,Nr  = 68,Nz = 64,rMin = .7, rMax = 1.5,ni0 =2.0e17,te0=10.0,Bz0 = .01,bphi0 = .1,Zmax=2.0,ti0 =1.0,ni_profile_type =1,te_profile_type =1 , ti_profile_type= 1
;read_uedata3, /s, d, /noref, /NOPLOTS,filename =
;"helimakv2.nc"

;small Helimak grid
;set_mesh_cyl,/export,Nr  = 20,Nz = 16,rMin = .7, rMax = 1.5,ni0
;=2.0e17,te0=10.0,Bz0 = .01,bphi0 = .1,Zmax=2.0,ti0
;=1.0,ni_profile_type =1,te_profile_type =1 , ti_profile_type= 1
;read_uedata3, /s, d, /noref, /NOPLOTS,filename="helimak_16x16_fixBxy.nc"

; for the CLM
; set_mesh_cyl,/export,Nr  = 68,Nz = 64,rMin = 0.001, rMax = .03,ni0 =
; 3.0e15,te0=15.,Bz0 = .10,bphi0 = 0,Zmax=6.0,ni_profile_type = 0,ti0 =
; 1.0, te_profile_type = 'CLM'
;read_uedata3, /s, d, /noref, /NOPLOTS,filename = "CLM.nc"

;small CLM grid for testing
;set_mesh_cyl,/export,Nr = 20, Nz = 16,rMin = 0.001, rMax = .03,ni0 =
;3.0e15,te0=15.,Bz0 = .10,bphi0 = 0,Zmax=6.0,ni_profile_type = 0,ti0 =
; 1.0, te_profile_type = 'CLM'

pro CLM_grid,filename=filename
  set_mesh_cyl,/export,Nr = 36, Nz = 32,rMin = 0.001, rMax = .03,ni0 = $
               3.0e15,te0=15.,Bz0 = .10,bphi0 = 0,Zmax=6.0,$
               ni_profile_type = 0,ti0 =1.0, te_profile_type = 1,$
               phi0V =5./10000., phi_profile_type = 0

  if not keyword_set(filename) then filename = "CLM.nc"
  read_uedata3, /s, d, /noref, /NOPLOTS,filename = filename
  spawn,"rm *.pdb"
  
end

pro Helimak_grid,expr_prof = expr_prof,grid_size = N,filename = filename

  
  restore,"helimak.sav"

  expr_data = (shot_data.set3)[5]


  if (not keyword_set(grid_size)) then N = 4
  
  Nz = 2^N
  Nr = 2^N + 4
  
  temp = ["Helimak_",string(2^N),"x",string(2^N),".nc"]
  if (not keyword_set(filename)) then filename = strcompress(strjoin(temp),/remove_all)

  set_mesh_cyl,/export,Nr = Nr, Nz = Nz,rMin = 0.7, rMax = 1.5,ni0 = $
               1.0e17,te0=10.,Bz0 = .01,bphi0 = .1,Zmax=2.0,$
               ni_profile_type = 1,ti0 =1.0, te_profile_type = 2,$
               ti_profile_type = 1,phi_profile_type = 0,expr_prof = expr_prof
  read_uedata3, /s, d, /noref, /NOPLOTS,filename = filename
  spawn,"rm *.pdb"
end

pro simple_cyl,grid_size = N,filename = filename
  set_mesh_cyl,/export,Nr = 36, Nz = 32,rMin = 0.001, rMax = .03,ni0 = $
               3.0e15,te0=15.,Bz0 = .10,bphi0 = 0,Zmax=6.0,$
               ni_profile_type = 0,ti0 =1.0, te_profile_type = 1,$
               phi0V =5./10000., phi_profile_type = 0

  if not keyword_set(filename) then filename = "CLM.nc"
  read_uedata3, /s, d, /noref, /NOPLOTS,filename = filename
  spawn,"rm *.pdb"
end

; ---------------------------------------------------------------------
function int3_gaussian, x, w
;  int3_delta_approx(x) Third integral of Gaussian (approximate delta-function)
;                       w -- width parameter

res = 0.125*( 2.*exp(-(x/w)^2)*w*x/sqrt(!PI) + (w*w+2.*x*x)*erf(x/w)) 

return, res

end
; ---------------------------------------------------------------------

function interp_and_smooth, input, x,r
; simple function to 
  temp = input/max(abs(input))
  output = interpol(smooth(temp,N_ELEMENTS(x)/9.0),x,r)
  return, output
end


function ni_radial_profile, r, p
; Density profile
; r in [m]


case p.ni_profile_type of ; Choose the density profile

1: begin ; Simple linear profile for debugging
   x = (r-p.ra)/(p.rb-p.ra)     ; 0<=x<=1
   ni = 1. - 0.5*x
end

2: begin                        ; Simple exp profile with some decay slope
   x = (r-p.ra)       ; 0<=x<=1
   ni = exp(-x/p.lam_n) ;lam_n must be given in meters here
   print,p.lam_n
end

3: begin                        ; Simple linear profile for debugging
   x = (r-p.ra)     ;unnormalized
   ni = 1. - p.slope_n*x
end

4: begin                        ; Simple linear profile for debugging
   

   x = (r-p.ra)                 ;unnormalized
   chi = 15+175*(.1-p.lam_n) ;only valid near lam_n ~= .1 m
   print,r,chi,p.lam_n
   ni = (atan(chi *(r-.9))+!PI/2.)/(atan(chi *(p.rb-.9))+!PI/2.)
   
end

9: begin
   ;import a custom profile
   
   ni = interp_and_smooth(p.expr_prof.density,p.expr_prof.R,r)
   end

else: ni=1.

endcase

return, ni
end

; ---------------------------------------------------------------------

function phi_radial_profile, r, p

case p.phi_profile_type of ; Choose phi profile

;----------------------------------------

0: begin ; Quadratic in r -> const omegaExB ( Er ~ d phi0/dr, omegaExB ~ 1/r Er )
   phi = (r/p.rb)^2   ; phi(rmax)=1
   end

1: begin
   phi = ((r-p.ra)/(p.rb-p.ra))^2   ; phi(rmax)=1
   end

2: begin
   x = (r-p.ra)/(p.rb-p.ra)*4  ; 0<=x<=4
   if (x le 1) then phi = 0.
   if ((x ge 1)and(x le 2)) then phi = (x-1.)^2
   if ((x ge 2)and(x le 3)) then phi = 2-(x-3.)^2
   if (x ge 3) then phi = 2.
   end
9: begin
   ;import a custom profile
   x = p.expr_prof.R
   phi= interp_and_smooth(p.expr_prof.vfloat,p.expr_prof.R,r)
   
   end

;----------------------------------------


else: phi=1.

endcase

return, phi
end

; ---------------------------------------------------------------------

function te_radial_profile, r, p

case p.te_profile_type of ; Choose Te profile

0: te = 1.

1: begin ; CLM profile
   ;print,r
   te = .5+ .48 *tanh(-(r-.0181)/.0033)
   help,r,te
end 
2: begin                        ; Simple exp profile with some decay slope
   x = (r-p.ra)       ; 0<=x<=1
   te = exp(-x/p.lam_e)
end
3: begin                        ; Simple linear profile for debugging
   x = (r-p.ra)     ;unnormalized
   print,"slope_te: ",p.slope_te
   te = 1. - p.slope_te*x
end


9: begin
   ;import a custom profile
   
   te = interp_and_smooth(p.expr_prof.te,p.expr_prof.R,r)
   end
else: te=1.

endcase

return, te
end



function ti_radial_profile, r, p

case p.ti_profile_type of ; Choose Te profile

0: ti = 1.

1: begin ; Simple linear profile for debugging
   x = (r-p.ra)/(p.rb-p.ra)     ; 0<=x<=1
   ti = 1. - 0.5*x
end


2: begin                        ; Simple exp profile with some decay slope
   x = (r-p.ra)       ; 0<=x<=1
   ti = exp(-x/p.lam_i)
end
3: begin
   x = (r-p.ra)
   ti = 1. - p.slope_ti*x
end
9: begin
   ;import a custom profile
   
   ti = interp_and_smooth(p.expr_prof.ti,p.expr_prof.R,r)
   end

else: ti=1.

endcase

return, ti
end
; ---------------------------------------------------------------------

pro set_mesh_cyl, d1, d2, d3, plotmesh=plotmesh, export=export, format=format,     $
                  noreset=noreset,   $
                  Nr=Nr, Nz=Nz,                                                    $
                  rMin=rMin, rMax=rMax, zMin=zMin, zMax=zMax,                      $
                  Bz0=Bz0, bphi0=bphi0,                                            $

                  ; Parameters for profiles (density, potential, temperature)
                  ra=ra, rb=rb,  $ ; common parameters for all profiles

                  ni0=ni0,  ni_profile_type=ni_profile_type,      $
                            pna=pna, pnb=pnb, slope_n = slope_n,                    $
                  
                  phi0V=phi0V, phi_profile_type=phi_profile_type, $
                            ppa=ppa, ppb=ppb,                     $

                  te0=te0,  te_profile_type=te_profile_type,      $
                            pta=pta, ptb=ptb,slope_te=slope_te,  $

                  ti0=ti0,   ti_profile_type=ti_profile_type,     $
                  slope_ti=slope_ti, $
                  nn0=nn0,       $
                  iysptrx=iysptrx, comment=comment, DEBUG=DEBUG, $
                  NOPLOTS=NOPLOTS,lam_i= lam_i, lam_n = lam_n, $
                  lam_e=lam_e,expr_prof = expr_prof
;
;-Set BOUT mesh for the shifted circle case
;----------------------------------------------------------------;

 if not keyword_set(comment) then comment=''

 IF NOT KEYWORD_SET(format) THEN format = "nc"  ; File output format

; Tesla
 if not keyword_set(Bz0) then Bz0=0.1
 if not keyword_set(bphi0) then bphi0=1e-10


 if not keyword_set(Nr) then Nr=20
 if not keyword_set(Nz) then Nz=30  ; number of intervals, open-end periodic conditions


; Plasma parameters
 if not keyword_set(te0) then te0=5.
 if not keyword_set(ti0) then ti0=1.e-2
 if not keyword_set(ni0) then ni0=1.e18 ;[m^-3]
 if not keyword_set(nn0) then nn0=1.e18 ;[m^-3]
 if not keyword_set(phi0V) then phi0V=0.  ; [V]

 

; When iysptrx=Nr:  periodic boundary condition
; When iysptrx=0:   sheath plate, material wall condition
 if not keyword_set(iysptrx) then iysptrx=Nr


; MKS, meters
 if not keyword_set(rMin) then rMin=0.1
 if not keyword_set(rMax) then rMax=1.
 if not keyword_set(zMin) then zMin=0.
 if not keyword_set(zMax) then zMax=17.
 
 if KEYWORD_SET(expr_prof) then begin
    rMin = min(expr_prof.R)
    rMax = max(expr_prof.R)
    
 endif

; Radial profiles, parameters common for all profiles
 if not keyword_set(ra) then ra=rMin
 if not keyword_set(rb) then rb=rMax


; Radial density profile, gradient region
 if not keyword_set(ni_profile_type) then ni_profile_type=0
 if not keyword_set(pna) then pna=5.  ; steepness of ni profile for ni_profile_type=2
 if not keyword_set(pnb) then pnb=0.
; if not keyword_set(N0_NLD) then N0_NLD=4.

; Radial phi0 profile
 if not keyword_set(phi_profile_type) then phi_profile_type=0
 if not keyword_set(ppa) then ppa=0.
 if not keyword_set(ppb) then ppb=0.

; Radial te0 profile
 if not keyword_set(te_profile_type) then te_profile_type=0 ; const Te
 if not keyword_set(pta) then pta=0.
 if not keyword_set(ptb) then ptb=0.


 if not keyword_set(ti_profile_type) then ti_profile_type=0 ; const Te
 
 if not keyword_set(lam_i) then lam_i = 1
 if not keyword_set(lam_e) then lam_e = 1
 if not keyword_set(lam_n) then lam_n = 1

 if not keyword_set(slope_te) then slope_te = 1
 if not keyword_set(slope_ti) then slope_ti = 1
 if not keyword_set(slope_n) then slope_n = 1
 

 eit_i = lam_n/lam_i
 eit_e = lam_n/lam_e

 if KEYWORD_SET(expr_prof) then begin
    which_p = tag_names(expr_prof)
    
    profs = ['r','vfloat','Te','Ti','density']
    prof_tgl = intarr(5)
    ;check for R,Te,Ti,phi and n profiles
     ; phi_profile_type = 9
    
;read in the avaliable profiles
    for i =0,N_ELEMENTS(profs)-1 do begin
       pro_i = WHERE(STRMATCH(which_p, profs[i], /FOLD_CASE) EQ 1)
       if pro_i NE -1 then begin 
          temp = which_p[pro_i]
          print,'importing experimental ',temp,'profiles'
          prof_tgl[i] = 1

       endif
    endfor

                                ;overwrite the profile type setting for all valid imported profiles
    

    r_imp = prof_tgl[0]
    phi_imp = prof_tgl[1]
    Te_imp = prof_tgl[2]
    Ti_imp = prof_tgl[3]
    ni_imp = prof_tgl[4]

                                ;we need to allow the user the user to
                                ;set the radial resolution with Nr
                                ;keyword

   ;;  if keyword_set(expr_prof) then begin
;;     te0 = max(expr_prof.Te)
;;     ti0 = max(expr_prof.Ti)
;;     ni0 = max(expr_prof.density)
;;     phi0V = max(abs(expr_prof.vfloat))
;;  endif


    if r_imp then r_profile_type = 9
    if ni_imp then begin
       ni_profile_type = 9
       ni0 = max(expr_prof.density)
    endif
    if phi_imp then begin
       phi_profile_type = 9 
       phi0V = max(expr_prof.vfloat)
    endif                                   
    if Te_imp then begin
       Te_profile_type = 9
       te0 = max(expr_prof.Te)
    endif
    if Ti_imp then begin
       Ti_profile_type = 9
       ti0 = max(expr_prof.Ti)
    endif
    

    ;set 
 endif  else expr_prof = 0

 

params = {$; set a structure with all radial profile parameters
         ni_profile_type   : ni_profile_type,          $
         phi_profile_type  : phi_profile_type,         $
         te_profile_type   : te_profile_type,$
         ti_profile_type   : ti_profile_type,$  
         ra: ra,   rb: rb,                             $
         pna:pna,  pnb:pnb,                            $
         ppa:ppa,  ppb:ppb,                            $
         pta:pta,  ptb:ptb,                            $
         lam_i:lam_i, lam_e: lam_e, lam_n:lam_n,  $
         slope_te:slope_te,slope_ti:slope_ti,$
         slope_n:slope_n,$
         expr_prof:expr_prof $
}
; lam_x is in meters or should be

;-theta mesh (this is the new (orthogonal) theta)
 Nz_wg=Nz+2 ;-with guards
 Nr_wg=Nr+2 ;-with guards

;-quantities defined on cell centers and corners
 rm=fltarr(Nz_wg,Nr_wg,5)
 zm=fltarr(Nz_wg,Nr_wg,5)

 ;magnetic field quantities
 psi=fltarr(Nz_wg,Nr_wg,5)
 br=fltarr(Nz_wg,Nr_wg,5)
 bz=fltarr(Nz_wg,Nr_wg,5)
 bpol=fltarr(Nz_wg,Nr_wg,5)
 bphi=fltarr(Nz_wg,Nr_wg,5)
 b=fltarr(Nz_wg,Nr_wg,5)

;-quantities defined on cell centers only
 ni=fltarr(Nz_wg,Nr_wg)
 te=fltarr(Nz_wg,Nr_wg)
 ti=fltarr(Nz_wg,Nr_wg)
 nn=fltarr(Nz_wg,Nr_wg)
 up=fltarr(Nz_wg,Nr_wg)
 phi=fltarr(Nz_wg,Nr_wg)
 dx=fltarr(Nz_wg,Nr_wg)

; Cell size (across)
 dr = (rMax-rMin)/(Nr-1)
 dz = (zMax-zMin)/(Nz)

 for jz=0, Nz_wg-1 do begin ; along the field line
         
   for ir=0, Nr_wg-1 do begin ;radial
     rm[jz,ir,0]=rMin + (rMax-rMin)*(ir-1)/(Nr-1)
     zm[jz,ir,0]=zMax - dz*(jz-0.5)
     
     rm[jz,ir,1]=rm[jz,ir,0]-dr/2
     zm[jz,ir,1]=zm[jz,ir,0]-dz/2

     rm[jz,ir,2]=rm[jz,ir,0]-dr/2
     zm[jz,ir,2]=zm[jz,ir,0]+dz/2

     rm[jz,ir,3]=rm[jz,ir,0]+dr/2
     zm[jz,ir,3]=zm[jz,ir,0]-dz/2

     rm[jz,ir,4]=rm[jz,ir,0]+dr/2
     zm[jz,ir,4]=zm[jz,ir,0]+dz/2


     bz[jz,ir,*]=Bz0
     br[jz,ir,*]=0.
     ;bphi[jz,ir,*]=bphi0
     print,(rMin + (rMax-rMin)/2.0)/(rm[jz,ir,*])
     
     bphi[jz,ir,*]=bphi0 * (rMin + (rMax-rMin)/2.0)/(rm[jz,ir,*]) ;1/R falloff
     bpol[jz,ir,*]=Bz0
     psi[jz,ir,*]=Bz0*(rm[jz,ir,*]^2-rMin^2)/2. ;makes sense

     ; Plasma profiles
     te[jz,ir]=te0 * te_radial_profile(rm[jz,ir,0],params)   ; eV
     ti[jz,ir]=ti0 * ti_radial_profile(rm[jz,ir,0],params)   ; eV
     ni[jz,ir]=ni0 * ni_radial_profile(rm[jz,ir,0],params)  
                    ; m^-3, L in meters

     nn[jz,ir]=nn0 ;[m^-3]

     phi[jz,ir]=phi0V * phi_radial_profile(rm[jz,ir,0],params) ; V
     up[jz,ir]=0.  ; m/s

     dx[jz,ir]=1.

   endfor
endfor



if not keyword_set(NOPLOTS) then begin
   !p.multi=[0,3,1,0,0]
   plot, ni[1,*], title='Density profile, m^-3', charsize=2
   plot, te[1,*], title='Te profile, eV', charsize=2
   plot, phi[1,*], title='Phi0 profile, V', charsize=2
   WAIT, 2.
   !p.multi=0
endif


;b=sqrt(br^2+bz^2) what about the phi(toroidal) component
b=sqrt(br^2+bz^2+bphi^2)

 if keyword_set(DEBUG) then stop

;-convert to MKS to simulate UEDGE output

  d1={$
      rm_com:rm,$
      zm_com:zm,$
      psi_com:psi,$
      br_com:br,$
      bz_com:bz,$
      bpol_com:bpol,$
      bphi_com:bphi,$
      b_com:b,$
      nx_com:Nz,$
      ny_com:Nr,$
      ixpt1_com:0,$
      ixpt2_com:Nz,$
      nxpt_com:0,$
      ixlb_com:0,$
      ixrb_com:Nz,$
      iysptrx1_com:iysptrx,$ 
      iysptrx2_com:iysptrx,$
      ix_lim_com:0,$
      runidg_grd:'  cylinder  '+comment}


    d2={$
      gy_com:1./(dx*1e-2),$ 
      xcs_com:0.0,$ ;-never used
      yyc_com:REFORM(rm[0,*,0]),$ 
      rmagx_flx:0.0,$ ;-never used       
      rseps_flx:0.0,$ ;-never used 
      zmagx_flx:0.0,$ ;-never used
      zseps_flx:0.0,$ ;-never used
      psi0max_flx:0.0,$ ;-never used
      psi0min1_flx:0.0,$ ;-never used
      psi0min2_flx:0.0,$ ;-never used
      sibdryg_grd:max(psi),$
      simagxg_grd:0.0}  
      
      
    d3={$
      ni___1__value:ni,$
      up___1__value:up,$
      te____ev_value:te,$
      ti____ev_value:ti,$
      nn___1__value:nn,$
      phi____value:phi}   


;   IF keyword_set(PLOTMESH) then begin
;     plot, rm[*,*,0], zm[*,*,0],psym=4,/iso
;     for i=0,Nrho_wg-1 do begin 
;        oplot, rm[*,i,0],zm[*,i,0]
;     endfor
;     for j=0,Ntheta_wg-1 do begin 
;        oplot, rm[j,*,0],zm[j,*,0]
;     endfor
;   ENDIF


   IF keyword_set(PLOTMESH) then begin
     plot, rm[*,*,0], zm[*,*,0],psym=3,/iso
     ind=[1,2,4,3,1] 

    for ir=0,Nr_wg-1 do begin 
     for jz=0,Nz_wg-1 do begin 
        oplot, rm[jz,ir,ind],zm[jz,ir,ind]
      endfor 
     endfor

   ENDIF


   IF keyword_set(EXPORT) then begin
     if not keyword_set(PATH) then PATH='.' 
     
     file1=PATH + '/gridue.'+format 
     file2=PATH + '/uedgegrd.'+format
     file3=PATH + '/uedgeout.'+format
     
     print, 'Writing files...'
     print, file1  
     print, file2  
     print, file3  
     
     status = file_export(file1, d1)
     status = file_export(file2, d2)
     status = file_export(file3, d3)

     print,"Done exporting"
   ENDIF

if keyword_set(DEBUG) then STOP
end
