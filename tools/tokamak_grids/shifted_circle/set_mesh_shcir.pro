pro set_rho_line, j, theta_now, rhoarr, deltarr, psiarr, rm, zm, psi, br, bz, bpol, bphi, b, $
nodetype=nodetype
;
; Set data on a line theta=const for a given node type:
; 0-center, 1-sw, 2-se, 3-nw, 4-ne in (theta,psi) coords
;----------------------------------------------------------------;

     rz=get_RZ(deltarr, rhoarr, theta_now, thetag=thetag)
     rm[j,*,nodetype]=rz.r
     zm[j,*,nodetype]=rz.z


     ;-set components of magnetic field
     bb=get_B(deltarr,rhoarr,thetag)
      bpol[j,*,nodetype]=bb.bp
      bphi[j,*,nodetype]=bb.bt
      br[j,*,nodetype]=bb.br
      bz[j,*,nodetype]=bb.bz
      b[j,*,nodetype]=bb.btot

      ;-set poloidal magnetic flux
      psi[j,*,nodetype]=psiarr

end


pro set_mesh_shcir, d1, d2, d3, plotmesh=plotmesh, export=export, noreset=noreset, $
                    Nrho=Nrho, Ntheta=Ntheta, thetaMin=thetaMin, thetaMax=thetaMax, $
                    rhoMin=rhoMin, rhoMax=rhoMax,DEBUG=DEBUG
;
;-Set BOUT mesh for the shifted circle case
;----------------------------------------------------------------;

 if not keyword_set(thetaMin) then thetaMin=0.
 if not keyword_set(thetaMax) then thetaMax=2*!PI

;-theta mesh (this is the new (orthogonal) theta)
 if not keyword_set(Ntheta) then Ntheta=32
 Ntheta_wg=Ntheta+2 ;-with guards
 
   eps=1e-3 ;-poloidal size of guard cells 
   delt=(ThetaMax-ThetaMin-2*eps)/Ntheta ;-poloidal size of regular cells
   theta=(ThetaMin+eps+delt/2) + delt*findgen(Ntheta)
   dThe=fltarr(Ntheta)+0.5*delt

   ;-add guard cells
   theta=[ThetaMin+eps/2, theta, ThetaMax-eps/2]
   dThe=[eps/2,dThe,eps/2]



;-rho mesh
 if not keyword_set(Nrho) then Nrho=20
 Nrho_wg=Nrho+2 ;-with guards
 if not keyword_set(rhoMin) then rhoMin=0.75
 if not keyword_set(rhoMax) then rhoMax=0.95
 drho=0.5*(rhomax-rhomin)/(Nrho_wg-1) ;-half cell size

 ;-rho for cell centers and corners up and down in rho
 rhoarrc=rhomin + (rhomax-rhomin)*findgen(Nrho_wg)/(Nrho_wg-1)  
 rhoarrp=rhoarrc+drho
 rhoarrm=rhoarrc-drho


;-quantities defined on cell centers and corners
 rm=fltarr(Ntheta_wg,Nrho_wg,5)
 zm=fltarr(Ntheta_wg,Nrho_wg,5)
 psi=fltarr(Ntheta_wg,Nrho_wg,5)
 br=fltarr(Ntheta_wg,Nrho_wg,5)
 bz=fltarr(Ntheta_wg,Nrho_wg,5)
 bpol=fltarr(Ntheta_wg,Nrho_wg,5)
 bphi=fltarr(Ntheta_wg,Nrho_wg,5)
 b=fltarr(Ntheta_wg,Nrho_wg,5)

;-quantities defined on cell centers only
 ni=fltarr(Ntheta_wg,Nrho_wg)
 te=fltarr(Ntheta_wg,Nrho_wg)
 ti=fltarr(Ntheta_wg,Nrho_wg)
 up=fltarr(Ntheta_wg,Nrho_wg)
 phi=fltarr(Ntheta_wg,Nrho_wg)
 dx=fltarr(Ntheta_wg,Nrho_wg)


;-Calculate the Shafranov shift for centers
  Get_delta, deltarr=deltarrc, psiarr=psiarrc, rhoarr=rhoarrc, noreset=noreset, qarr=qarr

;-Calculate the Shafranov shift for corners up in rho
  Get_delta, deltarr=deltarrp, psiarr=psiarrp, rhoarr=rhoarrp, noreset=noreset 

;-Calculate the Shafranov shift for corners down in rho
  Get_delta, deltarr=deltarrm, psiarr=psiarrm, rhoarr=rhoarrm, noreset=noreset




   FOR j=0,Ntheta_wg-1 do begin

     print, 'j=',j

     theta_now=theta[j]
     
      ;-set values at centers and corners
      for nodetype=0,4 do begin

             case nodetype of
               0: begin
                   rhoarr=rhoarrc
                   deltarr=deltarrc
                   psiarr=psiarrc
                   theta_now=theta[j]
                  end

               1: begin
                   rhoarr=rhoarrm
                   deltarr=deltarrm
                   psiarr=psiarrm
                   theta_now=theta[j]-dthe[j]
                  end

               2: begin
                   rhoarr=rhoarrm
                   deltarr=deltarrm
                   psiarr=psiarrm
                   theta_now=theta[j]+dthe[j]
                  end

               3: begin
                   rhoarr=rhoarrp
                   deltarr=deltarrp
                   psiarr=psiarrp
                   theta_now=theta[j]-dthe[j]
                  end

               4: begin
                   rhoarr=rhoarrp
                   deltarr=deltarrp
                   psiarr=psiarrp
                   theta_now=theta[j]+dthe[j]
                  end
              endcase

             Set_Rho_Line, j, theta_now, rhoarr, deltarr, psiarr, $
                      rm, zm, psi, br, bz, bpol, bphi, b, nodetype=nodetype             

             ;STOP
      endfor


     ;-set density, temperatures, and parallel velocity at cell centers
     pressure=npres(rhoarrc,/cgs)
     te[j,*]=100. ;-[eV]
     ti[j,*]=1. ;-[eV]
     ni[j,*]=pressure/((te[j,*]+ti[j,*])*1.6e-12);-[cm-3]
     up[j,*]=0.0
     phi[j,*]=0.0


     ;-calculate radial diameter of cell
     rbot=0.5*(rm[j,*,1]+rm[j,*,2])
     zbot=0.5*(zm[j,*,1]+zm[j,*,2])
     rtop=0.5*(rm[j,*,3]+rm[j,*,4])
     ztop=0.5*(zm[j,*,3]+zm[j,*,4])
     dx[j,*]=sqrt((rtop-rbot)^2+(ztop-zbot)^2)
     

   ENDFOR

;-convert to MKS to simulate UEDGE output

  d1={$
      rm_com:rm*1e-2,$
      zm_com:zm*1e-2,$
      psi_com:psi*1e-8,$
      br_com:br*1e-4,$
      bz_com:bz*1e-4,$
      bpol_com:bpol*1e-4,$
      bphi_com:bphi*1e-4,$
      b_com:b*1e-4,$
      qarr:qarr, $
      nx_com:Ntheta,$
      ny_com:Nrho,$
      ixpt1_com:0,$
      ixpt2_com:Ntheta,$
      nxpt_com:0,$
      ixlb_com:0,$
      ixrb_com:Ntheta,$
      iysptrx1_com:Nrho,$ 
      iysptrx2_com:Nrho,$
      ix_lim_com:0,$
      runidg_grd:'shifted circle'}


    d2={$
      gy_com:1./(dx*1e-2),$ 
      xcs_com:0.0,$ ;-never used
      yyc_com:rm[*,0],$ 
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
      ni___1__value:ni*1e6,$
      up___1__value:up*1e-2,$
      te____ev_value:te,$
      ti____ev_value:ti,$
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

     for i=0,Ntheta_wg-1 do begin 
      for j=0,Nrho_wg-1 do begin 
        oplot, rm[i,j,ind],zm[i,j,ind]
      endfor 
     endfor

   ENDIF


   IF keyword_set(EXPORT) then begin
    if not keyword_set(PATH) then PATH='.' 
    
      file1=PATH + '/gridue.nc' 
      file2=PATH + '/uedgegrd.nc'
      file3=PATH + '/uedgeout.nc'

      print, 'Writing files...'
      print, file1  
      print, file2  
      print, file3  
      
      status = file_export( file1, d1 )
      status = file_export(file2, d2)
      status = file_export(file3, d3)
   ENDIF

   if keyword_set(DEBUG) then STOP
end
