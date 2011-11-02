;.com wavelet,find_phase_ND,post_bout,basic_info,fftrot

function get_var,bout_str,var_name,default_val,string = string,$
                 subfield = subfield


  n_lines = (size(bout_str))[2]

  ;by default search the entire file
  ln_start = 0
  ln_stop = n_lines-1

  ;we can restrict the parameter value to certain subfield 
  if KEYWORD_SET(subfield) then begin
     sub_i = wc_where(bout_str[0,*],'[*]')
     sub_j = wc_where(bout_str[0,*],subfield)
     print,subfield,N_ELEMENTS(sub_j) ,sub_j
     

     ; if more that one or none subfield found return default value
     if N_ELEMENTS(sub_j) NE 1 then return, default_val
     if sub_j EQ -1 then return, default_val
     
     ln_start = sub_j

     next_i = where(sub_i GT round(dblarr(N_ELEMENTS(sub_i))+sub_j[0]))

     if next_i[0] NE -1 then ln_stop = sub_i[next_i[0]]-1 else ln_stop = n_lines-1
    
     ;val = (bout_str[2,wc_where(bout_str[0,sub_i: sub_i +3],var_name)])[0]
     
  endif 



                                ;if the relevant field is present in
                                ;the relevant range create a switch,
                                ;if not just return the default value 
  if (wc_where(bout_str[0,ln_start:ln_stop],var_name))[0] NE -1 then begin
     ;print,var_name
     choice = KEYWORD_SET(string) + 2* (wc_where(bout_str[0,ln_start:ln_stop],var_name) NE -1)

  endif else return, default_val

  print,var_name
  print,choice
  
  CASE choice OF
     2: begin
        val = (bout_str[2,ln_start+wc_where(bout_str[0,ln_start:ln_stop],var_name)])[0]
        ;; print,ln_start,ln_stop,val
        ;print,bout_str[*,ln_start+wc_where(bout_str[0,ln_start:ln_stop],var_name)]
        ;print,val
        ;;var = double(val)

;;         print,bout_str[2,ln_start+wc_where(bout_str[0,ln_start:ln_stop],var_name)]

        ;if ln_stop-ln_start LT 5 then begin
           ;print,bout_str[*,ln_start+wc_where(bout_str[0,ln_start:ln_stop],var_name)]
           
        ;endif
        
        if val EQ 'true' then var = 1
        if val EQ 'false' then var = 0
    
        if (val NE 'true' AND val NE 'false')  then var =  double(val)
     
     end

     3: begin
        var =  (bout_str[2,wc_where(bout_str[0,ln_start:ln_stop],var_name)])[0]
       
        var = strsplit(var,'"',/extract)
     end
     1: var = default_val
     0: var = default_val
  endcase
  
  print,var
  return,var
end

function comment, name,value
;out = [name,": ",strcompress(string(value))]
out = [name,": ",string(value)]
out = strjoin(out)
return, out
END


;example of typical usage:
; a = postbout(/clndr,path='/media/TACCWORK/meyerson/Helimak/data',/profiles)

function post_bout, tind = tind,cases = cases,path = path,$
                    save_dir=save_dir,clndr = clndr,$
                    profiles = profiles,remoteFS = remoteFS,$
                    rotate = rotate

  if KEYWORD_SET(path) EQ 0 then begin
     bout_inp = "BOUT.inp"
  endif else begin
     bout_inp = strcompress(strjoin([path,"/BOUT.inp"]),/remove_all)
  endelse
  

print,bout_inp

bout_inp = rd_tfile(bout_inp,/auto)

TWOPI = 2.0 * !PI

AA = get_var(bout_inp,'AA',2.0)   ; mass
ZZ = get_var(bout_inp,'ZZ',1.0)   ; charge
TIMESTEP = get_var(bout_inp,'TIMESTEP',1.0)
zeff = get_var(bout_inp,'zeff',1.0)
nu_perp = get_var(bout_inp,'nu_perp',0.0)
ShearFactor = get_var(bout_inp,'ShearFactor',1.0)
estatic = get_var(bout_inp,'estatic',0)
monitor = get_var(bout_inp,'monitor',0)


evolve_Ajpar = get_var(bout_inp,'evolve',1,subfield = '[Ajpar]')
evolve_Ni= get_var(bout_inp,'evolve',1,subfield = '[Ni]')
evolve_Vi = get_var(bout_inp,'evolve',1,subfield = '[Vi]')
evolve_Ti = get_var(bout_inp,'evolve',1,subfield = '[Ti]')
evolve_Te = get_var(bout_inp,'evolve',1,subfield = '[Te]')
evolve_rho = get_var(bout_inp,'evolve',1,subfield = '[rho]')

evolve_i = (findgen(6)+1)*[evolve_Ajpar,evolve_Vi,evolve_Ni,$
                           evolve_Ti,evolve_Te,evolve_rho]
print,evolve_i
print,"evolve: ",[evolve_Ajpar,evolve_Ni,evolve_Vi,evolve_Ti,evolve_Te,evolve_rho]

;return,0

evolve_i = evolve_i[where(evolve_i NE 0)] -1

print,evolve_i

evolve = (['Ajpar ','Vi ','Ni ','Ti ','Te ','rho '])[evolve_i]
print,evolve
;return, 0

MXG = get_var(bout_inp,'MXG',0)
MYG = get_var(bout_inp,'MYG',2)

ys_opt = get_var(bout_inp,'ys_opt',2)
zs_opt = get_var(bout_inp,'zs_opt',2)
xs_opt = get_var(bout_inp,'xs_opt',2)
opt = round([xs_opt,ys_opt,zs_opt]);

ys_mode = get_var(bout_inp,'ys_mode',0)
zs_mode = get_var(bout_inp,'zs_mode',0)
xs_mode = get_var(bout_inp,'xs_mode',0)

nt = get_var(bout_inp,'NOUT',10)

mode = round([xs_mode,ys_mode,zs_mode])

nz =  get_var(bout_inp,'MZ',0)
grid = (get_var(bout_inp,'grid',0,/string))[0]

if KEYWORD_SET(remoteFS) then begin
  print,"grid: ",grid
  substrings = strsplit(grid,' ',/extract) 
  print,substrings
  substrings = repstr(substrings,'/work/01523','/media/TACCWORK')
  print,'using sshfs to access TACC:',strjoin(substrings)
  grid = strjoin(substrings)
endif

grid_long = grid

;; if not KEYWORD_SET(path) then begin
;;    ;grid_long = strcompress(strjoin(["../../",grid]),/remove_all)
;;    grid_long = strcompress(strjoin(["/../",grid[0]]),/remove_all) 
;; endif else begin
;;    grid_long = strcompress(strjoin([path,"/../",grid]),/remove_all) 
;;    ;;grid_long = grid
;; endelse

d = file_import(grid_long)
d.Ni0 = d.Ni0*1e14



zmin = get_var(bout_inp,'ZMIN',0)
zmax = get_var(bout_inp,'ZMAX',0)

zeta_range = zmax-zmin
dt = get_var(bout_inp,'TIMESTEP',1)

time_x = indgen(nt+1,/double) * dt


d.bmag = d.bmag*1.0e4           ;to cgs
d.Ni_x = d.Ni_x*1.0e14          ;cm^-3
rho_s = 1.02e2*sqrt(AA*d.Te_x)/ZZ/d.bmag   ; ion gyrorad at T_e, in cmsubaru wrx -  
rho_i = 1.02e2*sqrt(AA*d.Ti_x)/ZZ/d.bmag   ; ion gyrorad at T_i
rho_e = 2.38*sqrt(d.Te_x)/d.bmag   ; elec gyrorad at T_e

fmei  = 1./1836.2/AA            ;

lambda_ei = 24.-alog(sqrt(d.Ni_x)/d.Te_x) ;
lambda_ii = 23.-alog(ZZ*ZZ*ZZ*sqrt(2.*d.Ni_x)/(d.Ti_x^1.5)) ;

wci       = 9.58e3*ZZ*d.bmag/AA ; ion gyrofrteq
wpi       = 1.32e3*ZZ*sqrt(d.Ni_x/AA) ; ion plasma freq 

wce       = 1.78e7*d.bmag ;electron gyrofreq
wpe       = 5.64e4*sqrt(d.Ni_x);electron plasma freq

v_the    = 4.19e7*sqrt(d.Te_x);cm/s
v_thi    = 9.79e5*sqrt(d.Ti_x/AA) ;cm/s
c_s      = 9.79e5*sqrt(5.0/3.0 * ZZ * d.Te_x/AA);
v_A      = 2.18e11*sqrt(1.0/(AA * d.Ni_x))

nueix     = 2.91e-6*d.Ni_x*lambda_ei/d.Te_x^1.5 ;
nuiix     = 4.78e-8*ZZ^4.*d.Ni_x*lambda_ii/d.Ti_x^1.5/sqrt(AA) ;
nu_hat    = zeff*nueix/wci      ;
  

L_d       = 7.43e2*sqrt(d.Te_x/d.Ni_x)
L_i_inrt  = 2.28e7*sqrt(AA/d.Ni_x)/ZZ ;ion inertial length in cm
L_e_inrt  = 5.31e5*sqrt(d.Ni_x) ;elec inertial length in cm

wce = 1.76e7*d.bmag 
Ve_x = 4.19e7*d.Te_x

  Vi_x = wci * rho_s 
  Vi_xx = wci * rho_i
  L_n=MEAN(ABS(d.Ni0[*,d.ny/2]/DERIV(d.Rxy[*,d.ny/2]*1e2,d.Ni0[*,d.ny/2]))) ;-Ni scale length [cm]
  L_te=MEAN(ABS(d.Te0[*,d.ny/2]/DERIV(d.Rxy[*,d.ny/2]*1e2,d.Te0[*,d.ny/2]))) ;-Te scale length [cm]
  L_ti=MEAN(ABS(d.Ti0[*,d.ny/2]/DERIV(d.Rxy[*,d.ny/2]*1e2,d.Te0[*,d.ny/2]$
                                     ))) ;-Te scale length [cm]


  norm_s =  get_var(bout_inp,'norm_s','i')
  norm_s = 'i'
  if norm_s EQ 'e' then begin 
     w_n= wce
     v_n= v_the
     rho_n = rho_e
  endif else begin 
     w_n= wci
     v_n= v_thi
     rho_n = rho_i
  endelse

if KEYWORD_SET(clndr) then d.R0 = (max(d.Rxy)+min(d.Rxy))/2.0 

  L_z = 1e2 * TWOPI * d.R0 *(zmax - zmin) ; in cm toroidal range
  print,"L_z: ",L_z
 
  lbNorm=L_z*(d.BPXY[0,d.ny/2]/d.BXY[0,d.ny/2])     ;-binormal coord range [cm]
  zPerp=lbNorm*findgen(nz)/(nz-1)                   ;-binormal coordinate [cm]
  ;L_z = lbNorm
  help,lbNorm


  Btor = mean(d.BTXY[*,d.ny/2])*1e4
  Bpol = mean(d.BPXY[*,d.ny/2])*1e4

  lpar=1e2*total((d.Bxy/d.Bpxy)*d.dlthe)/d.nx ;-[cm], average over flux surfaces, parallel length
                                ;dlthe is just the spacing between the
                                ;points along the y direction;simply
                                ;dy*hthe
 
  ycoord = 1e2*total((d.Bxy/d.Bpxy)*d.dlthe,2,/cumulative)
  
  d = struct_addtags(d,{lpar:lpar})
  d = struct_addtags(d,{ycoord:ycoord})
  d = struct_addtags(d,{Nz:nz,DZ:zeta_range,Nt:nt})

  kpar=2*!pi/(lpar)
  
  dzeta = ((d.Btxy/d.Bpxy)*(d.dlthe/d.Rxy))

  zeta = total((d.Btxy/d.Bpxy)*(d.dlthe/d.Rxy),2,/cumulative) mod (TWOPI)  ;toroidal angle in rads 
 ; zeta = zeta -rebin(zeta[*,0],8,16)

  zeta0 = total((mean(d.Btxy)/d.Bpxy)*(d.dlthe/mean(d.Rxy)),2,/cumulative) mod (TWOPI)  
  d = struct_addtags(d,{zeta:zeta})
  d = struct_addtags(d,{zeta0:zeta0})

;spar=(kpar/kperp)^2 * wci * wce / (0.51 * nuei) ;[1/s]
  
  ;R0 is determined by 
  ; 
  ;not sure about this one - look
  L_y = sqrt((1.0e2*TWOPI *d.R0)^2.0) ; this is the lower bound, y is field aligned, along the covar vector, coord here, 
                                ;for q >0 this is wrong

  L_x = 1.0e2 *(max(d.Rxy)-min(d.Rxy));
   
  L = [L_x,lpar,L_z]
  L_str = {x:L_x,b:lpar,zeta:L_z}
  if KEYWORD_SET(save_dir) NE 0 then begin
     save_dir = save_dir
  endif else save_dir = "./"

  path_i = strsplit(grid,"/",/extract)
  path_n = N_ELEMENTS(path_i)
  if path_n GE 2 then grid = path_i[path_n-1]

  save_name = strcompress(strjoin([save_dir,grid,'.ps']),/remove_all)

 
 

  var_values = [AA,ZZ,zeff,estatic,mode]

  var_values = [{name:'evolved vars ', val:ptr_new(evolve),units:''},$
                {name:textoidl('m_{i}',FONT=-1),val:ptr_new(round(AA)),units:textoidl('m_p')},$
                {name:textoidl('q_{i}'),val:ptr_new(round(ZZ)),units:''},$
                {name:textoidl('Z_{eff}'),val:ptr_new(zeff),units:''},$
                {name:textoidl('estatic'),val:ptr_new(round(estatic)),units:''},$
                {name:textoidl('{n_x, n_y, n_z}'),val:ptr_new(mode),units:''},$
                {name:textoidl('opt'),val:ptr_new(opt),units:''},$
                {name:textoidl('grid'),val:ptr_new(grid),units:''},$
                {name:textoidl('{\phi_{min}, \phi_{max}}'),$
                 val:ptr_new(TWOPI *[zmin,zmax]),units:''},$
                {name:textoidl('\omega_{ci}/\omega_{pi}'),val:ptr_new(wci/wpi),units:''},$
                {name:textoidl('\omega_{ce}/\omega_{pe}'),val:ptr_new(wce/wpe),units:''},$
                {name:textoidl('L/\rho_i'),val:ptr_new(round(L/[rho_i,rho_i,rho_i])),units:''},$
                {name:textoidl('L/\rho_s'),val:ptr_new(round(L/[rho_s,rho_s,rho_s])),units:''},$
                {name:textoidl('\nu_{ei}/\omega_{ci}'),val:ptr_new(nueix/wci),units:''},$
                {name:textoidl('\nu_{ii}/\omega_{ci}'),val:ptr_new(nuiix/wci),units:''},$
                {name:textoidl('B_0'),val:ptr_new(d.bmag),units:'G'},$
                {name:textoidl('Btor and Bpol'),val:ptr_new([Btor,Bpol]),units:'G'},$
                {name:textoidl('\omega_{ci}'),val:ptr_new(wci),units:'Hz'},$
                {name:textoidl('{\omega_{ce}, \omega_{pe}}'),val:ptr_new([wce,wpe]),units:'Hz'},$
                {name:textoidl('Te_x'),val:ptr_new(d.Te_x),units:'eV'}$
                ,{name:textoidl('Ti_x'),val:ptr_new(d.Ti_x),units:'eV'},$
                {name:textoidl('v_{th,e}'),val:ptr_new(v_the),units:'cm/s'},$
                {name:textoidl('v_{th,i}'),val:ptr_new(v_thi),units:'cm/s'},$
                {name:textoidl('c_{s,i}'),val:ptr_new(c_s),units:'cm/s'},$
                {name:textoidl('v_{A}'),val:ptr_new(v_A),units:'cm/s'},$
                {name:textoidl('Ni_x'),val:ptr_new(d.Ni_x),units:textoidl('cm^{-3}')},$
                {name:textoidl('\DeltaT'),val:ptr_new(TIMESTEP),units:''},$
                {name:textoidl('\lambda_{De}'),val:ptr_new(L_d),units:'cm'},$
                {name:textoidl('N_x^{1/3}'),val:ptr_new((d.Ni_x)^(-(1.0/3.0))),units:'cm'},$
                {name:textoidl('{L_{i,inertial}, L_{e,inertial}}'),$
                 val:ptr_new([L_i_inrt,L_e_inrt]),units:'cm'}, $
                {name:textoidl('{\rho_{e}, \rho_{i}, \rho_s}'),$
                 val:ptr_new([rho_e,rho_i,rho_s]),units:'cm'}, $
                {name:textoidl('L'),val:ptr_new(round(L)),units:'cm'}, $
                {name:textoidl('L_n'),$
                 val:ptr_new(round(L_n)),units:'cm'},$
                {name:textoidl('L_{te}'),$
                 val:ptr_new(round(L_te)),units:'cm'},$
                {name:textoidl('L_{ti}'),$
                 val:ptr_new(round(L_ti)),units:'cm'},$
                {name:textoidl('normalization'),val:ptr_new(norm_s),units:''}]
  


  N_names = N_ELEMENTS(var_values)
  
  row_size = 1.0/(N_names+1)

  ;plot,[0,1],[0,1],psym=2
 

  !p.multi = [0,1,1]
  set_plot,'ps'
  device,filename=save_name,/color,YSIZE=25,YOFFSET =1,landscape =0
  
  
  !p.charsize = 1

  for i =0,N_names-1 do begin 

     ;if mean(double(*(var_values.val)[i])) GE 10 then sci = 1 else sci = 0
     sci = 1
     val_str = sigfig(*(var_values.val)[i],2,sci = sci)
     ;print,val_str
     

     if (var_values.units)[i] NE '' then begin
        ;; xyouts,.1,(i+1)*row_size,comment(strjoin([(var_values.name)[i],"[",(var_values.units)[i],"]"]),sigfig(*(var_values.val)[i],2,sci=sci)),CHARSIZE = 1,/NORMAL
        xyouts,.1,(i+1)*row_size,(var_values.name)[i] + '['+ $
               (var_values.units)[i]+']'+': '+$
               strjoin(sigfig(*(var_values.val)[i],2,sci=1),'  '),CHARSIZE = 1,/NORMAL
        
     endif else begin
        
        ;if size(*(var_values.val)[i],/type) LE 5 then $
     
          
        xyouts,.1,(i+1)*row_size,(var_values.name)[i] + ': '+strjoin(sigfig(*(var_values.val)[i],2,sci=1),'  '),CHARSIZE = 1,/NORMAL

        ;print,(var_values.name)[i]+ val_str
        ;;xyouts,.1,(i+1)*row_size,(var_values.name)[i]+
        ;;val_str,CHARSIZE = 1,/NORMAL
        ;xyouts,.1,(i+1)*row_size,val_str,CHARSIZE = 1,/NORMAL
               
     endelse
  endfor


                                ;look again at the bout.inp file to
                                ;see which variables were actually
                                ;evolved

  ;see which variable have to be pulled up

  help,evolve
  print,N_ELEMENTS(evolve)
  print,strcmp(evolve[0],"Ajpar")
  print,evolve[0],"Ajpar"
  
 
  ;pull up simulation results for all evolved fields
if KEYWORD_SET(tind) then begin
   for i = 0,N_ELEMENTS(evolve) -1 do begin
      void = execute(strcompress(string(evolve[i]),/remove_all)+ '_d = temporary(collect(var = "'+strcompress(evolve[i],/remove_all)+'",tind = tind,path= path,/debug))')
   endfor
endif else begin   
   for i = 0,N_ELEMENTS(evolve) -1 do begin
      void = execute(strcompress(string(evolve[i]),/remove_all)+ '_d = temporary(collect(var = "'+strcompress(evolve[i],/remove_all)+'",path= path))')
   endfor
endelse

if KEYWORD_SET(profiles) then begin
   ni0 = temporary(collect(var = "Ni0",path= path))
   Te0 = temporary(collect(var = "Te0",path= path))
   Ti0 = temporary(collect(var = "Ti0",path= path))
endif 

help,Ni0

if KEYWORD_SET(monitor) then begin
   vEB = temporary(collect(var = "vEB",path= path))
   vGradP =  temporary(collect(var = "vGradP",path= path))
   
   vEB = temporary(basic_info(vEB))
   vGradP = temporary(basic_info(vGradP))

   F_Ni = temporary(collect(var = "F_Ni",path= path))
   F_Ajpar = temporary(collect(var = "F_Ajpar",path= path))
   F_Te = temporary(collect(var = "F_Te",path= path))

   F_Ni = temporary(basic_info(F_Ni))
   F_Ajpar = temporary(basic_info(F_Ajpar))
   F_Te = temporary(basic_info(F_Te))
   
   
endif
help,Ni_d

print,"printing to : ",save_name
!p.multi = [0,1,2]
!x.thick=1
!y.thick=1
!p.thick=1
!p.charthick=1
!p.charsize =1


nx = round(d.nx)
ny = round(d.ny)
plot,d.Rxy,d.Zxy,psym = 4,ystyle = 2,xstyle = 2,TITLE = "R vs Z cross section",$
     xtitle = "radius [cm]",ytitle = "Z [cm]"
oplot,d.Rxy[0:MXG-1,*],d.Zxy[0:MXG-1,*],psym = 1,color = 75
oplot,d.Rxy[nx-MXG:nx-1,*],d.Zxy[nx-MXG:nx-1,*],psym = 1,color = 75 

for i = 0,ny-1 do begin
   oplot, d.Rxy[*,i],d.Zxy[*,i]
endfor

for i = 0,nx-1 do begin
   oplot, d.Rxy[i,*],d.Zxy[i,*]
endfor

;return,d
; read in all the evolved variables
for i = 0,N_ELEMENTS(evolve) -1 do begin
   void = execute(strcompress(evolve[i],/remove_all)+ '_str = temporary(basic_info('+strcompress(evolve[i],/remove_all)+'_d,time_x=time_x,info= ptr_new(d),L=L_str,/rotate))')
   
   help,strcompress(evolve[i],/remove_all)+'_str'
endfor

ni_str = Ni_str

time = indgen(round(ni_str.nt))

;; nx = round(ni_str.nx)
;; ny = round(ni_str.ny)


;spawn,"date +%M%S",date
;post12save=strjoin(["post_bout_",strtrim(string(runid),1),".ps"])
;postboutsave=strjoin(["post_bout.ps"])
;surface,(total(d.hthe,1,/cumulative)),R,Z
;plot,(total(d.hthe,1,/cumulative))



;vertex_list = [reform(d.Rxy,nx*ny)nx-MXG-1:nx-1,reform(d.Zxy,nx*ny))]


  
loadct,40
!x.thick=2
!y.thick=2
!p.thick=2
!p.charthick=2
!p.charsize =1
!p.multi = [0,1,3]




if KEYWORD_SET(profiles) then begin
   plot,d.Rxy[*,d.ny/2]*1e2,d.Ni0[*,d.ny/2],$
        TITLE = "Ni0",xtitle = "radius [cm]",ytitle = textoidl('cm^{-3}'),$
        ystyle = 2
   plot,d.Rxy[*,d.ny/2]*1e2,d.Ti0[*,d.ny/2],TITLE = "Ti0 [eV]",xtitle = "radius [cm]",$
        ytitle = 'eV',ystyle = 2
   
   plot,d.Rxy[*,d.ny/2]*1e2,d.Te0[*,d.ny/2],TITLE = "Te0 [eV]",xtitle = "radius [cm]",$
        ytitle = 'eV',ystyle = 2

   ;L_te
   plot,d.Rxy[*,d.ny/2]*1e2,abs(d.Te0[*,d.ny/2]/deriv(d.Rxy[*,d.ny/2]*1e2,d.Te0[*,d.ny/2])),$
        TITLE = textoidl(' \lambda_{Te}'),xtitle = "radius [cm]",$
        ytitle = 'cm',ystyle = 2,/YLOG
   l_indx = where(abs(deriv(d.Te0[*,d.ny/2])) EQ max(abs(deriv(d.Te0[*,d.ny/2]))))
   lam_te = (abs(d.Te0[*,d.ny/2]/deriv(d.Rxy[*,d.ny/2]*1e2,d.Te0[*,d.ny/2])))
   
   oplot,d.Rxy[l_indx,d.ny/2]*1e2,lam_te[l_indx],psym = 1

   xyouts,d.Rxy[l_indx,d.ny/2]*1e2,2*lam_te[l_indx],$
          sigfig(lam_te[l_indx],2)+' cm',charsize = 1

   ;L_n
   plot,d.Rxy[*,d.ny/2]*1e2,abs(d.Ni0[*,d.ny/2]/deriv(d.Rxy[*,d.ny/2]*1e2,d.Ni0[*,d.ny/2])),$
        TITLE = textoidl(' \lambda_{N}'),xtitle = "radius [cm]",$
        ytitle = 'cm',ystyle = 2,/YLOG
   l_indx = where(abs(deriv(d.Ni0[*,d.ny/2])) EQ max(abs(deriv(d.Ni0[*,d.ny/2]))))
   lam_ni = (abs(d.Ni0[*,d.ny/2]/deriv(d.Rxy[*,d.ny/2]*1e2,d.Ni0[*,d.ny/2])))
   
   oplot,d.Rxy[l_indx,d.ny/2]*1e2,lam_ni[l_indx],psym = 1

   xyouts,d.Rxy[l_indx,d.ny/2]*1e2,2*lam_ni[l_indx],$
          sigfig(lam_ni[l_indx],2)+' cm',charsize = 1
endif

if KEYWORD_SET(monitor) then begin 
plot,vEB.ave.amp_gt, TITLE = 'vEB',/YLOG
plot,vGradP.ave.amp_gt, TITLE = 'vGradP',/YLOG

plot,F_Ni.ave.amp_gt,TITLE = "F_Ni"
plot,F_Ajpar.ave.amp_gt,TITLE = "F_Ajpar"
plot,F_Te.ave.amp_gt,TITLE = "F_Te"

endif


;looping over evolved fields for the purpose of producing visuals

for i =0,N_ELEMENTS(evolve) -1 do begin
   void = execute('current_str = ptr_new('+ strcompress(string(evolve[i]),/remove_all)+ '_str)')
   help,*current_str
   
   ;; loadct,0
;;    shade_surf,((*current_str).ave.data)[round(nx/2),*,*,0]
   loadct,40
   !p.multi = [0,2,2]
   ;polar_contour,transpose(reform(((*current_str).ave.data)[*,d.ny/2,*,0])),(TWOPI/nz)*indgen(nz-1,/double),d.Rxy[*,d.ny/2],nlevels = 20,/cell_fill,/isotropic
   
   
   surface,dist(10),/save,/nodata,xrange = [-max(d.Rxy),max(d.Rxy)],$
           yrange = [-max(d.Rxy),max(d.Rxy)],$
           zrange = [min(d.Zxy),max(d.Zxy)]

   axis,yaxis = 1, 1.0, 0.0, 0.0, /T3D
   axis,xaxis = 1, 0.0, 1.0, 0.0, /T3D

   for ii = 0,d.nx-1 do begin
      plots,d.Rxy[ii,*]*cos(d.zeta[ii,*]),d.Rxy[ii,*]*sin(d.zeta[ii,*]),d.Zxy[ii,*],linestyle = 0, /T3D,color = ii*30 + 100
      plots,d.Rxy[ii,*]*cos(d.zeta[ii,*]),d.Rxy[ii,*]*sin(d.zeta[ii,*]),d.Zxy[ii,*],psym = 2, /T3D
   endfor

   ;; for ii = 0,d.ny-1 do begin
   ;;    plots,d.Rxy[*,ii]*cos(d.zeta[*,ii]),d.Rxy[*,ii]*sin(d.zeta[*,ii]),d.Zxy[*,ii],linestyle = 2, /T3D
   ;;    plots,d.Rxy[*,ii]*cos(d.zeta[*,ii]),d.Rxy[*,ii]*sin(d.zeta[*,ii]),d.Zxy[*,ii],psym = 2, /T3D
   ;; endfor


   surface,dist(10),/save,/nodata,xrange = [-max(d.Rxy),max(d.Rxy)],$
           yrange = [-max(d.Rxy),max(d.Rxy)],$
           zrange = [min(d.Zxy),max(d.Zxy)]

   axis,yaxis = 1, 1.0, 0.0, 0.0, /T3D
   axis,xaxis = 1, 0.0, 1.0, 0.0, /T3D

   dqint =  TWOPI*zmax*(1.0/nz)

   for jj = 0,nz/5-1 do begin ; loop over z slices
      for ii = 0,0 do begin
         ;;plots,d.Rxy[ii,*]*cos(d.qinty[ii,*]+ 5*dqint*jj),d.Rxy[ii,*]*sin(d.qinty[ii,*]+5*dqint*jj),d.Zxy[ii,*],linestyle = 0, /T3D,color = ii*30 + 100
         ;;plots,d.Rxy[ii,*]*cos(d.qinty[ii,*]+ 5*dqint*jj),d.Rxy[ii,*]*sin(d.qinty[ii,*]+5*dqint*jj),d.Zxy[ii,*],psym=1, /T3D,color = ii*30 + 10
         ;;plots,d.Rxy[ii,*]*cos(d.qinty[ii,*]+dqint),d.Rxy[ii,*]*sin(d.qinty[ii,*]+dqint*jj),d.Zxy[ii,*],psym = 1, /T3D,color = ii*30 + 100 
         
      endfor
   endfor
   
   ;surface,dist(10),/save,/nodata,xrange = [-max(d.Rxy),max(d.Rxy)],$
   ;        yrange = [-max(d.Rxy),max(d.Rxy)],$
   ;        zrange = [min(d.Zxy),max(d.Zxy)]
   
   ;axis,yaxis = 1, 1.0, 0.0, 0.0, /T3D
   ;axis,xaxis = 1, 0.0, 1.0, 0.0, /T3D
  
   
   for jj = 0,d.ny-1 do begin     ; along field line
      for ii = 0,0 do begin ;; flux surfaces

         if (jj mod 5) eq 0 then begin
            plots,d.Rxy[ii,jj]*cos(d.qinty[ii,jj]+ dqint*indgen(nz)),d.Rxy[ii,jj]*sin(d.qinty[ii,jj]+dqint*indgen(nz)),d.Zxy[ii,jj],psym=1, /T3D,color = 0
            plots,d.Rxy[ii,jj]*cos(d.qinty[ii,jj]+ dqint*indgen(nz)),d.Rxy[ii,jj]*sin(d.qinty[ii,jj]+dqint*indgen(nz)),d.Zxy[ii,jj],linestyle = 0, /T3D,color = ii*30 + 200,thick  = 3
         endif
         ;connect to the prev point along field line
         
         ;; plots,d.Rxy[ii,*]*cos(d.qinty[ii,*]+dqint),d.Rxy[ii,*]*sin(d.qinty[ii,*]+dqint*jj),d.Zxy[ii,*],psym = 1, /T3D,color = ii*30 + 100 
         
      endfor
   endfor
   
   for ii = 0,0 do begin
      plots,d.Rxy[0,*]*cos(d.qinty[0,*]+ dqint*ii),d.Rxy[0,*]*sin(d.qinty[0,*]+ dqint*ii),d.Zxy[0,*],linestyle = 0, /T3D,color = 10,thick = 3
   endfor
   
   ;; plot,d.Rxy[0,*]*cos(d.zeta[0,*]),d.Rxy[0,*]*sin(d.zeta[0,*]),psym = 2,/isotropic,color = 0*180 + 100
   ;; for ii = 0,d.nx-1 do begin
   ;;    oplot,d.Rxy[ii,*]*cos(d.zeta[ii,*]),d.Rxy[ii,*]*sin(d.zeta[ii,*] ),psym = 2,color = ii*180 + 100
   ;; endfor

   a = mg_cc_demo()
   clipboard = Obj_New("IDLgrClipboard", Dimensions=[4,3], Units=1, $
                       Resolution=[2.54/600., 2.54/600.])
 
  ;clipboard->Draw, oView, Filename='objgfx.eps', /PostScript
   clipboard->Draw, a,Filename='objgfx.eps', /PostScript
   print,strcompress(strjoin(['psmerge ','-o',grid,'.ps ',grid,'.ps',' objgfx.eps']))
   spawn,strcompress(strjoin(['psmerge ','-o','output.ps ',grid,'.ps',' objgfx.eps']))
   spawn,strcompress(strjoin(['psmerge ','-o','output2.ps ','output.ps',' objgfx.eps']))
   ;spawn,strcompress(strjoin(['psmerge ','-o','output.ps ',grid,'.ps',' objgfx.eps']))
   print,strcompress(strjoin(['cp ','output.ps ',grid,'.ps']))
   ;spawn,strcompress(strjoin(['cp ','output.ps ',grid,'.ps']))

  ;spawn,strjoin(["ps2pdf ","objgfx.eps", " ","objgfx.pdf"])
   
   ;; plot,d.Rxy[0,*],psym = 2,color = 0*30 + 100
   ;; for ii = 0,d.nx-1 do begin
   ;;    oplot,d.Rxy[ii,*],psym = 2,color = ii*30 + 100
   ;; endfor
   
   ;plot,d.ycoord[*,10],psym = 2
   
   ;for ii = 0,d.nx-1 do begin
   ;   oplot,d.ycoord[ii,*],psym = 2,color = ii*30 + 100
   ;endfor
   
  
   ;; surface,dist(10),/save,/nodata,xrange = [-max(d.Rxy),max(d.Rxy)],$
   ;;         yrange = [-max(d.Rxy),max(d.Rxy)],$
   ;;         zrange = [min(d.Zxy),max(d.Zxy)]

   ;; axis,yaxis = 1, 1.0, 0.0, 0.0, /T3D
   ;; axis,xaxis = 1, 0.0, 1.0, 0.0, /T3D
   
   ;; plots,d.Rxy*cos(d.zeta),d.Rxy*sin(d.zeta),d.Zxy,linestyle = 2, /T3D,color = ii*30 + 100
   ;; plots,d.Rxy*cos(d.zeta),d.Rxy*sin(d.zeta),d.Zxy,psym = 2, /T3D,color = 74
   
   ;polar_contour,transpose(reform(((*current_str).ave.data)[*,d.ny/2,*,nt-1])),(TWOPI/nz)*indgen(nz-1,/double),d.Rxy[*,d.ny/2],nlevels = 20,/cell_fill,/isotropic
  ;;  plot,d.Rxy[*,d.ny/2]
   
;;    contour,transpose(reform(((*current_str).ave.data)[*,d.ny/2,*,0])),/fill,nlevels = 30

   kz_i = (TWOPI/L_z)*(indgen((*current_str).fft.nz+1))*rho_i
   kz_e = kz_i*(rho_e/rho_i)
   kz_s = kz_i*(rho_s/rho_i)

   knorm = (TWOPI/lbNorm)*(indgen((*current_str).fft.nz+1))*rho_s

   ky = (TWOPI/lpar)*indgen(round(ni_str.fft.ny)+1)*rho_s
   ;;later to be determed from BOUT.inp or the source file
   kz = kz_s
   plot,kz,title ='kz'
   kz_max = N_ELEMENTS(kz)-1
   if N_ELEMENTS(kz) EQ 32 then begin
      kz_max = 25
      kz = kz[0:kz_max]
   endif
   
   if N_ELEMENTS(time_x) GT N_ELEMENTS((*current_str).ave.amp_gt) then begin
      time_x = time_x[0: N_ELEMENTS((*current_str).ave.amp_gt)-1]
      status = textoidl('crash at t_{'+ N_ELEMENTS((*current_str).ave.amp_gt)+'}')
      nt = N_ELEMENTS(time_x)-1
   endif else begin
      status = 'complete'
   endelse
   
   

   gamma_est = abs(deriv(time_x,(*current_str).ave.amp_gt)/(*current_str).ave.amp_gt)
   
   if total(finite(gamma_est)) LT 3 then gamma_est =(*current_str).ave.amp_gt * 0

   plot,time_x,(*current_str).ave.amp_gt,TITLE = strcompress(evolve[i],/remove_all)+'.amp_gt',/YLOG,ystyle = 1,XTITLE = textoidl('[T \omega_{cs}]')
   plot,time_x,gamma_est,$
        YTITLE = textoidl('\gamma /\omega_{cs}'),$
        XTITLE = textoidl('[T \omega_{cs}]'),/YLOG,ystyle =8,xstyle= 9
   axis,/xaxis,xrange = double([0,nt*dt])/w_n
   axis,/yaxis,yrange = [min(gamma_est),max(gamma_est)] * w_n,$
        ytitle = textoidl('\gamma [1/s]')
   

   ;plot,abs(ni_str.DC),TITLE= 'Ni.DC',/YLOG
   
   
   nlevels = N_ELEMENTS(ky)
   
   contour,alog((*current_str).fft.pow.kyt),ky,time_x,/fill,NLEVELS = nlevels,$
           TITLE = strcompress(evolve[i],/remove_all)+'_ky',$
           YTITLE = textoidl('[T \omega_{cs}]'),$
           XTITLE = textoidl('k_{y}')

   contour,(*current_str).fft.pow.kyt,ky,time_x,/fill,NLEVELS = nlevels,$
           TITLE = strcompress(evolve[i],/remove_all)+'_ky',$
           YTITLE = textoidl('[T \omega_{cs}]'),$
           XTITLE = textoidl('k_{y}')
   
   contour,alog((*current_str).fft.pow.kyt),ky,time_x,/overplot
   ;nt = N_ELEMENTS(time_x)-1

   plot,ky,((*current_str).fft.pow.kyt)[*,nt],/YLOG,TITLE =textoidl("K_{\parallel} spectrum at the first and last time step")
   arrow,ky[1],((*current_str).fft.pow.kyt)[1,0],ky[1],((*current_str).fft.pow.kyt)[1,nt],/data
   
   oplot,ky,((*current_str).fft.pow.kyt)[*,0],color = 75

   ;plot the phase of the dominant ky mode as a function of x and z

   ;plot,(*current_str).fft*

   temp = max((*current_str).fft.pow.ky,DIMENSION = 1)
   ;plot,temp
   
   contour,alog(((*current_str).fft.pow.kzt)[0:kz_max,*]),kz,time_x,/fill,NLEVELS = 50,$
           TITLE = strcompress(evolve[i],/remove_all)+'_kz',$
           YTITLE = textoidl('[T \omega_{cs}]'),$
           XTITLE = textoidl('k_{\zeta} \rho_s'),$
           ystyle =8,xstyle= 9
   
   axis,/xaxis,xrange = [.0,N_ELEMENTS(kz)-1]
   axis,/yaxis,yrange = [min(time_x),max(time_x)]/ w_n,$
        ytitle = textoidl('[s]')
   
   ;radial power 

   
   
   contour,((*current_str).fft.pow.kzt)[0:kz_max,*],kz,time_x,/fill,NLEVELS = 50,$
           TITLE = strcompress(evolve[i],/remove_all)+ textoidl('k_{\zeta}'),$
           YTITLE = textoidl('[T \omega_{cs}]'),$
           XTITLE = textoidl('k_{\zeta} \rho_s'),$
           ystyle =8,xstyle= 9
   
   axis,/xaxis,xrange = [.0,N_ELEMENTS(kz)-1]
   axis,/yaxis,yrange = [min(time_x),max(time_x)]/ w_n,$
        ytitle = textoidl('[s]')
   

 ;;   contour,(*current_str).fft.pow.kz_full,/fill,NLEVELS = 50,$
;;            TITLE = strcompress(evolve[i],/remove_all)+'_kz_full',$
;;            YTITLE = textoidl('[T \omega_{cs}]'),$
;;            XTITLE = textoidl('k_{z} \rho')

;;    contour,alog((*current_strs).pow.kz),kz,time_x,/overplot

   plot,kz,((*current_str).fft.pow.kzt)[0:kz_max,nt],/YLOG,TITLE =textoidl('K_{\zeta} spectrum at the first and  last time step'),XTITLE = textoidl('k_{\zeta} \rho_s'),YRANGE =[min((*current_str).fft.pow.kzt),max((*current_str).fft.pow.kzt)],ystyle = 2
   arrow,kz[1],((*current_str).fft.pow.kzt)[1,0],kz[1],((*current_str).fft.pow.kzt)[1,nt],/data

   arrow,kz[10],((*current_str).fft.pow.kzt)[10,0],kz[10],((*current_str).fft.pow.kzt)[10,nt],/data
   

   oplot,kz,((*current_str).fft.pow.kzt)[0:kz_max,0],color = 75, linestyle = 2,thick = 3
    
   plot,knorm,((*current_str).fft.pow.kzt)[0:kz_max,nt],/YLOG,TITLE =textoidl('K_{\norm (sloppy)} '),XTITLE = textoidl('k_{\zeta} \rho_s'),YRANGE =[min((*current_str).fft.pow.kzt),max((*current_str).fft.pow.kzt)],ystyle = 2
   arrow,knorm[1],((*current_str).fft.pow.kzt)[1,0],knorm[1],((*current_str).fft.pow.kzt)[1,nt],/data

   arrow,knorm[10],((*current_str).fft.pow.kzt)[10,0],knorm[10],((*current_str).fft.pow.kzt)[10,nt],/data
   

   oplot,kz,((*current_str).fft.pow.kzt)[0:kz_max,0],color = 75, linestyle = 2,thick = 3
   ;plot the phase of the dominant kz mode as a function of x and y

   
   if keyword_set(rotate) then begin
      ;; plot,rho_s*(*(((*current_str).rot.kzz)[0]))[0:nz/2]
      ;; plot,kz[0:nz/2-1]
      
   contour,alog(((*current_str).rot_fft.pow.kzt)[0:nz,*]),rho_s*(*(((*current_str).rot.kzz)[0]))[0:nz],time_x,/fill,NLEVELS = 50,$
           TITLE = strcompress(evolve[i],/remove_all)+'_kz',$
           YTITLE = textoidl('[T \omega_{cs}]'),$
           XTITLE = textoidl('k_{\zeta} \rho_s'),$
           ystyle =8,xstyle= 9

   contour,alog(((*current_str).rot_fft.pow.kyt)[0:ny/2,*]),rho_s*(*(((*current_str).rot.kyy)[0]))[0:ny/2],time_x,/fill,NLEVELS = 50,$
           TITLE = strcompress(evolve[i],/remove_all)+'_kz',$
           YTITLE = textoidl('[T \omega_{cs}]'),$
           XTITLE = textoidl('k_{b \times \nabla r} \rho_s'),$
           ystyle =8,xstyle= 9
    
   knorm = rho_s*(*(((*current_str).rot.kzz)[0]))[0:nz]
   
   plot,knorm,((*current_str).rot_fft.pow.kzt)[0:nz,nt],/YLOG,TITLE =textoidl('K_{\norm} '),XTITLE = textoidl('k_{\zeta} \rho_s'),YRANGE =[min((*current_str).rot_fft.pow.kzt),max((*current_str).rot_fft.pow.kzt)],ystyle = 2
   arrow,knorm[1],((*current_str).rot_fft.pow.kzt)[1,0],knorm[1],((*current_str).rot_fft.pow.kzt)[1,nt],/data
   
   arrow,knorm[10],((*current_str).rot_fft.pow.kzt)[10,0],knorm[10],((*current_str).rot_fft.pow.kzt)[10,nt],/data
   

   oplot,knorm,((*current_str).rot_fft.pow.kzt)[0:nz,0],color = 75, linestyle = 2,thick = 3
   
;; axis,/xaxis,xrange = [.0,N_ELEMENTS(kz)-1]
   ;; axis,/yaxis,yrange = [min(time_x),max(time_x)]/ w_n,$
   ;;      ytitle = textoidl('[s]')
endif
endfor

  device,/close
  ;;postboutsave_pdf=strjoin(["post_bout.pdf"])
  ;; postboutsave_pdf = " post_bout.pdf"
  ;;help,postboutsave,postboutsave_pdf
  ;spawn,strjoin(["ps2pdf ",postboutsave,postboutsave_pdf])
  print, strcompress(strjoin([save_dir,grid,sigfig(dt,2),'.pdf']),/remove_all)
  print,sigfig(dt,2),save_dir,grid
  pdf_save_name = strcompress(strjoin([grid,sigfig(dt,2),'.pdf']),/remove_all)
  spawn,strjoin(["ps2pdf ",save_name, " ",pdf_save_name])
  ;spawn,strjoin(["rm ",save_name])
  
  !p.multi = [0,1,1]
  !x.thick=1
  !y.thick=1
  !p.thick=1
  !p.charthick=1
  !p.charsize =1

  set_plot,'X'
  ;return,d
  ;;gamma = {ni:1.0/dt * mean((deriv(ni_str.amp))[0:10]), phi:1.0/dt * mean((deriv(phi_str.amp))[0:10]),$
           ;;rho:1.0/dt * mean((deriv(rho_str.amp))[0:10]),Te:1.0/dt * mean((deriv(Te_str.amp))[0:10])}
  ;;return, ptr_new(gamma)
  return,0
end




