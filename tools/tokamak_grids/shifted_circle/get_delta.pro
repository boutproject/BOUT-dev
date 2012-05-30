; New version of shifted circle equilibrium
;=========================================================;

function str2str, s, debug=debug
;
;-Convert structure to a string suitable for printing with xyouts
;-----------------------------------------------------------------;

 ntags=n_tags(s)
 names=tag_names(s)

  str=''

  for i=0,ntags-1 do begin
   str=str+ names[i]+' = '+ STRTRIM(STRING(s.(i)),2) + '!c'
   if keyword_set(DEBUG) then STOP
  endfor  

return,str
end


function powlaw, x, a, b, c, deriv=deriv

 if keyword_set(DERIV) then begin
  res=b*c*x^(c-1.)
 endif else begin
  res=a + b*x^c
 endelse 

return, res
end


function tanhlaw, x, a, b, c, deriv=deriv

 if keyword_set(DERIV) then begin
  res=0.5*((c-1)/b)*(1./cosh((x-a)/b))^2
 endif else begin
  res=0.5*(1-tanh((x-a)/b))*(1-c)+c
 endelse 

return, res
end


function explaw, x, a, b, c, deriv=deriv

 if keyword_set(DERIV) then begin
     res=-((1-c)/b)*exp(-(x-a)/b)
 endif else begin
     res=exp(-(x-a)/b)*(1-c)+c
 endelse 

return, res
end


function npres, rho, deriv=deriv, cgs=cgs
;
; normalized pressure profile as a function of rho=r/a
;----------------------------------------------------------;

 ;gamma=4.0
 ;nu=-0.15

 COMMON cmn, p

 case p.pfun of 
     "exp": res=ExpLaw(rho,p.rho0,p.nu,p.gamma, deriv=deriv)
     "tanh": res=TanhLaw(rho,p.rho0,p.nu,p.gamma, deriv=deriv)
     else: STOP, 'unknown value of pfun'
 endcase


 if keyword_set(CGS) then begin
   ;-return actual plasma pressure in CGS units
   p0=(p.betap*p.B0^2*p.eps^2/(8*!PI))
   res=res*p0
 endif

return, res
end


function qprof, rho
;
; q profile as a function of rho=r/a
;---------------------------------------------------------;

  ;q0=1.5
  ;q1=1.0
  ;delta=2.0

 COMMON cmn, p
  q0=p.q0
  q1=p.q1
  delta=p.delta  

  res=PowLaw(rho,q0,q1,delta)

return, res
end

function integrand1, x
;
; integrand involved in calculation of Delta'
;-------------------------------------------;

 COMMON cmn, p
 res=x*(p.betap*x*npres(x,/D) - (x/qprof(x))^2)

return, res
end


function integrand2, rho
;
; integrand involved in calculation of Delta
;-------------------------------------------;

  COMMON cmn, p

  res=(p.eps*(qprof(rho))^2/(rho+1e-9)^3)*QROMB( 'integrand1', 0, rho)

return, res
end


function integrand3, r
;
; integrand involved in calculation of polodial flux
;---------------------------------------------------;

  COMMON cmn, p

  res=r/qprof(r/p.amin)

return, res
end


function differential, r, theta
;
; RHS of ODE involved in calculation of theta(r,\tilde(\theta))
;--------------------------------------------------------------;

  COMMON cmn, p

  res=integrand2(r/p.amin)*sin(theta)/r 

return,res
end


function get_theta, r, tnew, debug=debug, ORDERED=ORDERED, OLDORIG=OLDORIG
;
; Calculate theta(rho) for theta_new=const
;
; If set keyword OLDORIG theta_new goes counter-clockwise [0,2Pi] starting 
; at the outer midplane; ELSE ttheta_new goes clockwise from -Pi to Pi, 
; where the point -P<=>Pi corresponds to inner midplane location; otherwise
;--------------------------------------------------------------------------;

  COMMON cmn, p

  if NOT keyword_set(OLDORIG) then begin
    tnew=!Pi-tnew
  endif 

  ;-add phase shift
  ;tnew=tnew+!pi/4

  rinit=p.amin ;-initial r value
  tinit=[-tnew]   ;-initial theta value

  nr=n_elements(r)
  res=fltarr(nr)

  if keyword_set(DEBUG) then STOP

  IF keyword_set(ORDERED) then begin
  ;-this is for an ordered vector of r

    res[nr-1]=tnew

    for i=nr-2,0,-1 do begin
     rstep=r[i]-r[i+1]
     rinit=r[i+1]
     tinit=[-[res[i+1]]]      
     res[i]=-LSODE(tinit,rinit,rstep,'differential') 
   endfor

  ENDIF ELSE BEGIN
   ;-for arbitrary vector of r - but slow!

    for i=0,nr-1 do begin
     rstep=r[i]-rinit
     res[i]=-LSODE(tinit,rinit,rstep,'differential') 
    endfor

 ENDELSE

;;stop
return,res
end


function get_RZ, delta, rho, theta_new, thetag=thetag, debug=debug
;
; Calculate R,Z coords for given delta, rho, theta_new
; delta and rho can be vectors, theta_new is a scalar
;-----------------------------------------------------;

  COMMON cmn, p

  ;-calculate the [vector of] geometric theta
  thetag=get_theta(p.amin*rho,theta_new,/ordered)

  Rmaj=p.R0+delta + (p.amin*rho)*cos(thetag)
  Z=(p.amin*rho)*sin(thetag)

if keyword_set(DEBUG) then STOP
return,{R:Rmaj,Z:Z}
end



function get_B, delta, rho, thetag
;
; Calculate components of B field for given delta, rho, theta_new
; delta, rho, and thetag can be vectors
;-----------------------------------------------------;

  COMMON cmn, p

  Rmaj=p.R0+delta+p.amin*rho*cos(thetag)
  Z=p.amin*rho*sin(thetag)

  bpol=p.B0*((p.amin*rho)/(Rmaj*qprof(rho)))/(1.+ (integrand2(rho))*cos(thetag))
  btor=p.B0*p.R0/Rmaj
  br=bpol*sin(thetag)
  bz=-bpol*cos(thetag)
  btot=sqrt(bpol^2+btor^2)

return,{bp:bpol,bt:btor,br:br,bz:bz,btot:btot}
end


pro get_delta, rarr=rarr, deltaparr=deltaparr, deltarr=deltarr, psiarr=psiarr, rhoarr=rhoarr, $
set=set, show=show, noreset=noreset, qarr=qarr
;
;-calculate shafranov shift Delta(rho)
;-------------------------------------------;

  Set_params, set=set, noreset=noreset
  COMMON cmn, p

  if keyword_set(rhoarr) then begin
    nn=n_elements(rhoarr)
    rarr=rhoarr*p.amin
  endif else begin
    nn=101
    rarr=p.amin*findgen(nn)/(nn-1) ;-array of r values
  endelse

  deltarr=fltarr(nn)
  deltaparr=fltarr(nn)
  psiarr=fltarr(nn)

  for i=0,nn-1 do begin
    ;
    deltaparr[i]=integrand2(rarr[i]/p.amin) 
    ;
    deltarr[i]=p.Delta0 + p.amin*QROMB( 'integrand2', 0, rarr[i]/p.amin) 
    ;
    psiarr[i]=p.B0*QROMB('integrand3',0,rarr[i])
    ;
  endfor
  
  qarr = qprof(rarr/p.amin)
  

  if keyword_set(SHOW) then begin

    print,'Showing contours'
    !p.multi=[0,2,2,0,0]


    ;-plot a few r=const lines
    IF (1) then begin
      plot,/nod, [p.R0-p.Delta0-p.amin,p.R0+p.Delta0+p.amin],[-1,1]*p.amin,/iso,$
        tit='With theta=const lines', xtit='R [cm]', ytit='Z [cm]'
      for i=0,nn-1,10 do begin
        ntheta=101
        theta=2*!PI*findgen(ntheta)/(ntheta-1)
        ;rnow=p.amin*i/(nn-1)
        rnow=rarr[i]
        rline_r=p.R0+deltarr[i]+rnow*cos(theta)
        rline_z=rnow*sin(theta)
        oplot, rline_r,rline_z
      endfor
    ENDIF

    ;-overplot a few theta=const lines
    IF (1) then begin
      ntheta=21
      for i=0,ntheta-1 do begin
        ;nr=101
        ;rr=p.amin*findgen(nr)/(nr-1)
        ;rr=rarr
        theta_now=2*!pi*i/(ntheta-1)
        ;tline_r=p.R0+delta[FIX((nn-1)*rarr/p.amin)]+rarr*cos(theta_now)
        tline_r=p.R0 + deltarr + rarr*cos(theta_now)
        tline_z=rarr*sin(theta_now)
        oplot, tline_r, tline_z
      endfor
    ENDIF

    oplot, p.R0*[1,1],P.amin*[-2,2], lin=2, thick=2
    oplot, (p.R0+p.Delta0)*[1,1],P.amin*[-2,2], lin=2, thick=2
    oplot, [0,(p.R0+2*p.amin)], [0,0], lin=2, thick=2


    IF (1) then begin
      ;-plot a few r=const lines
      plot,/nod, [p.R0-p.Delta0-p.amin,p.R0+p.Delta0+p.amin],[-1,1]*p.amin,/iso,$
        tit='With theta_new=const lines', xtit='R [cm]', ytit='Z [cm]'
      for i=0,nn-1,10 do begin
        ntheta=101
        theta=2*!PI*findgen(ntheta)/(ntheta-1)
        rnow=rarr[i]            ;;p.amin*i/(nn-1)
        rline_r=p.R0+deltarr[i]+rnow*cos(theta)
        rline_z=rnow*sin(theta)
        oplot, rline_r,rline_z
      endfor
    ENDIF



    IF (1) then begin 
      ;-overplot a few theta_new=const lines
      ntheta=21
      for i=0,ntheta-1 do begin 
        theta_now=2*!pi*i/(ntheta-1)
        th=get_theta(rarr,theta_now,/ord)
        oplot, p.r0+deltarr+rarr*cos(th), rarr*sin(th)
      endfor
    ENDIF


    oplot, p.R0*[1,1],P.amin*[-2,2], lin=2, thick=2
    oplot, (p.R0+p.Delta0)*[1,1],P.amin*[-2,2], lin=2, thick=2
    oplot, [0,(p.R0+2*p.amin)], [0,0], lin=2, thick=2

    plot, rarr,  qarr, thick=2, tit='q(r) and p1(r)', yr=[0,3], xtit='r [cm]'
    oplot, rarr, npres(rarr/p.amin), thick=2


    plot, rarr, deltarr, tit="Delta(r) and Delta'(r)", thick=2, yr=[-1.,10.], xtit='r [cm]'
    oplot, rarr, deltaparr, thick=2

    PRINT, "Shift of outer surface: ", deltarr[N_ELEMENTS(deltarr)-1]

    ;-show all parameters
    str=str2str(p)
    xyouts,10,8,str, col=2

    !p.multi=0
  endif

end
