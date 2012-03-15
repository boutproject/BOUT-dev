function set_structure, s

ntags=n_tags(s)
names=tag_names(s)

 for i=0,ntags-1 do begin
  prompt= names[i]+ '=' + STRING(s.(i))+', ok? [y/n]'
  ans='y'
  read, ans, prompt=prompt

  if (STRLOWCASE(ans) ne 'y') then begin
   val=s.(i)
   read, val, prompt='Setting '+names[i]+ '='
   s.(i)=val
  endif

 endfor

return, s
end

pro set_params, set=set, noreset=noreset, $
                R0=R0,amin=amin,Del0=Del0,B0=B0,betap=betap,$
                qfun=qfun,q0=q0,q1=q1,delta=delta,$
                pfun=pfun,rho0=rho0,gamma=gamma,nu=nu

;-set parameters

 COMMON cmn, p

    IF keyword_set(SET) then begin	
        
        IF NOT KEYWORD_SET(noreset) THEN Set_params              ;-set defaults first
        p=Set_Structure(p)
        p.eps=p.amin/p.R0       ;-aspect ratio
    ENDIF ELSE BEGIN

     if not keyword_set(NORESET) then begin
       ;print, 'resetting to defaults...'

	  if not keyword_set(R0) then R0=2e2      ;-[cm]
	  if not keyword_set(amin) then amin=50.0   ;-[cm]
	  if not keyword_set(Del0) then Delta0=2.5 else Delta0=Del0 ;-[cm]
	  if not keyword_set(B0) then B0=1e4      ;-[Gauss]
	  if not keyword_set(betap) then betap=0.1   ;-poloidal beta
	  eps=amin/R0 ;-aspect ratio

	  ;-power law for q(rho): q0+q1*rho^delta
          if not keyword_set(qfun) then qfun="pow"
	  if not keyword_set(q0) then q0=1.5      
	  if not keyword_set(q1) then q1=1.0
	  if not keyword_set(delta) then delta=2.0   

          ;-tanh law for p1(rho): 0.5*(1-tanh((rho-rho0)/nu))*(1-gamma)+gamma
          ;-exp law for p1(rho): exp(-(rho-rho0)/nu)*(1-gamma)+gamma

          if not keyword_set(pfun) then pfun="tanh"
	  if not keyword_set(rho0) then rho0=0.85   
	  if not keyword_set(gamma) then gamma=0.3   
	  if not keyword_set(nu) then nu=0.1    
  
          p={R0:R0,amin:amin,Delta0:Delta0,B0:B0,eps:eps,betap:betap, $
             qfun:qfun,q0:q0,q1:q1,delta:delta, $
             pfun:pfun,rho0:rho0,gamma:gamma,nu:nu}

     endif else begin
      ;print, 'no resetting ...'
     endelse
   
   ENDELSE

end
