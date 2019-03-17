pro dct2d, sig, fsig, nsig, inverse=inverse
;
;-calculate 2D discrete cosine transform
;for a square array of size nsig by nsig
;----------------------------------------

nx=nsig
ny=nsig

IF NOT KEYWORD_SET(INVERSE) THEN BEGIN ;---direct transform---

  fsig=dblarr(nx,ny)

  for iu=0,Nx-1 do begin
   for jv=0,Ny-1 do begin

        if (iu eq 0) then cu=0.707107 else cu=1.
        if (jv eq 0) then cv=0.707107 else cv=1.

          sum=0.0
 
          for jy=0,Ny-1 do begin
           for ix=0,Nx-1 do begin
           ;     
            sum=sum + sig[ix,jy]*cos(jv*!PI*(2*jy+1)/(2*Ny))*$
                                  cos(iu*!PI*(2*ix+1)/(2*Nx))
           ;
           endfor
         endfor
        
      fsig[iu,jv]=(2./nsig)*cv*cu*sum

   endfor
  endfor

ENDIF ELSE BEGIN ;---inverse transform---

  sig=dblarr(nx,ny)

  for ix=0,Nx-1 do begin
   for jy=0,Ny-1 do begin

        sum=0.0

        for iu=0,Nx-1 do begin
           for jv=0,Ny-1 do begin
           ;     
           if (iu eq 0) then cu=0.707107 else cu=1.
           if (jv eq 0) then cv=0.707107 else cv=1.

            sum=sum + cv*cu*fsig[iu,jv]*cos(jv*!PI*(2*jy+1)/(2*Ny))*$
                                         cos(iu*!PI*(2*ix+1)/(2*Nx))

           ;
           endfor
        endfor

       sig[ix,jy]=(2./nsig)*sum
   
     ;STOP
   endfor
  endfor


ENDELSE
 

;
;
;
end


function EvalCosP, fsig, nsig, x0=x0,y0=y0
;
;-evaluate the cos polinomial at point x0,y0
;and the partial derivatives
;--------------------------------------------

nx=nsig
ny=nsig

        sum=0.0
        sumx=0.0
        sumy=0.0

        for iu=0,Nx-1 do begin
           for jv=0,Ny-1 do begin
           ;     
           if (iu eq 0) then cu=0.707107 else cu=1.
           if (jv eq 0) then cv=0.707107 else cv=1.

            sum=sum + cv*cu*fsig[iu,jv]*$
              COS(jv*!PI*(2*y0+1)/(2*Ny))*COS(iu*!PI*(2*x0+1)/(2*Nx))

            sumx=sumx + cv*cu*fsig[iu,jv]*$
               COS(jv*!PI*(2*y0+1)/(2*Ny))*SIN(iu*!PI*(2*x0+1)/(2*Nx))*$
                (-iu*!PI/Nx)

            sumy=sumy + cv*cu*fsig[iu,jv]*$
              SIN(jv*!PI*(2*y0+1)/(2*Ny))*COS(iu*!PI*(2*x0+1)/(2*Nx))*$
               (-jv*!PI/Ny)

           ;
           endfor
        endfor

   res=(2./nsig)*[sum,sumx,sumy]   

;
;
;
return, res
end


function EvalCosPfast, fsig, nsig, x0=x0,y0=y0
;
;-evaluate the cos polinomial at point x0,y0
;and the partial derivatives
;--------------------------------------------

  cuvec=dblarr(nsig)+1.
   cuvec[0]=0.707107

  cvvec=dblarr(nsig)+1.
   cvvec[0]=0.707107
 

     uvec=findgen(nsig)
     uvec=COS(!PI*findgen(nsig)*(x0+0.5)/nsig)
     uvex=(-findgen(nsig)*!PI/nsig)*SIN(!PI*findgen(nsig)*(x0+0.5)/nsig)

     vvec=findgen(nsig)
     vvec=COS(!PI*vvec*(y0+0.5)/nsig)
     vvey=(-findgen(nsig)*!PI/nsig)*SIN(!PI*findgen(nsig)*(y0+0.5)/nsig)


       ;-value
       res=(2./nsig) * TOTAL(((cuvec # cvvec) * fsig) * (uvec # vvec))  

       ;d/dx
       rex=(2./nsig) * TOTAL(((cuvec # cvvec) * fsig) * (uvex # vvec))  

       ;d/dy
       rey=(2./nsig) * TOTAL(((cuvec # cvvec) * fsig) * (uvec # vvey))  

;
;
;
return, [res,rex,rey]
end
