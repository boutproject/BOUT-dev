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

        if (iu eq 0) then cu=0.707107D else cu=1.D
        if (jv eq 0) then cv=0.707107D else cv=1.D

          sum=0.0D
 
          for jy=0,Ny-1 do begin
           for ix=0,Nx-1 do begin
           ;     
            sum=sum + sig[ix,jy]*cos(jv*!DPI*(2*jy+1)/(2*Ny))*$
                                  cos(iu*!DPI*(2*ix+1)/(2*Nx))
           ;
           endfor
         endfor
        
      fsig[iu,jv]=(2.D/nsig)*cv*cu*sum

   endfor
  endfor

ENDIF ELSE BEGIN ;---inverse transform---

  sig=dblarr(nx,ny)

  for ix=0,Nx-1 do begin
   for jy=0,Ny-1 do begin

        sum=0.0D

        for iu=0,Nx-1 do begin
           for jv=0,Ny-1 do begin
           ;     
           if (iu eq 0) then cu=0.707107D else cu=1.D
           if (jv eq 0) then cv=0.707107D else cv=1.D

            sum=sum + cv*cu*fsig[iu,jv]*cos(jv*!DPI*(2*jy+1)/(2*Ny))*$
                                         cos(iu*!DPI*(2*ix+1)/(2*Nx))

           ;
           endfor
        endfor

       sig[ix,jy]=(2.D/nsig)*sum
   
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

        sum=0.0D
        sumx=0.0D
        sumy=0.0D

        for iu=0,Nx-1 do begin
           for jv=0,Ny-1 do begin
           ;     
           if (iu eq 0) then cu=0.707107D else cu=1.D
           if (jv eq 0) then cv=0.707107D else cv=1.D

            sum=sum + cv*cu*fsig[iu,jv]*$
              COS(jv*!DPI*(2*y0+1)/(2*Ny))*COS(iu*!DPI*(2*x0+1)/(2*Nx))

            sumx=sumx + cv*cu*fsig[iu,jv]*$
               COS(jv*!DPI*(2*y0+1)/(2*Ny))*SIN(iu*!DPI*(2*x0+1)/(2*Nx))*$
                (-iu*!DPI/Nx)

            sumy=sumy + cv*cu*fsig[iu,jv]*$
              SIN(jv*!DPI*(2*y0+1)/(2*Ny))*COS(iu*!DPI*(2*x0+1)/(2*Nx))*$
               (-jv*!DPI/Ny)

           ;
           endfor
        endfor

   res=(2.D/nsig)*[sum,sumx,sumy]   

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

  cuvec=dblarr(nsig)+1.D
   cuvec[0]=0.707107D

  cvvec=dblarr(nsig)+1.D
   cvvec[0]=0.707107D
 

     uvec=findgen(nsig)
     uvec=COS(!DPI*findgen(nsig)*(x0+0.5D)/nsig)
     uvex=(-findgen(nsig)*!DPI/nsig)*SIN(!DPI*findgen(nsig)*(x0+0.5D)/nsig)

     vvec=findgen(nsig)
     vvec=COS(!DPI*vvec*(y0+0.5D)/nsig)
     vvey=(-findgen(nsig)*!DPI/nsig)*SIN(!DPI*findgen(nsig)*(y0+0.5D)/nsig)


       ;-value
       res=(2.D/nsig) * TOTAL(((cuvec # cvvec) * fsig) * (uvec # vvec))  

       ;d/dx
       rex=(2.D/nsig) * TOTAL(((cuvec # cvvec) * fsig) * (uvex # vvec))  

       ;d/dy
       rey=(2.D/nsig) * TOTAL(((cuvec # cvvec) * fsig) * (uvec # vvey))  

;
;
;
return, [res,rex,rey]
end
