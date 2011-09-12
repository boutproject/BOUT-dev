;user can change these inputs here
nstart=50
nmap=300
tend=75
method='bicube_fft'
mc_method='bicubic'
nthreads=!CPU.HW_NCPU-1
;end user inputs

CD,'./data/'
CD,current=path
nst=STRTRIM(STRING(nstart),2)
nm=STRTRIM(STRING(nmap),2)
mt_obj=objarr(nthreads)
ibegin=tend-nthreads+1
print,'Begin initialization loop'
FOR i=ibegin,tend DO BEGIN $
  index=fmodulo(i,nthreads) & $
  str=STRTRIM(STRING(i),2) & $
  cmd="poincare,path=path,tindex="+str+",nstart="+nst+",nmap="+nm+ $
      ",method=meth,mc_method=mc_meth,quiet=2" & $
  print,cmd & $
  mt_obj[index]=OBJ_NEW('IDL_IDLBridge') & $
  mt_obj[index]->execute,"@"+PREF_GET('IDL_STARTUP') & $
  mt_obj[index]->SetVar,"path",path & $
  mt_obj[index]->SetVar,"meth",method & $
  mt_obj[index]->SetVar,"mc_meth",mc_method & $
  mt_obj[index]->execute, cmd, /nowait & $
ENDFOR

print,'End initialization loop'
sum=nthreads
WHILE sum GT 0 DO BEGIN $
  sum=0 & $
  FOR index=0,nthreads-1 DO BEGIN & $
    sum=sum+mt_obj[index]->status() & $
  ENDFOR & $
  print,'Not finished, waiting for 10 minutes' & $
  WAIT,600 & $
ENDWHILE

print,'all done!'
EXIT
