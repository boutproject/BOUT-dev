module global_data_g

  implicit none
  INTEGER*4 i,j, idum
  INTEGER*4 iunit, ios
  INTEGER*4 mrfac, nx4, ny4, nlim, nbdry, kyord, nyefit, nxefit 
  INTEGER*4 kxord, nworkf
  INTEGER*4 nbndry

  DOUBLE PRECISION xdim, zdim, rcentr, rgrid1, zmid
  DOUBLE PRECISION rmagx, zmagx,simagx,sibdry,bcentr
  DOUBLE PRECISION cpasma, xdum

  DOUBLE PRECISION, dimension(:), allocatable :: fpol, pres, &
       workk1, workk2, qpsi, rbdry, zbdry, xlim, ylim
  DOUBLE PRECISION, dimension(:,:), allocatable :: fold

  character*80 runid
  character*8 label(6)

end module global_data_g



module global_data_a

  implicit none

  character header*80
  
  character uday*10
  INTEGER*4 vmonth,vday,vyear
  INTEGER*4 eshot,ktime1
  DOUBLE PRECISION etime
  INTEGER*4 jflag,lflag,limloc,mco2v,mco2r,qmflag
  
  DOUBLE PRECISION tsaisq,rcencm,bcentr,pasmat
  DOUBLE PRECISION cpasma,rout,zout,aout
  DOUBLE PRECISION eout,doutu,doutl,vout
  DOUBLE PRECISION rcurrt,zcurrt,qsta,betat
  DOUBLE PRECISION betap,ali,oleft,oright
  DOUBLE PRECISION otop,obott,qpsi95,vertn

  DOUBLE PRECISION, dimension(:), allocatable :: rco2v
  DOUBLE PRECISION, dimension(:), allocatable :: dco2v
  DOUBLE PRECISION, dimension(:), allocatable :: rco2r
  DOUBLE PRECISION, dimension(:), allocatable :: dco2r
  
  DOUBLE PRECISION shearb,bpolav,s1,s2
  DOUBLE PRECISION s3,qout,olefs,orighs
  DOUBLE PRECISION otops,sibdry,areao,wplasm
  DOUBLE PRECISION terror,elongm,qqmagx,cdflux
  DOUBLE PRECISION alpha,rttt,psiref,xndnt
  DOUBLE PRECISION rseps1,zseps1,rseps2,zseps2

  DOUBLE PRECISION sepexp,obots,btaxp,btaxv
  DOUBLE PRECISION aaq1,aaq2,aaq3,seplim
  DOUBLE PRECISION rmagx,zmagx,simagx,taumhd
  
  DOUBLE PRECISION betapd,betatd,wplasmd,fluxx
  DOUBLE PRECISION vloopt,taudia,qmerci,tavem
  
  INTEGER*4 nsilop, magpri, nfcoil, nesum
  
  DOUBLE PRECISION, dimension(:), allocatable :: csilop
  DOUBLE PRECISION, dimension(:), allocatable :: cmpr2
  DOUBLE PRECISION, dimension(:), allocatable :: ccbrsp
  DOUBLE PRECISION, dimension(:), allocatable :: eccurt
  
  DOUBLE PRECISION pbinj,rvsin,zvsin,rvsout
  DOUBLE PRECISION zvsout,vsurfa,wpdot,wbdot
  DOUBLE PRECISION slantu,slantl,zuperts,chipre
  DOUBLE PRECISION cjor95,pp95,ssep,yyy2
  DOUBLE PRECISION xnnc,cprof,oring,cjor0
  
  DOUBLE PRECISION fexpan,qqmin,chigamt,ssi01
  DOUBLE PRECISION fexpvs,sepnose,ssi95,rqqmin
  DOUBLE PRECISION cjor99,cj1ave,rmidin,rmidout
  DOUBLE PRECISION psurfa,xdum

end module global_data_a


subroutine readg 


  use global_data_g
  INTEGER*4 stat


  !!c *************** read the output of the EFIT code ******************

  print *, 'Reading neqdsk file .......'


  iunit=1
  open (iunit, file='neqdsk', form='formatted', iostat=ios, status='old')

  if (ios .ne. 0) then
     print *, "**** neqdsk file not found, exiting..."
     STOP
  else

     read(iunit,2000) (label(i),i=1,6),idum,nxefit,nyefit
     print *, "   nxefit=", nxefit, ", nyefit=", nyefit

     runid=label(1)//label(2)//label(3)//label(4)//label(5)//label(6)
     print *, "   runid=", runid


     read(iunit,2020) xdim,zdim,rcentr,rgrid1,zmid
     read(iunit,2020) rmagx,zmagx,simagx,sibdry,bcentr
     read(iunit,2020) cpasma,simagx,xdum,rmagx,xdum
     read(iunit,2020) zmagx,xdum,sibdry,xdum,xdum
  
  
     !!-read arrays

     !!-do it in case of a repetitive call to this function
     deallocate (fpol, stat=stat)
     deallocate (pres, stat=stat)
     deallocate (workk1, stat=stat)
     deallocate (workk2, stat=stat)
     deallocate (fold, stat=stat)
     deallocate (qpsi, stat=stat)
     deallocate (rbdry, stat=stat)
     deallocate (zbdry, stat=stat)
     deallocate (xlim, stat=stat)
     deallocate (ylim, stat=stat)


     allocate (fpol(1:nxefit))
     read(iunit,2020) (fpol(i),i=1,nxefit)
     
     allocate (pres(1:nxefit))
     read(iunit,2020) (pres(i),i=1,nxefit)
     
     allocate (workk1(1:nxefit))
     read(iunit,2020) (workk1(i),i=1,nxefit)

     allocate (workk2(1:nxefit))     
     read(iunit,2020) (workk2(i),i=1,nxefit)
     
     allocate (fold(1:nxefit,1:nyefit))
     read(iunit,2020) ((fold(i,j),i=1,nxefit),j=1,nyefit)
     
     allocate (qpsi(1:nxefit))
     read(iunit,2020) (qpsi(i),i=1,nxefit)
     
     read(iunit,2022) nbdry,nlim
     allocate (rbdry(1:nbdry))
     allocate (zbdry(1:nbdry))
     read(iunit,2020) (rbdry(i),zbdry(i),i=1,nbdry)

     allocate (xlim(1:nlim))
     allocate (ylim(1:nlim))
     read(iunit,2020) (xlim(i),ylim(i),i=1,nlim)

     close (iunit)

2000 format(6a8,3i4)
2020 format(5e16.9)
2022 format(2i5)
     
  endif
  
  print *, '.......finished'
end subroutine readg



subroutine writeg 

  use global_data_g
  INTEGER*4 stat

  !!c *************** write the output of the EFIT code ******************

  print *, 'Writing neqdsk2 file .......'


  iunit=1
  open (iunit, file='neqdsk2', form='formatted', iostat=ios, status='unknown')

  if (ios .ne. 0) then
     print *, "**** cannot open neqdsk2"
     STOP
  endif

  write(iunit,2000) (label(i),i=1,6),idum,nxefit,nyefit
  print *, "   nxefit=", nxefit, ", nyefit=", nyefit

  runid=label(1)//label(2)//label(3)//label(4)//label(5)//label(6)
  print *, "   runid=", runid


  write(iunit,2020) xdim,zdim,rcentr,rgrid1,zmid
  write(iunit,2020) rmagx,zmagx,simagx,sibdry,bcentr
  write(iunit,2020) cpasma,simagx,xdum,rmagx,xdum
  write(iunit,2020) zmagx,xdum,sibdry,xdum,xdum



  !!-write arrays

  write(iunit,2020) (fpol(i),i=1,nxefit)
  write(iunit,2020) (pres(i),i=1,nxefit)
  write(iunit,2020) (workk1(i),i=1,nxefit)
  write(iunit,2020) (workk2(i),i=1,nxefit)
  write(iunit,2020) ((fold(i,j),i=1,nxefit),j=1,nyefit)
  write(iunit,2020) (qpsi(i),i=1,nxefit)
  write(iunit,2022) nbdry,nlim
  write(iunit,2020) (rbdry(i),zbdry(i),i=1,nbdry)
  write(iunit,2020) (xlim(i),ylim(i),i=1,nlim)

  close (iunit)

2000 format(6a8,3i4)
2020 format(5e16.9)
2022 format(2i5)

  
  deallocate (fpol, stat=stat)
  deallocate (pres, stat=stat)
  deallocate (workk1, stat=stat)
  deallocate (workk2, stat=stat)
  deallocate (fold, stat=stat)
  deallocate (qpsi, stat=stat)
  deallocate (rbdry, stat=stat)
  deallocate (zbdry, stat=stat)
  deallocate (xlim, stat=stat)
  deallocate (ylim, stat=stat)
  

  print *, '.......finished'
end subroutine writeg


subroutine reada 


  use global_data_a
  INTEGER*4 stat



  !!c *************** read the output of the EFIT code ******************

  print *, 'Reading aeqdsk file .......'


  iunit=1
  open (iunit, file='aeqdsk', form='formatted', iostat=ios, status='old')

  if (ios .ne. 0) then
     print *, "**** aeqdsk file not found, exiting..."
     STOP
  else

     read (iunit,1056) uday,vmonth,vday,vyear
     print *, 'uday=', uday

     read (iunit,*) eshot,ktime1
     read (iunit,1040) etime
     read (iunit,1060) etime,jflag,lflag,limloc,mco2v,mco2r,qmflag

     read (iunit,1040) tsaisq,rcencm,bcentr,pasmat
     read (iunit,1040) cpasma,rout,zout,aout
     read (iunit,1040) eout,doutu,doutl,vout
     read (iunit,1040) rcurrt,zcurrt,qsta,betat
     read (iunit,1040) betap,ali,oleft,oright
     read (iunit,1040) otop,obott,qpsi95,vertn

     deallocate (rco2v, stat=stat)
     deallocate (dco2v, stat=stat)
     deallocate (rco2r, stat=stat)
     deallocate (dco2r, stat=stat)

     allocate (rco2v(1:mco2v))
     allocate (dco2v(1:mco2v))
     allocate (rco2r(1:mco2r))
     allocate (dco2r(1:mco2r))

     read (iunit,1040) (rco2v(k),k=1,mco2v)
     read (iunit,1040) (dco2v(k),k=1,mco2v)
     read (iunit,1040) (rco2r(k),k=1,mco2r)
     read (iunit,1040) (dco2r(k),k=1,mco2r)

     read (iunit,1040) shearb,bpolav,s1,s2
     read (iunit,1040) s3,qout,olefs,orighs
     read (iunit,1040) otops,sibdry,areao,wplasm
     read (iunit,1040) terror,elongm,qqmagx,cdflux
     read (iunit,1040) alpha,rttt,psiref,xndnt
     read (iunit,1040) rseps1,zseps1,rseps2,zseps2
     print *, 'rseps1,zseps1,rseps2,zseps2=', rseps1,zseps1,rseps2,zseps2

     read (iunit,1040) sepexp,obots,btaxp,btaxv
     read (iunit,1040) aaq1,aaq2,aaq3,seplim
     read (iunit,1040) rmagx,zmagx,simagx,taumhd

     read (iunit,1040) betapd,betatd,wplasmd,fluxx
     read (iunit,1040) vloopt,taudia,qmerci,tavem


     !-old efit files don't have this line (???)
     read (iunit,1044) nsilop,magpri,nfcoil,nesum

     deallocate (csilop, stat=stat)
     allocate (csilop(1:nsilop))

     deallocate (cmpr2, stat=stat)
     allocate (cmpr2(1:magpri))

     deallocate (ccbrsp, stat=stat)
     allocate (ccbrsp(1:nfcoil))

     deallocate (eccurt, stat=stat)
     allocate (eccurt(1:nesum))

     read (iunit,1040) (csilop(k),k=1,nsilop),(cmpr2(k),k=1,magpri)
     read (iunit,1040) (ccbrsp(k),k=1,nfcoil)
     read (iunit,1040) (eccurt(k),k=1,nesum)

     read (iunit,1040) pbinj,rvsin,zvsin,rvsout
     read (iunit,1040) zvsout,vsurfa,wpdot,wbdot
     read (iunit,1040) slantu,slantl,zuperts,chipre
     read (iunit,1040) cjor95,pp95,ssep,yyy2
     read (iunit,1040) xnnc,cprof,oring,cjor0

     read (iunit,1040) fexpan,qqmin,chigamt,ssi01
     read (iunit,1040) fexpvs,sepnose,ssi95,rqqmin
     read (iunit,1040) cjor99,cj1ave,rmidin,rmidout
     read (iunit,1040) psurfa,xdum,xdum,xdum
 
     read (iunit,1042) header
     print *, 'header=', header

     close(iunit)

1040 format (1x,4e16.9)
1042 format (1x,a80)
1044 format (1x,4i5)
1050 format (1x,i5,11x,i5)
1053 format (1x,i6,11x,i5)
1055 format (1x,a10,2a5)
1056 format (1x,a10,i2,1x,i2,1x,i4)
!!1056 format (1x,a10,i2,"/",i2,"/",i4)
1057 format (1x,i6.6,14x,i2)
1060 format (1x,f7.2,10x,i5,11x,i5,1x,a3,1x,i3,1x,i3,1x,a3)
!!1060 format (1h*,f7.2,10x,i5,11x,i5,1x,a3,1x,i3,1x,i3,1x,a3)
     
  endif
  
  print *, '.......finished'
end subroutine reada


subroutine writea 


  use global_data_a
  INTEGER*4 stat



  !!c *************** read the output of the EFIT code ******************

  print *, 'Writing aeqdsk2 file .......'


  iunit=1
  open (iunit, file='aeqdsk2', form='formatted', iostat=ios, status='replace')

  if (ios .ne. 0) then
     print *, "**** Cannot open aeqdsk2"
     STOP
  else


     write (iunit,1056) uday,vmonth,vday,vyear
     write (iunit,1057) eshot,ktime1
     write (iunit,1040) etime
     write (iunit,1060) etime,jflag,lflag,limloc,mco2v,mco2r,qmflag

     write (iunit,1040) tsaisq,rcencm,bcentr,pasmat
     write (iunit,1040) cpasma,rout,zout,aout
     write (iunit,1040) eout,doutu,doutl,vout
     write (iunit,1040) rcurrt,zcurrt,qsta,betat
     write (iunit,1040) betap,ali,oleft,oright
     write (iunit,1040) otop,obott,qpsi95,vertn

     write (iunit,1040) (rco2v(k),k=1,mco2v)
     write (iunit,1040) (dco2v(k),k=1,mco2v)
     write (iunit,1040) (rco2r(k),k=1,mco2r)
     write (iunit,1040) (dco2r(k),k=1,mco2r)

     write (iunit,1040) shearb,bpolav,s1,s2
     write (iunit,1040) s3,qout,olefs,orighs
     write (iunit,1040) otops,sibdry,areao,wplasm
     write (iunit,1040) terror,elongm,qqmagx,cdflux
     write (iunit,1040) alpha,rttt,psiref,xndnt
     write (iunit,1040) rseps1,zseps1,rseps2,zseps2

     write (iunit,1040) sepexp,obots,btaxp,btaxv
     write (iunit,1040) aaq1,aaq2,aaq3,seplim
     write (iunit,1040) rmagx,zmagx,simagx,taumhd

     write (iunit,1040) betapd,betatd,wplasmd,fluxx
     write (iunit,1040) vloopt,taudia,qmerci,tavem


     !-old efit files don't have this line (???)
     write (iunit,1044) nsilop,magpri,nfcoil,nesum




     write (iunit,1040) (csilop(k),k=1,nsilop),(cmpr2(k),k=1,magpri)
     write (iunit,1040) (ccbrsp(k),k=1,nfcoil)
     write (iunit,1040) (eccurt(k),k=1,nesum)

     write (iunit,1040) pbinj,rvsin,zvsin,rvsout
     write (iunit,1040) zvsout,vsurfa,wpdot,wbdot
     write (iunit,1040) slantu,slantl,zuperts,chipre
     write (iunit,1040) cjor95,pp95,ssep,yyy2
     write (iunit,1040) xnnc,cprof,oring,cjor0

     write (iunit,1040) fexpan,qqmin,chigamt,ssi01
     write (iunit,1040) fexpvs,sepnose,ssi95,rqqmin
     write (iunit,1040) cjor99,cj1ave,rmidin,rmidout
     write (iunit,1040) psurfa,xdum,xdum,xdum

     write (iunit,1042) header

     close(iunit)

1040 format (1x,4e16.9)
1042 format (1x,a80)
1044 format (1x,4i5)
1050 format (1x,i5,11x,i5)
1053 format (1x,i6,11x,i5)
1055 format (1x,a10,2a5)
1056 format (1x,a10,i2.2,"/",i2.2,"/",i4)
1057 format (1x,i6.6,14x,i2)
1060 format (1h*,f7.2,10x,i5,11x,i5,1x,a3,1x,i3,1x,i3,1x,a3)


     deallocate (csilop, stat=stat)
     deallocate (cmpr2, stat=stat)
     deallocate (ccbrsp, stat=stat)
     deallocate (eccurt, stat=stat)
     deallocate (rco2v, stat=stat)
     deallocate (dco2v, stat=stat)
     deallocate (rco2r, stat=stat)
     deallocate (dco2r, stat=stat)

  endif
  
  print *, '.......finished'
end subroutine writea


subroutine get_nxy (nxefit_, nyefit_)
!
! return dimensions so that memory can be allocated at main (i.e. IDL) level
!----------------------------------------------------------------------------!

  use global_data_g
  INTEGER*4, intent (OUT) :: nxefit_, nyefit_

  !-can I just return-by-value???
  print *, "***getting g-dims"
  nxefit_=nxefit
  nyefit_=nyefit

end subroutine get_nxy



subroutine get_psi (&
     nxefit_, nyefit_, fold_, xdim_, zdim_,&
     rcentr_, rgrid1_, zmid_, rmagx_, zmagx_,&
     simagx_, sibdry_, bcentr_)
!
! return a copy of array fold, the memory is allocated by the caller
!----------------------------------------------------------------------------!

  use global_data_g

  INTEGER*4, intent (IN) :: nxefit_, nyefit_
  DOUBLE PRECISION, dimension(nxefit_,nyefit_), intent (IN OUT)  :: fold_
  DOUBLE PRECISION, intent (IN OUT)  :: xdim_, zdim_, rcentr_, rgrid1_, zmid_, &
       rmagx_, zmagx_, simagx_, sibdry_, bcentr_

  print *, "***getting g-data"
  fold_=fold !!-copy arrays
  xdim_=xdim 
  zdim_=zdim
  rcentr_=rcentr
  rgrid1_=rgrid1
  zmid_=zmid
  rmagx_=rmagx 
  zmagx_=zmagx 
  simagx_=simagx
  sibdry_=sibdry
  bcentr_=bcentr
  
end subroutine get_psi



subroutine put_psi(&
     nxefit_, nyefit_, fold_, xdim_, zdim_,&
     rcentr_, rgrid1_, zmid_, rmagx_, zmagx_,&
     simagx_, sibdry_, bcentr_)
!
! modify the values of array fold stored as global data
!----------------------------------------------------------------------------!

  use global_data_g

  INTEGER*4, intent (IN) :: nxefit_, nyefit_
  DOUBLE PRECISION, dimension(nxefit_,nyefit_), intent (IN)  :: fold_
  DOUBLE PRECISION, intent (IN) :: xdim_, zdim_, rcentr_, rgrid1_, zmid_, &
       rmagx_, zmagx_, simagx_, sibdry_, bcentr_

  print *, "***putting g-data"
  fold=fold_ !!-copy arrays
  xdim=xdim_
  zdim=zdim_
  rcentr=rcentr_
  rgrid1=rgrid1_
  zmid=zmid_
  rmagx=rmagx_
  zmagx=zmagx_
  simagx=simagx_
  sibdry=sibdry_
  bcentr=bcentr_

end subroutine put_psi



subroutine get_rzseps (rseps_1, zseps_1, rseps_2, zseps_2) 
!
! get rseps, zseps values stored in the aeqdsk file
!----------------------------------------------------------------------------!

  use global_data_a
  DOUBLE PRECISION, intent (INOUT) :: rseps_1, zseps_1, rseps_2, zseps_2
  
  rseps_1=rseps1
  zseps_1=zseps1
  rseps_2=rseps2
  zseps_2=zseps2
  print *, '***getting a-data, ', rseps_1, zseps_1, rseps_2, zseps_2

end subroutine get_rzseps 



subroutine put_rzseps (rseps_1, zseps_1, rseps_2, zseps_2) 
!
! get rseps, zseps values stored in the aeqdsk file
!----------------------------------------------------------------------------!

  use global_data_a
  DOUBLE PRECISION, intent (IN) :: rseps_1, zseps_1, rseps_2, zseps_2
  
  print *, '***putting a-data, ', rseps_1, zseps_1, rseps_2, zseps_2
  rseps1 = rseps_1
  zseps1 = zseps_1
  rseps2 = rseps_2
  zseps2 = zseps_2

end subroutine put_rzseps 



subroutine writedata (&
     nxm,nym,ixpt1,ixpt2,iysptrx1,&
     rm,zm,psi,br,bz,bpol,bphi,b)


  implicit none

  character*80 fname !!="gridue"
  character*80 runidg !!="snow-flake-plus"

  INTEGER*4 nxm, nym, ixpt1,ixpt2,iysptrx1
  double precision, dimension(nxm+2,nym+2,5), intent (IN) :: rm,zm,psi,br,bz,bpol,bphi,b

  INTEGER*4 nunit, ix, iy, n
  INTEGER*4 iunit, ios

  fname="gridue"
  runidg="snow-flake-plus"


  iunit=1
  open (iunit, file='gridue', form='formatted', iostat=ios, status='replace')

  if (ios .ne. 0) then
     print *, "**** Cannot open gridue"
     STOP
  else

     write(iunit,1999) nxm,nym,ixpt1,ixpt2,iysptrx1
     write(iunit,2000)
     write(iunit,2001) (((rm(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     write(iunit,2000)
     write(iunit,2001) (((zm(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     write(iunit,2000)
     write(iunit,2001) (((psi(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     write(iunit,2000)
     write(iunit,2001) (((br(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     write(iunit,2000)
     write(iunit,2001) (((bz(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     write(iunit,2000)
     write(iunit,2001) (((bpol(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     write(iunit,2000)
     write(iunit,2001) (((bphi(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     write(iunit,2000)
     write(iunit,2001) (((b(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     write(iunit,2002) runidg

     close (iunit)

     write(*,*) 'Wrote file "', fname, '" with runidg:  ', runidg
     write(*,*)
          
  endif
     


1999 format(5i4)  
2000 format()

     !!  ifelse([WORDSIZE],64,\
     !!2001 format(1p3e23.15)
     !!  ,\
     !!2001 format(1p3d23.15)
     !!  )\
     
2001 format(1p3d23.15)
2002 format(a60)
     
end subroutine writedata
