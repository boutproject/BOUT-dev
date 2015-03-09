#include <bout.hxx>
#include <boutmain.hxx>
#include <bout/constants.hxx>
#include <bout/sys/timer.hxx>
#include <cmath>
#include <vector>
#include <utils.hxx>
#include <fft.hxx>
#include <dcomplex.hxx>

Field3D fld1,fld1Bak,fldShift,fldShift2d;
Field3D fft_1, fft_2, fft_2T, fft_3, fft_3T;

Field3D shift1d(const Field3D fin, const BoutReal zangle);
Field3D shift2d(const Field3D fin, const BoutReal zangle, bool trans);
Field3D shift3d(const Field3D fin, const BoutReal zangle, bool trans);
BoutReal fracErr(const Field3D ans, const Field3D val);
int nrep;

int physics_init(bool restarting){
  /////////////////
  // OPTIONS
  /////////////////

  //Get options root, create field etc.
  Options *globalOptions = Options::getRoot();
  OPTION(globalOptions, nrep, 100);
  mesh->geometry();
  solver->add(fld1,"fld1");
  fld1Bak=fld1;
  fft_1=fft_2=fft_2T=fft_3=fft_3T=fldShift=fldShift2d=0.;
  BoutReal zangle=0.5;
  
  //Now do  shiftZ different ways
  for (int irep=0;irep<nrep;irep++){
    Timer timer("2d");
    fft_2=shift2d(fld1,zangle,false);
  }

  for (int irep=0;irep<nrep;irep++){
    Timer timer("2dT");
    fft_2T=shift2d(fld1,zangle,true);
  }

  for (int irep=0;irep<nrep;irep++){
    Timer timer("1d");
    fft_1=shift1d(fld1,zangle);
  }

  for (int irep=0;irep<nrep;irep++){
    Timer timer("3d");
    fft_3=shift3d(fld1,zangle,false);
  }

  for (int irep=0;irep<nrep;irep++){
    Timer timer("3dT");
    fft_3T=shift3d(fld1,zangle,true);
  }

  for (int irep=0;irep<nrep;irep++){
    Timer timer("sz2d");
    fldShift2d=fld1.shiftZ2D(zangle);
  }

  for (int irep=0;irep<nrep;irep++){
    Timer timer("sz2b");
    fldShift2d=fld1;
    fldShift2d.shiftZ2D(zangle,true);
  }

  for (int irep=0;irep<nrep;irep++){
    Timer timer("sz");
    fldShift=fld1.shiftZ(zangle);
  }


  BoutReal time1d=Timer::resetTime("1d");
  BoutReal time2d=Timer::resetTime("2d");
  BoutReal time2dT=Timer::resetTime("2dT");
  BoutReal time3d=Timer::resetTime("3d");
  BoutReal time3dT=Timer::resetTime("3dT");
  BoutReal timesz2d=Timer::resetTime("sz2d");
  BoutReal timesz2b=Timer::resetTime("sz2b");
  BoutReal timesz=Timer::resetTime("sz");
  output<<"######################################"<<endl;
  output<<"Average times:"<<endl;
  output<<"\tsz     : "<<timesz/nrep<<endl;
  output<<"\t1d     : "<<time1d/nrep<<endl;
  output<<"\tMax error% : "<<fracErr(fldShift,fft_1)<<endl;
  output<<"\t2d     : "<<time2d/nrep<<endl;
  output<<"\tMax error% : "<<fracErr(fldShift,fft_2)<<endl;
  output<<"\t2dT    : "<<time2dT/nrep<<endl;
  output<<"\tMax error% : "<<fracErr(fldShift,fft_2T)<<endl;
  output<<"\t3d  : "<<time3d/nrep<<endl;
  output<<"\tMax error% : "<<fracErr(fldShift,fft_3)<<endl;
  output<<"\t3dT : "<<time3dT/nrep<<endl;
  output<<"\tMax error% : "<<fracErr(fldShift,fft_3T)<<endl;
  output<<"\tsz2d b : "<<timesz2b/nrep<<endl;
  output<<"\tsz2d   : "<<timesz2d/nrep<<endl;
  output<<"\tMax error% : "<<fracErr(fldShift,fldShift2d)<<endl;
  output<<"######################################"<<endl;
  
  //Save result
  SAVE_ONCE(fft_1);
  SAVE_ONCE(fft_2); 
  SAVE_ONCE(fft_2T);
  SAVE_ONCE(fft_3);
  SAVE_ONCE(fft_3T);
  SAVE_ONCE(fldShift);
  SAVE_ONCE(fldShift2d); 

  return 1;
};

BoutReal fracErr(const Field3D ans, const Field3D val){
  return 100*max(abs((ans-val)/ans),true);
};

//Currently this routine causes a crash, looks like memory corruption
Field3D shift3d(const Field3D fin, const BoutReal zangle, bool trans){
  const int ncz = mesh->ngz-1;
  const int nx = mesh->ngx;
  const int ny = mesh->ngy;
  const int ntot=nx*ny;
  int nkz = ncz/2+1;
  dcomplex **v2d;
  BoutReal **b2d;
  v2d = cmatrix(ntot,nkz);
  b2d = rmatrix(ntot,ncz);

  Field3D result=fin;

  //First pack array
  for (int jx=0;jx<nx;jx++){
    for (int jy=0;jy<ny;jy++){
      for (int jz=0;jz<ncz;jz++){
	result.getData(jx,jy,jz,&b2d[jx+jy*nx][jz]);
      };
    };
  };

  //Now R->C
  rfft(b2d,ntot,ncz,v2d,trans);

  //Now apply phase
  BoutReal fac=2.0*PI/mesh->zlength;
  for (int jx=0;jx<nx;jx++){
    for (int jy=0;jy<ny;jy++){
      for (int jz=0;jz<ncz;jz++){
  	BoutReal kwave=jz*fac;
  	dcomplex phase(cos(kwave*zangle),-sin(kwave*zangle));
  	v2d[jx+nx*jy][jz] *= phase;
      };
    };
  };

  //C->R
  irfft(v2d,ntot,ncz,b2d,trans);
  
  //Pack result
  for (int jx=0;jx<nx;jx++){
    for (int jy=0;jy<ny;jy++){
      for (int jz=0;jz<ncz;jz++){
	result.setData(jx,jy,jz,&b2d[jx+nx*jy][jz]);
      };
      result.setData(jx,jy,ncz,&b2d[jx+nx*jy][0]);
    };
  };
  free_cmatrix(v2d);
  free_rmatrix(b2d);
  return result;
};

Field3D shift2d(const Field3D fin, const BoutReal zangle, bool trans){
  const int ncz = mesh->ngz-1;
  const int nx = mesh->ngx;
  int nkz = ncz/2+1;
  dcomplex **v2d;
  BoutReal **b2d;
  v2d = cmatrix(nx,nkz);
  b2d = rmatrix(nx,ncz);

  Field3D result=fin;
  for (int jy=0;jy<mesh->ngy;jy++){
    //FieldPerp tmp=fin.slice(jy);

    //R->C
    for (int jx=0;jx<nx;jx++){
      for (int jz=0;jz<ncz;jz++){
	fin.getData(jx,jy,jz,&b2d[jx][jz]);
      };
    };
    //BoutReal **t=tmp.getData();
    
    rfft(b2d,nx,ncz,v2d,trans);
    //rfft(t,nx,ncz,v2d);
    
    //Phase
    for (int jz=0;jz<nkz;jz++){
      BoutReal kwave=jz*2.0*PI/mesh->zlength;
      dcomplex phase(cos(kwave*zangle),-sin(kwave*zangle));
      for (int jx=0;jx<nx;jx++){
	v2d[jx][jz] *= phase;
      };
    };

    //C->R
    irfft(v2d,nx,ncz,b2d,trans);
    //irfft(v2d,nx,ncz,t);
    
    for (int jx=0;jx<nx;jx++){
      for (int jz=0;jz<ncz;jz++){
	result.setData(jx,jy,jz,&b2d[jx][jz]);
	//result.setData(jx,jy,jz,&t[jx][jz]);
      };
      result.setData(jx,jy,ncz,&b2d[jx][0]);
      //result.setData(jx,jy,ncz,&t[jx][0]);
    };
  };
  free_cmatrix(v2d);
  free_rmatrix(b2d);
  return result;
};

Field3D shift1d(const Field3D fin, const BoutReal zangle){
  const int ncz = mesh->ngz-1;
  const int nx = mesh->ngx;
  int nkz = ncz/2+1;
  dcomplex *v1d;
  v1d = new dcomplex[nkz];

  Field3D result=fin;
  
  BoutReal *t = new BoutReal[ncz];
  for (int jy=0;jy<mesh->ngy;jy++){
    for (int jx=0;jx<mesh->ngx;jx++){
      //R->C
      for (int jz=0;jz<ncz;jz++){
	fin.getData(jx,jy,jz,&t[jz]);
      };
      rfft(t,ncz,v1d);
    
      //Phase
      for (int jz=0;jz<nkz;jz++){
	BoutReal kwave=jz*2.0*PI/mesh->zlength;
	dcomplex phase(cos(kwave*zangle),-sin(kwave*zangle));
	v1d[jz] *= phase;
      };

      //C->R
      irfft(v1d,ncz,t);
      for (int jz=0;jz<ncz;jz++){
	result.setData(jx,jy,jz,&t[jz]);
      };
      result.setData(jx,jy,ncz,&t[0]);
    };
  };
  delete[] t;
  delete[] v1d;
  return result;
};

int physics_run(BoutReal time){
  ddt(fld1)=0.;
  return 0;
};
