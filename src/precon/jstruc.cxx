/*calculate locations of nonzero elelments in the Jacobian*/

#include <stdio.h>
#include <petsc.h>
#include <globals.hxx>
#include <interpolation.hxx> // Cell interpolation

//#define NR_END 1

#define IJKth_index(ivar,ix,jy,kz,NVARS,NSMX,NSMXY) \
 ((ivar-1) + (ix)*NVARS + (jy)*NSMX + (kz)*NSMXY)

//-neighbors on 3D stencil
#if 0
//-1st and 2nd order neighbors
enum neib_type {xp1,xp2,xm1,xm2, yp1,yp2,ym1,ym2, zp1,zp2,zm1,zm2};
int NNEIB=12; //-number of neighbors on stencil
#else
//-1st order neighbors
enum neib_type {xp1,xm1,yp1,ym1,zp1,zm1};
int NNEIB=6; //-number of neighbors on stencil
#endif

#if 1

/*--------Parameters set in BOUT.inp---------*/
//long NVARS=2;  //-number of plasma field variables
//long mesh->NXPE=1;   //-number of radial subdomains
//long MXSUB=4;  //40; //-number of radial polongs per subdomain
//long NYPE=1;   //-number of poloidal subdomains
//long MYSUB=8;  //64; //-number of poloidal polongs per subdomain
//long MZ=5;     //65;    // toroidal grid size +1
//long MYG=2;    //-poloidal guard cells
//long MXG=2;    //-radial guard cells

/*auxiliary parameters*/
//long NSMX, NSMXY, mesh->ngx, mesh->ngy, mesh->ngz, neq, local_N;


void nrerror(string error_text)
{
  msg_stack.push("Run-time error...\n%s\n", error_text.c_str());
  bout_error("...now exiting to system...");
}


int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  int **m;
  int NR_END=1;

  /* allocate pointers to rows */
  m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}



int Map2sv(int ivar, long ixgrid, long iygrid, long izgrid, long NVARS, long NSMX, long NSMY)
/*
  Find 1D index for mapping to the state vector for given:
  iVar-index of fluid variable, (0,1,2,3,4,5) for (rho,te,ti,ni,up,ajpar),
  ixgrid,iygrid,izgrid - indices of grid node
*/
{
  int indx1d;
  //long NVARS=2; // should get from BOUT.inp!!!
  //long NSMX, NSMXY;

  indx1d=IJKth_index(ivar,ixgrid,iygrid,izgrid,NVARS,NSMX,NSMY);

  return indx1d;
}



void Neighbor(enum neib_type neib, long ix1, long iy1, long iz1, long *ix2, long *iy2, long *iz2)
/*
  For given grid point (ix1,iy1,iz1) calculate indices of 
  neighbor grid point (ix2,iy2,iz2)
*/
{

  //-set default
  *ix2=ix1;
  *iy2=iy1;
  *iz2=iz1;


  //-next point in positive x-direction
  if(neib==xp1)
  {
    if (ix1 == mesh->ngx-1)
    {
      *ix2=0; //-open-periodic ([0,1,..mesh->ngx-1]) where 0<=>mesh->ngx
    }
    else
    {
      *ix2=ix1+1;
    }
  }


  //-next point in negative x-direction
  if(neib== xm1)
  {
    if (ix1 == 0)
    {
      *ix2=mesh->ngx-1; //-open-periodic ([0,1,..mesh->ngx-1]) where 0<=>mesh->ngx
    }
    else
    {
      *ix2=ix1-1;
    }
  }


  //-next point in positive y-direction
  if(neib==yp1)
  {
    if (iy1 == mesh->ngy-1)
    {
      *iy2=0; //-open-periodic ([0,1,..mesh->ngy-1]) where 0<=>mesh->ngy
    }
    else
    {
      *iy2=iy1+1;
    }
  }


  //-next point in negative y-direction
  if(neib==ym1)
  {
    if (iy1 == 0)
    {
      *iy2=mesh->ngy-1; //-open-periodic ([0,1,..mesh->ngy-1]) where 0<=>mesh->ngy
    }
    else
    {
      *iy2=iy1-1;
    }
  }


  //-next point in positive z-direction
  if(neib==zp1)
  {
    if (iz1 == mesh->ngz-1)
    {
      *iz2=1; //-closed-periodic ([0,1,..mesh->ngz-1]) where 0<=>mesh->ngz-1
    }
    else
    {
      *iz2=iz1+1;
    }
  }


  //-next point in negative z-direction
  if(neib==zm1)
  {
    if (iz1 == 0)
    {
      *iz2=mesh->ngz-2; //-closed-periodic ([0,1,..mesh->ngz-1]) where 0<=>mesh->ngz-1
    }
    else
    {
      *iz2=iz1-1;
    }
  }



  output << "neib=" << neib << " ";
  output << ">> ix1=" << ix1 << ", iy1=" << iy1 << ", iz1=" << iz1 << "; ";
  output << "<< ix2=" << *ix2 << ", iy2=" << *iy2 << ", iz2=" << *iz2 << endl;
}



int jstruc(int NVARS, int NXPE, int MXSUB, int NYPE, int MYSUB, int MZ, int MYG, int MXG)
{
  int **jmatr;

  int ivar1, ivar2;
  long ix1, iy1, iz1;
  long ix2, iy2, iz2;
  long ij1, ij2;
  enum neib_type neib;
  int neib_tmp;

  //long NVARS=2; // get from BOUT.inp!!!
  long neq,local_N;


  /*calculate auxiliary parametes*/
  mesh->ngx=MXSUB*mesh->NXPE+4; //-guard cells
  mesh->ngy=MYSUB*NYPE;
  mesh->ngz=MZ;
  int ncz=mesh->ngz-1;
  int NSMX = NVARS*mesh->ngx;
  int NSMXY = NVARS*mesh->ngx*MYSUB;
  neq = NVARS*mesh->ngx*mesh->ngy*ncz;
  //local_N = NVARS*MXSUB*MYSUB*mesh->ngz;


  output << endl << "jstruc(): " << endl;
  output << "NVARS=" << NVARS << endl;
  output << "mesh->ngx" << mesh->ngx << endl;
  output << "mesh->ngy" << mesh->ngy << endl;
  output << "mesh->ngz" << mesh->ngz << endl;
  output << "neq=" << neq << endl;

  //-allocate matrix for the jacobian (use list instead)
  output << "Allocating Jacobian matrix, [" << neq << "x" << neq << "]" << endl;
  jmatr=imatrix(0,neq-1,0,neq-1);


  /*====================Local interaction on the stencil====================*/

  for (ivar1=1;ivar1<=NVARS;ivar1++)
  {
    for (ivar2=1;ivar2<=NVARS;ivar2++)
  	{
  	  for (ix1=0;ix1<mesh->ngx;ix1++)
	    {
        for (iy1=0;iy1<MYSUB;iy1++)
    		{
    		  for (iz1=0;iz1<MZ;iz1++) 
  		    {
	     	      
  		      //-find 1D index for 1st variable at this grid point
  		      ij1=Map2sv(ivar1,ix1,iy1,iz1, NVARS,NSMX,NSMXY);

		      
  		      //-loop over the stencil
  		      //for (neib=0;neib<NNEIB;neib++)
            for (neib_tmp=0;neib_tmp<NNEIB;neib_tmp++)
      			{

      			  //-find 3D indices for 2nd variable for neighbor grid points
              neib = (neib_type)neib_tmp;
      			  Neighbor(neib, ix1, iy1, iz1, &ix2, &iy2, &iz2);

			  
      			  //if (ij2>=0){ //-neighbor exists... 
	     
              //-find 1D index for 2nd variable at this neighbor point
              ij2=Map2sv(ivar2,ix2,iy2,iz2, NVARS,NSMX,NSMXY);
              if (ij2>=0) //-neighbor exists..
              {
    			      //-record an entry in the Jacobian
    			      output << "Coupling of " << ij1 << " to " << ij2 << endl;
                if (ij1<0 || ij2 <0) SETERRQ2(1,"1st: ij1 %d, ij2 %d",ij1, ij2);
    			      jmatr[ij1][ij2]=1; 
    			    }
                        
	   			  
            }



  		    } //-iz1
    		} //-iy1
      } //-ix1
  	} //-ivar2
  } //-ivar1
 



  /*===Interaction between (ix1,iy1,iz1) and (ix2=ix1+-1,2;iy2=iy1;all iz)===*/
  /*due to the Fourier interpolation*/

  for (ivar1=1;ivar1<=NVARS;ivar1++)
  {
    for (ivar2=1;ivar2<=NVARS;ivar2++)
    {
      for (ix1=0;ix1<mesh->ngx;ix1++)
      {
        for (iy1=0;iy1<MYSUB;iy1++)
        {
          for (iz1=0;iz1<MZ;iz1++) 
          {
	     	      
  		      //-find 1D index for 1st variable at this grid point
  		      ij1=Map2sv(ivar1,ix1,iy1,iz1, NVARS,NSMX,NSMXY);


		      
  		      //-loop over radial neigbors on the stencil (2 or 4 first items in the list)
  		      for (neib_tmp=0;neib_tmp<2;neib_tmp++)
      			{
      			  //-find 3D indices for 2nd variable for neighbor grid points
              neib = (neib_type)neib_tmp;
      			  Neighbor(neib, ix1, iy1, iz1, &ix2, &iy2, &iz2);

      			  //-instead of this iz2 use whole toroidal line for this (ix2,iy2)
      			  for (iz2=0;iz2<MZ;iz2++)
    			    {

    			      //if (ij2>=0) //-neighbor exists... 
        				//{
        				  //-find 1D index for 2nd variable at this neighbor point
                  ij2=Map2sv(ivar2,ix2,iy2,iz2, NVARS,NSMX,NSMXY);
                  if (ij2>=0) { //-neighbor exists...
                    //-record an entry in the Jacobian
                    output << "Coupling of " << ij1 << " to " << ij2 << endl;
                    if (ij1<0 || ij2 <0) SETERRQ2(1,"2nd: ij1 %d, ij2 %d",ij1, ij2);
                    jmatr[ij1][ij2]=1; 
                  }
			      
    			    }
	   			  
      			}



  		    } //-iz1
    		} //-iy1
	    } //-ix1
  	} //-ivar2
  } //-ivar1




  /*====Interaction due to the vorticity inversion====================*/
  /*between each variable at (ix1,iy1,iz1) and vorticity for (all ix,iy2=iy1,all iz)*/

  for (ivar1=1;ivar1<=NVARS;ivar1++)
  {
    for (ix1=0;ix1<mesh->ngx;ix1++)
    {
      for (iy1=0;iy1<MYSUB;iy1++)
      {
        for (iz1=0;iz1<MZ;iz1++) 
        {
		  
    		  //-find 1D index for 1st variable at this grid point
    		  ij1=Map2sv(ivar1,ix1,iy1,iz1, NVARS,NSMX,NSMXY);
		  
		    
    		  for (ix2=0;ix2<mesh->ngx;ix2++)
  		    {
  		      for (iz2=0;iz2<MZ;iz2++) 
      			{
      			  //-find 1D index for 2nd variable at this neighbor point
      			  ivar2=1; //-vorticity
      			  ij2=Map2sv(ivar2,ix2,iy2,iz2, NVARS,NSMX,NSMXY);
			  
      			  //-record an entry in the Jacobian
              if (ij2>=0){
                output << "Coupling of " << ij1 << " to " << ij2 << endl;
                if (ij1<0 || ij2 <0) SETERRQ2(1,"3rd: ij1 %d, ij2 %d",ij1, ij2);
                jmatr[ij1][ij2]=1; 
              }
      			}
  		    }			
    		} //-iz1
	    } //-iy1
  	} //-ix1
  } //-ivar1

  return 0;
}

#endif
