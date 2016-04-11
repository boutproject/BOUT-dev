/**************************************************************************
 * Perpendicular Laplacian inversion. 
 *                           Using Geometrical Multigrid Solver
 *
 * Equation solved is:
 d*\nabla^2_\perp x + (1/c1)\nabla_perp c2\cdot\nabla_\perp x + a x = b
 *
 **************************************************************************
 * Copyright 2015 K.S. Kang
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 * 
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#include "multigrid_laplace.hxx"

void LaplaceMultigrid::solveMG(BoutReal *sol,BoutReal *rhs,int level) {
  int i,m,MAXIT = 50;
  BoutReal ini_e,perror,error,rederr;
  BoutReal *y,*r;
  
  int size = zNP*xNP;
  int ldim = (lnx[level]+2)*(lnz[level]+2);

  for(i = 0;i<ldim;i++) sol[i] = 0.0;
  int num = 0;

  communications(rhs,level);
  ini_e = sqrt(vectorProd(level,rhs,rhs));
  //  if(xProcI == 0) 
  //  printf("%d \n  In MGsolve ini = %24.18f\n",size,ini_e);
  y = new BoutReal[ldim];
  r = new BoutReal[ldim];
  for(i = 0;i<ldim;i++) r[i] = rhs[i];

  perror = ini_e;
  for(m=0;m<MAXIT;m++) {
    for(i = 0;i<ldim;i++) y[i] = 0.0;
    cycleMG(level,y,r);
    for(i = 0;i<ldim;i++) sol[i] = sol[i]+y[i];
    residualVec(level,sol,rhs,r);
    error = sqrt(vectorProd(level,r,r));
    if(error < rtol*ini_e+atol) break;
    if((fabs(perror-error)/error <rtol) || (error > dtol)) {
      throw BoutException("In MG Limited Error %10.4e \n",error);
      m-= 1;
      break;
    }
    perror = error;
  }

  if((xProcI == 0) & (pcheck == 1)) {
    rederr = log(error/ini_e)/((double)m+1.0);
    rederr = exp(rederr);  
    printf("The average error reduction of %d: %14.8f\n",m+1,rederr);
    fflush(stdout);
  }

  delete [] r;
  delete [] y;
  return;
}

void LaplaceMultigrid::cycleMG(int level,BoutReal *sol,BoutReal *rhs)
{
  int i,sdof,ldim;
  MPI_Comm mycomm;

  if(level == 0) {
    pGMRES(sol,rhs,level,0);
  }
  else {
    BoutReal *r,*pr,*y,*iy;
    r = new BoutReal[(lnx[level]+2)*(lnz[level]+2)];
    pr = new BoutReal[(lnx[level-1]+2)*(lnz[level-1]+2)];
    y = new BoutReal[(lnx[level-1]+2)*(lnz[level-1]+2)];
    iy = new BoutReal[(lnx[level]+2)*(lnz[level]+2)];
    smoothings(level,sol,rhs);
    residualVec(level,sol,rhs,r);
    projection(level,r,pr);

    for(i=0;i<(lnx[level-1]+2)*(lnz[level-1]+2);i++) y[i] = 0.0;
    cycleMG(level-1,y,pr);

    prolongation(level-1,y,iy);
    for(i=0;i<(lnx[level]+2)*(lnz[level]+2);i++) 
       sol[i] += iy[i];

    smoothings(level,sol,rhs);

    delete [] iy;
    delete [] y;
    delete [] pr;
    delete [] r;

  }
  return;
}

void LaplaceMultigrid::projection(int level,BoutReal *r,BoutReal *pr) 
{

  int nn,n0,n1,n2,n3;
  communications(r,level);
  for(int i=0;i<(lnx[level-1]+2)*(lnz[level-1]+2);i++) pr[i] = 0.;
  for (int i=1; i<lnx[level-1]+1; i++) {
    int i2 = 2*i-1;
    for (int k=1; k<lnz[level-1]+1; k++) {
      int k2 = 2*k-1;
      nn = i*(lnz[level-1]+2)+k;
      n0 = i2*(lnz[level]+2)+k2;
      n1 = n0 + 1;
      n2 = n0 + lnz[level]+2;
      n3 = n2 + 1;
      pr[nn] = (r[n0]+r[n1]+r[n2]+r[n3])/4.0;
    }
  }
  communications(pr,level-1);
  return;
}

void LaplaceMultigrid::prolongation(int level,BoutReal *x,BoutReal *ix) {

  int nn,n0,n1,n2,n3;
  communications(x,level);
  for(int i=0;i<(lnx[level+1]+2)*(lnz[level+1]+2);i++) ix[i] = 0.;
  for (int i=1; i<lnx[level]+1; i++) {
    int i2 = 2*i-1;
    for (int k=1; k<lnz[level]+1; k++) {
      int k2 = 2*k-1;
      nn = i*(lnz[level]+2)+k;
      n0 = i2*(lnz[level+1]+2)+k2;
      n1 = n0 + 1;
      n2 = n0 + lnz[level+1]+2;
      n3 = n2 +1;
      ix[n0] = x[nn];
      ix[n1] = x[nn];
      ix[n2] = x[nn];
      ix[n3] = x[nn];
    }
  }
  communications(ix,level+1);
  return;
}

void LaplaceMultigrid::smoothings(int level, BoutReal *x, BoutReal *b) {

  BoutReal val,*x0;
  int nn,dim;
  int mm = lnz[level]+2;
  dim = mm*lnx[level]+2;
  if(mgsm == 0) {
    x0 = new BoutReal[dim];
    communications(x,level);
    for(int num =0;num < 2;num++) {
      for(int i = 0;i<dim;i++) x0[i] = x[i];    
      for(int i = 1;i<lnx[level]+1;i++)
        for(int k=1;k<lnz[level]+1;k++) {
          nn = i*mm+k;
          val = b[nn] - matmg[level][nn*9+3]*x0[nn-1]
	    - matmg[level][nn*9+5]*x0[nn+1] - matmg[level][nn*9+1]*x0[nn-mm]
            - matmg[level][nn*9+7]*x0[nn+mm] - matmg[level][nn*9]*x0[nn-mm-1]
            - matmg[level][nn*9+2]*x0[nn-mm+1] - matmg[level][nn*9+6]*x0[nn+mm-1]
            - matmg[level][nn*9+8]*x0[nn+mm+1];
          if(fabs(matmg[level][nn*9+4]) <atol)
            throw BoutException("Error at matmg(%d-%d)",level,nn);

          x[nn] = 0.3*x[nn]+0.7*val/matmg[level][nn*9+4];
        } 
      communications(x,level);
    }
    delete [] x0;
  }
  else {
    communications(x,level);    
    for(int i = 1;i<lnx[level]+1;i++)
      for(int k=1;k<lnz[level]+1;k++) {
        nn = i*mm+k;
        val = b[nn] - matmg[level][nn*9+3]*x[nn-1]
	    - matmg[level][nn*9+5]*x[nn+1] - matmg[level][nn*9+1]*x[nn-mm]
            - matmg[level][nn*9+7]*x[nn+mm] - matmg[level][nn*9]*x[nn-mm-1]
            - matmg[level][nn*9+2]*x[nn-mm+1] - matmg[level][nn*9+6]*x[nn+mm-1]
            - matmg[level][nn*9+8]*x[nn+mm+1];
          if(fabs(matmg[level][nn*9+4]) <atol)
            throw BoutException("Error at matmg(%d-%d)",level,nn);
        x[nn] = val/matmg[level][nn*9+4];
      } 
    communications(x,level);
    for(int i = lnx[level];i>0;i--)
      for(int k= lnz[level];k>0;k--) {
        nn = i*mm+k;
        val = b[nn] - matmg[level][nn*9+3]*x[nn-1]
	    - matmg[level][nn*9+5]*x[nn+1] - matmg[level][nn*9+1]*x[nn-mm]
            - matmg[level][nn*9+7]*x[nn+mm] - matmg[level][nn*9]*x[nn-mm-1]
            - matmg[level][nn*9+2]*x[nn-mm+1] - matmg[level][nn*9+6]*x[nn+mm-1]
            - matmg[level][nn*9+8]*x[nn+mm+1];
          if(fabs(matmg[level][nn*9+4]) <atol)
            throw BoutException("Error at matmg(%d-%d)",level,nn);
        x[nn] = val/matmg[level][nn*9+4];
      } 
    communications(x,level);
  }
  return;
}
