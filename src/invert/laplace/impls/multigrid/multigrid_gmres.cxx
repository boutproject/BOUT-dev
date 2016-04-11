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

void LaplaceMultigrid::pGMRES(BoutReal *sol,BoutReal *rhs,int level,int iplag) {
  int i,k,ii,it,etest = 1,MAXIT;
  BoutReal ini_e,error,a0,a1,rederr,perror,EPS,AbsE, RHO;
  BoutReal **v,*p,*q,*r;
  BoutReal c[MAXGM+1],s[MAXGM+1],y[MAXGM+1],g[MAXGM+1],h[MAXGM+1][MAXGM+1];

  if((level == 0) || (iplag == 0)) MAXIT = 10000;
  else MAXIT = 200;

  int size = zNP*xNP;
  int ldim = (lnx[level]+2)*(lnz[level]+2);
  v = new BoutReal *[MAXGM+1];
  for(i=0;i<MAXGM+1;i++) v[i] = new BoutReal[ldim];
  p = new BoutReal[ldim];
  q = new BoutReal[ldim];
  r = new BoutReal[ldim];
  for(i = 0;i<ldim;i++) sol[i] = 0.0;
  int num = 0;

  communications(rhs,level);
  ini_e = sqrt(vectorProd(level,rhs,rhs));
  //  if(xProcI == 0) 
  //  printf("%d --  In gmres ini = %24.18f\n",size,ini_e);
  for(i = 0;i<ldim;i++) r[i] = 0.0;
  if(level > 0)  cycleMG(level,r,rhs); 
  else smoothings(level,r,rhs);
  for(i = 0;i < ldim;i++) v[0][i] = r[i];
  perror = ini_e;
  do{
    a1 = vectorProd(level,v[0],v[0]);
    a1 = sqrt(a1);
    if(fabs(a1) < atol)
      throw BoutException("a1 in GMRES is wrong \n");
    a0 = 1.0/a1;
    for(i=0;i<ldim;i++) v[0][i] *= a0;
    g[0] = a1;
    for(i=1;i<MAXGM+1;i++) g[i] = 0.0;
    for(it = 0;it<MAXGM;it++) {
      multiAVec(level,v[it],q);
      for(i=0;i<ldim;i++) v[it+1][i] = 0.0;
      if(level >0) cycleMG(level,v[it+1],q); 
      else smoothings(level,v[it+1],q);
      for(i=0;i<it+1;i++) h[i][it] = vectorProd(level,v[it+1],v[i]);
      for(i=0;i<it+1;i++) {
        a0 = -h[i][it];
        for(k=0;k<ldim;k++) v[it+1][k] += a0*v[i][k]; 
      }
      a1 = vectorProd(level,v[it+1],v[it+1]);
      a1 = sqrt(a1);
      if(fabs(a1) < atol)
        throw BoutException("a1 in GMRES is wrong \n");
      a0 = 1.0/a1;
      for(i=0;i<ldim;i++) v[it+1][i] *= a0;
      h[it+1][it] = a1;

      for(i=0;i<it;i++) {
        a0 = c[i]*h[i][it] -s[i]*h[i+1][it];
        a1 = s[i]*h[i][it] +c[i]*h[i+1][it];
        h[i][it] = a0;
        h[i+1][it] = a1;
      }
      a0 = h[it][it]*h[it][it] + h[it+1][it]*h[it+1][it];
      if(a0 < 0.0)
        throw BoutException("a0 in GMRES is negative \n");
      a0 = sqrt(a0);   
      if(fabs(a0) < atol)
        throw BoutException("a0 in GMRES is wrong \n");
      c[it] = h[it][it]/a0;
      s[it] = -h[it+1][it]/a0;
      h[it][it] = a0;
      h[it+1][it] = 0.0;
      a0 = c[it]*g[it]-s[it]*g[it+1];
      a1 = s[it]*g[it]+c[it]*g[it+1];
      g[it] = a0;
      g[it+1] = a1;
    
    /* Get solution y and x_m*/
      for(i=it;i>=0;i--) {
        y[i] = g[i];
        for(int j=i+1;j<=it;j++) y[i] -= h[i][j]*y[j];
        y[i] = y[i]/h[i][i];
      }
      for(i=0;i<ldim;i++) p[i] = sol[i];
      for(i=0;i<=it;i++) { 
        for(k=0;k<ldim;k++) p[k] += y[i]*v[i][k]; 
      }
    
      /* Get r_m and test convergence.*/
      residualVec(level,p,rhs,r);
      error = sqrt(vectorProd(level,r,r));
      num += 1;
      if((error >dtol) || num > MAXIT) {
        throw BoutException("Error in GMRES %16.10f (%d)\n",error,num);
        etest = 0; 
        break;
      }
      if(error <= rtol*ini_e+atol) {
        etest = 0;
        break;
      }
      if(fabs(perror-error)/error <rtol) {
        if(it == 0) etest = 0;
        num -= 1;
        break;
      }
      perror = error;
    }
  /* Restart with new initial */
    for(i = 0;i<ldim;i++) v[0][i] = 0.0;
    if(level > 0)  cycleMG(level,v[0],r); 
    else smoothings(level,v[0],r);
    for(i = 0;i<ldim;i++) sol[i] = p[i];
    if(num>MAXIT) {
      throw BoutException(" GMRES Iteration limit.\n");
      etest = 0;
    }
    //    if((etest == 1) & (xProcI == 0)) 
    //  printf("Restart GMRES  %d | %20.14f\n",num,error/ini_e);    
  } while(etest == 1); 

  if((xProcI == 0) & (pcheck == 1)) {
    rederr = log(error/ini_e)/((double)num);
    rederr = exp(rederr);      
    printf("The average error reduction of GMRES %d: %14.8f\n",num,rederr);
    fflush(stdout);
  }

  delete [] p;
  delete [] q;
  delete [] r;
  for(i=0;i<MAXGM+1;i++) {
    delete [] v[i];
  }
  delete [] v;

  return; 
}




BoutReal LaplaceMultigrid::vectorProd(int level,BoutReal* x,BoutReal* y) {
  // nx does not include guard cells
  
  BoutReal val;
  BoutReal ini_e = 0.0;
  for(int i= 1;i<lnx[level]+1;i++)
    for(int k=1;k<lnz[level]+1;k++) {
      int ii = i*(lnz[level]+2)+k;
      ini_e += x[ii]*y[ii];
    }
  if(xNP*zNP > 1) 
    MPI_Allreduce(&ini_e,&val,1,MPI_DOUBLE,MPI_SUM,commXZ);
  else val = ini_e;

  return(val);  
}

void LaplaceMultigrid::multiAVec(int level, BoutReal *x, BoutReal *b) {

  BoutReal val;
  int nn;
  communications(x,level);
  int mm = lnz[level]+2;
  for(int i = 1;i<lnx[level]+1;i++)
    for(int k=1;k<lnz[level]+1;k++) {
      nn = i*mm+k;
      b[nn] = matmg[level][nn*9+4]*x[nn] + matmg[level][nn*9+3]*x[nn-1]
	+matmg[level][nn*9+5]*x[nn+1] + matmg[level][nn*9+1]*x[nn-mm]
        +matmg[level][nn*9+7]*x[nn+mm] +matmg[level][nn*9]*x[nn-mm-1]
        +matmg[level][nn*9+2]*x[nn-mm+1] + matmg[level][nn*9+6]*x[nn+mm-1]
        +matmg[level][nn*9+8]*x[nn+mm+1];
    } 
  communications(b,level);
}

void LaplaceMultigrid::residualVec(int level, BoutReal *x, BoutReal *b,
BoutReal *r) {

  BoutReal val;
  int nn;
  communications(x,level);
  int mm = lnz[level]+2;
  for(int i = 1;i<lnx[level]+1;i++)
    for(int k=1;k<lnz[level]+1;k++) {
      nn = i*mm+k;
      val = matmg[level][nn*9+4]*x[nn] + matmg[level][nn*9+3]*x[nn-1]
	+matmg[level][nn*9+5]*x[nn+1] + matmg[level][nn*9+1]*x[nn-mm]
        +matmg[level][nn*9+7]*x[nn+mm] +matmg[level][nn*9]*x[nn-mm-1]
        +matmg[level][nn*9+2]*x[nn-mm+1] + matmg[level][nn*9+6]*x[nn+mm-1]
        +matmg[level][nn*9+8]*x[nn+mm+1];
      r[nn] = b[nn]-val;
    } 
  communications(r,level);

}
