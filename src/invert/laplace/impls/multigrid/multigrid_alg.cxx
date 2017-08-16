/**************************************************************************
 * Perpendicular Laplacian inversion. 
 *                           Using Geometrical Multigrid Solver
 *
 * Equation solved is:
 d*\nabla^2_\perp x + (1/c1)\nabla_perp c2\cdot\nabla_\perp x + a x = b
 *
 **************************************************************************
 * Copyright 2016 K.S. Kang
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


// Define basic multigrid algorithm

MultigridAlg::MultigridAlg(int level,int lx,int lz, int gx, int gz, 
			   MPI_Comm comm,int check) {

  mglevel = level;
  commMG = comm;
  pcheck = check;

  if(pcheck > 0) output<<"Construct MG "<<level<<endl; 

  /* Momory allocate for Multigrid */
  gnx = new int[mglevel];
  gnz = new int[mglevel];
  lnx = new int[mglevel];
  lnz = new int[mglevel];
  gnx[mglevel-1] = gx;
  gnz[mglevel-1] = gz;
  lnx[mglevel-1] = lx;
  lnz[mglevel-1] = lz;
  if(mglevel > 1) {
    for(int i=mglevel-1;i>0;i--) {
      gnx[i-1] = gnx[i]/2;
      gnz[i-1] = gnz[i]/2;
      lnx[i-1] = lnx[i]/2;
      lnz[i-1] = lnz[i]/2;
    }
  }

  matmg = new BoutReal *[mglevel];
  for(int i = 0;i<mglevel;i++) {
    matmg[i] = new BoutReal[(lnx[i]+2)*(lnz[i]+2)*9];
  }
}

MultigridAlg::~MultigridAlg() {
  output<<"End deconstruction Malg AAAA "<<numP<<endl;
}

void MultigridAlg::cleanMem() {
  // Finalize, deallocate memory, etc.
  for(int i = 0;i<mglevel;i++) delete [] matmg[i];
  delete [] matmg;
  delete [] lnz;
  delete [] lnx;
  delete [] gnz;
  delete [] gnx;
}

void MultigridAlg::getSolution(BoutReal *x,BoutReal *b,int flag) {

  if(flag == 0) {
    //Solve exaclty
    if(mglevel == 1) pGMRES(x,b,mglevel-1,1);
    else if(mgplag == 1) pGMRES(x,b,mglevel-1,1);
    else solveMG(x,b,mglevel-1);
  }
  else {
    cycleMG(mglevel-1,x,b);
    if(flag > 1) {
      BoutReal *y,*r;
      int level = mglevel-1;
      int ldim = (lnx[level]+2)*(lnz[level]+2);
      y = new BoutReal[ldim];
      r = new BoutReal[ldim];
      for(int n = 1;n<flag;n++) {
        residualVec(level,x,b,r);
#pragma omp parallel default(shared)
#pragma omp for
        for(int i = 0;i<ldim;i++) y[i] = 0.0;
        cycleMG(level,y,r);
#pragma omp parallel default(shared)
#pragma omp for
        for(int i = 0;i<ldim;i++) x[i] = x[i]+y[i];
      }
      delete [] r;
      delete [] y;
    }
  }
}


void MultigridAlg::cycleMG(int level,BoutReal *sol,BoutReal *rhs)
{
  int i;

  if(level == 0) {
    lowestSolver(sol,rhs,0);
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

#pragma omp parallel default(shared) private(i)
#pragma omp for
    for(i=0;i<(lnx[level-1]+2)*(lnz[level-1]+2);i++) y[i] = 0.0;
  
    cycleMG(level-1,y,pr);

    prolongation(level-1,y,iy);
#pragma omp parallel default(shared) private(i)
#pragma omp for
    for(i=0;i<(lnx[level]+2)*(lnz[level]+2);i++) 
       sol[i] += iy[i];

    smoothings(level,sol,rhs);

    delete [] iy;
    delete [] y;
    delete [] pr;
    delete [] r;

  }
  communications(sol,level);
  return;
}

void MultigridAlg::projection(int level,BoutReal *r,BoutReal *pr) 
{

  int nn,n0,n1,n2,n3;
  communications(r,level);
  for(int i=0;i<(lnx[level-1]+2)*(lnz[level-1]+2);i++) pr[i] = 0.;
  for (int i=1; i<lnx[level-1]+1; i++) {
    int i2 = 2*i-1;
#pragma omp parallel default(shared) private(nn,n0,n1,n2,n3)
#pragma omp for
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

void MultigridAlg::prolongation(int level,BoutReal *x,BoutReal *ix) {

  int nn,n0,n1,n2,n3;
  communications(x,level);
#pragma omp parallel default(shared)
#pragma omp for
  for(int i=0;i<(lnx[level+1]+2)*(lnz[level+1]+2);i++) ix[i] = 0.;
  for (int i=1; i<lnx[level]+1; i++) {
    int i2 = 2*i-1;
#pragma omp parallel default(shared) private(nn,n0,n1,n2,n3)
#pragma omp for
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

void MultigridAlg::smoothings(int level, BoutReal *x, BoutReal *b) {

  BoutReal val,*x0;
  int nn,dim;
  int mm = lnz[level]+2;
  dim = mm*(lnx[level]+2);
  if(mgsm == 0) {
    x0 = new BoutReal[dim];
    communications(x,level);
    for(int num =0;num < 2;num++) {
#pragma omp parallel default(shared)
#pragma omp for
      for(int i = 0;i<dim;i++) x0[i] = x[i];    
      for(int i = 1;i<lnx[level]+1;i++)
#pragma omp parallel default(shared) private(nn)
#pragma omp for
        for(int k=1;k<lnz[level]+1;k++) {
          nn = i*mm+k;
          val = b[nn] - matmg[level][nn*9+3]*x0[nn-1]
	   - matmg[level][nn*9+5]*x0[nn+1] - matmg[level][nn*9+1]*x0[nn-mm]
           - matmg[level][nn*9+7]*x0[nn+mm] - matmg[level][nn*9]*x0[nn-mm-1]
           - matmg[level][nn*9+2]*x0[nn-mm+1] - matmg[level][nn*9+6]*x0[nn+mm-1]
           - matmg[level][nn*9+8]*x0[nn+mm+1];
          if(fabs(matmg[level][nn*9+4]) <atol)
            throw BoutException("Error at matmg(%d-%d)",level,nn);

          x[nn] = (1.0-omega)*x[nn] + omega*val/matmg[level][nn*9+4];
        } 
      communications(x,level);
    }
    delete [] x0;
  }
  else {
    communications(x,level);    
    for(int i = 1;i<lnx[level]+1;i++)
#pragma omp parallel default(shared) private(nn)
#pragma omp for
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
#pragma omp parallel default(shared) private(nn)
#pragma omp for
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

void MultigridAlg::pGMRES(BoutReal *sol,BoutReal *rhs,int level,int iplag) {
  int k,it,etest = 1,MAXIT;
  BoutReal ini_e,error,a0,a1,rederr,perror;
  BoutReal **v,*p,*q,*r;
  BoutReal c[MAXGM+1],s[MAXGM+1],y[MAXGM+1],g[MAXGM+1],h[MAXGM+1][MAXGM+1];

  if((level == 0) || (iplag == 0)) MAXIT = 40000;
  else MAXIT = 500;

  int ldim = (lnx[level]+2)*(lnz[level]+2);
  v = new BoutReal *[MAXGM+1];
  for(int i=0;i<MAXGM+1;i++) v[i] = new BoutReal[ldim];
  p = new BoutReal[ldim];
  q = new BoutReal[ldim];
  r = new BoutReal[ldim];
#pragma omp parallel default(shared)
#pragma omp for
  for(int i = 0;i<ldim;i++) sol[i] = 0.0;
  int num = 0;

  communications(rhs,level);
  ini_e = sqrt(vectorProd(level,rhs,rhs));
  if((pcheck == 1) && (rProcI == 0)) {
    output<<numP<<"--In GMRES ini "<<ini_e<<endl;
  }
  if(ini_e < atol*rtol) {
    if((pcheck == 1) && (rProcI == 0)) {
      output<<numP<<"Don't need to solve. E= "<<ini_e<<endl;
    }
    return;
  }
#pragma omp parallel default(shared)
#pragma omp for
  for(int i = 0;i<ldim;i++) r[i] = 0.0;
  if(iplag ==  0)  smoothings(level,r,rhs);
  else  cycleMG(level,r,rhs); 
#pragma omp parallel default(shared)
#pragma omp for
  for(int i = 0;i < ldim;i++) v[0][i] = r[i];
  perror = ini_e;
  do{
    a1 = vectorProd(level,v[0],v[0]);
    a1 = sqrt(a1);
    if(fabs(a1) < atol*rtol) {
      output<<num<<"First a1 in GMRES is wrong "<<a1<<":"<<level<<endl;
    }
    a0 = 1.0/a1;
#pragma omp parallel default(shared)
#pragma omp for
    for(int i=0;i<ldim;i++) v[0][i] *= a0;
    g[0] = a1;
    for(int i=1;i<MAXGM+1;i++) g[i] = 0.0;
    for(it = 0;it<MAXGM;it++) {
      multiAVec(level,v[it],q);
#pragma omp parallel default(shared)
#pragma omp for
      for(int i=0;i<ldim;i++) v[it+1][i] = 0.0;

      if(iplag == 0)  smoothings(level,v[it+1],q);
      else cycleMG(level,v[it+1],q); 

      for(int i=0;i<it+1;i++) h[i][it] = vectorProd(level,v[it+1],v[i]);
      for(int i=0;i<it+1;i++) {
        a0 = -h[i][it];
        for(k=0;k<ldim;k++) v[it+1][k] += a0*v[i][k]; 
      }
      a1 = vectorProd(level,v[it+1],v[it+1]);
      a1 = sqrt(a1);
      if(fabs(a1) < atol*rtol) {
        output<<num<<"In Second a1 in GMRES is wrong "<<a1<<endl;
      }
      a0 = 1.0/a1;
#pragma omp parallel default(shared)
#pragma omp for
      for(int i=0;i<ldim;i++) v[it+1][i] *= a0;
      h[it+1][it] = a1;

      for(int i=0;i<it;i++) {
        a0 = c[i]*h[i][it] -s[i]*h[i+1][it];
        a1 = s[i]*h[i][it] +c[i]*h[i+1][it];
        h[i][it] = a0;
        h[i+1][it] = a1;
      }
      a0 = h[it][it]*h[it][it] + h[it+1][it]*h[it+1][it];
      if(a0 < 0.0)
        throw BoutException("a0 in GMRES is negative \n");
      a0 = sqrt(a0);   
      if(fabs(a0) < atol*rtol)
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
      for(int i=it;i>=0;i--) {
        y[i] = g[i];
        for(int j=i+1;j<=it;j++) y[i] -= h[i][j]*y[j];
        y[i] = y[i]/h[i][i];
      }
#pragma omp parallel default(shared)
#pragma omp for
      for(int i=0;i<ldim;i++) p[i] = sol[i];
      for(int i=0;i<=it;i++) { 
#pragma omp parallel default(shared)
#pragma omp for
        for(int k=0;k<ldim;k++) p[k] += y[i]*v[i][k]; 
      }
    
      /* Get r_m and test convergence.*/
      residualVec(level,p,rhs,r);
      error = sqrt(vectorProd(level,r,r));
      num += 1;
      if((error > dtol) || num > MAXIT) {
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
#pragma omp parallel default(shared)
#pragma omp for
    for(int i = 0;i<ldim;i++) v[0][i] = 0.0;
    if(iplag ==  0)  smoothings(level,v[0],r);
    else cycleMG(level,v[0],r); 
    
#pragma omp parallel default(shared)
#pragma omp for
    for(int i = 0;i<ldim;i++) sol[i] = p[i];
    if(num>MAXIT) {
      throw BoutException(" GMRES Iteration limit.\n");
      etest = 0;
    }
    //    if((etest == 1) & (xProcI == 0)) 
    //  printf("Restart GMRES  %d | %20.14f\n",num,error/ini_e);    
  } while(etest == 1); 

  if(level == 0) {
    if((rProcI == 0) && (pcheck == 1)) {
      rederr = log(error/ini_e)/((double)num);
      rederr = exp(rederr);      
      printf("The average error reduction of GMRES %d: %14.8f(%18.10f)\n",num,rederr,error);
      fflush(stdout);
    }
  }
  else {
    if((rProcI == 0) && (pcheck == 1)) {
      rederr = log(error/ini_e)/((double)num);
      rederr = exp(rederr);      
      printf("The average error reduction of PGMRES %d: %14.8f(%18.10f)\n",num,rederr,error);
      fflush(stdout);
    }
  }
  delete [] p;
  delete [] q;
  delete [] r;
  for(int i=0;i<MAXGM+1;i++) {
    delete [] v[i];
  }
  delete [] v;

  return; 
}

void MultigridAlg::setMultigridC(int UNUSED(plag)) {

  int level = mglevel - 1;
  for(int n = level;n>0;n--) {
    setMatrixC(n);
    if(pcheck == 2) {
      output<<n<<"matrix in Alg = "<<lnx[n]<<","<<lnz[n]<<endl;
      output<<gnx[n]<<","<<gnz[n]<<endl;
    }
  }
}

void MultigridAlg::lowestSolver(BoutReal *x, BoutReal *b, int UNUSED(plag)) {
  pGMRES(x, b, 0, 0);
}

BoutReal MultigridAlg::vectorProd(int level,BoutReal* x,BoutReal* y) {
  // nx does not include guard cells
  
  BoutReal val;
  BoutReal ini_e = 0.0;
  for(int i= 1;i<lnx[level]+1;i++)
#pragma omp parallel default(shared) 
#pragma omp for reduction(+:ini_e)
    for(int k=1;k<lnz[level]+1;k++) {
      int ii = i*(lnz[level]+2)+k;
      ini_e += x[ii]*y[ii];
    }
  if(numP > 1) 
    MPI_Allreduce(&ini_e,&val,1,MPI_DOUBLE,MPI_SUM,commMG);
  else val = ini_e;

  return(val);  
}

void MultigridAlg::multiAVec(int level, BoutReal *x, BoutReal *b) {

  communications(x,level);
  int mm = lnz[level]+2;
#pragma omp parallel default(shared)
#pragma omp for
  for(int i = 0;i<mm*(lnx[level]+2);i++) b[i] = 0.0;
  for(int i = 1;i<lnx[level]+1;i++)
#pragma omp parallel default(shared)
#pragma omp for
    for(int k=1;k<lnz[level]+1;k++) {
      int nn = i*mm+k;
      b[nn] = matmg[level][nn*9+4]*x[nn] + matmg[level][nn*9+3]*x[nn-1]
	+matmg[level][nn*9+5]*x[nn+1] + matmg[level][nn*9+1]*x[nn-mm]
        +matmg[level][nn*9+7]*x[nn+mm] +matmg[level][nn*9]*x[nn-mm-1]
        +matmg[level][nn*9+2]*x[nn-mm+1] + matmg[level][nn*9+6]*x[nn+mm-1]
        +matmg[level][nn*9+8]*x[nn+mm+1];
    } 
  communications(b,level);
}

void MultigridAlg::residualVec(int level, BoutReal *x, BoutReal *b,
BoutReal *r) {

  BoutReal val;
  int mm;
  communications(x,level);
  mm = lnz[level]+2;
#pragma omp parallel default(shared)
#pragma omp for
  for(int i = 0;i<mm*(lnx[level]+2);i++) r[i] = 0.0;
  for(int i = 1;i<lnx[level]+1;i++)
#pragma omp parallel default(shared)
#pragma omp for
    for(int k=1;k<lnz[level]+1;k++) {
      int nn = i*mm+k;
      val = matmg[level][nn*9+4]*x[nn] + matmg[level][nn*9+3]*x[nn-1]
	+matmg[level][nn*9+5]*x[nn+1] + matmg[level][nn*9+1]*x[nn-mm]
        +matmg[level][nn*9+7]*x[nn+mm] +matmg[level][nn*9]*x[nn-mm-1]
        +matmg[level][nn*9+2]*x[nn-mm+1] + matmg[level][nn*9+6]*x[nn+mm-1]
        +matmg[level][nn*9+8]*x[nn+mm+1];
      r[nn] = b[nn]-val;
    } 
  communications(r,level);

}

void MultigridAlg::setMatrixC(int level) {

  BoutReal ratio = 8.0; 

#pragma omp parallel default(shared)
#pragma omp for
  for(int i=0;i<(lnx[level-1]+2)*(lnz[level-1]+2)*9;i++) { 
    matmg[level-1][i] = 0.0;
  }
  for(int i = 1;i<lnx[level-1]+1;i++) {
    int i2 = 2*i-1;
#pragma omp parallel default(shared)
#pragma omp for
    for(int k = 1;k<lnz[level-1]+1;k++) {
      int k2 = 2*k-1;
      int mm = i*(lnz[level-1]+2)+k;
      int m0 = i2*(lnz[level]+2)+k2;
      int m1 = i2*(lnz[level]+2)+k2+1;
      int m2 = (i2+1)*(lnz[level]+2)+k2;
      int m3 = (i2+1)*(lnz[level]+2)+k2+1;
      BoutReal val = matmg[level][m0*9+4]+matmg[level][m1*9+4];
      val += matmg[level][m2*9+4] + matmg[level][m3*9+4];
      val += matmg[level][m0*9+5] + matmg[level][m1*9+3];
      val += matmg[level][m2*9+5] + matmg[level][m3*9+3];
      val += matmg[level][m0*9+7] + matmg[level][m2*9+1];
      val += matmg[level][m1*9+7] + matmg[level][m3*9+1];
      val += matmg[level][m0*9+8] + matmg[level][m3*9];
      val += matmg[level][m1*9+6] + matmg[level][m2*9+2];
      matmg[level-1][mm*9+4] = val/ratio;
      val = matmg[level][m0*9+1]+matmg[level][m1*9+1];
      val += matmg[level][m0*9+2]+matmg[level][m1*9];
      matmg[level-1][mm*9+1] = val/ratio;
      val = matmg[level][m0*9+3]+matmg[level][m2*9+3];
      val += matmg[level][m0*9+6]+matmg[level][m2*9];
      matmg[level-1][mm*9+3] = val/ratio;
      val = matmg[level][m1*9+5]+matmg[level][m3*9+5];
      val += matmg[level][m1*9+8]+matmg[level][m3*9+2];
      matmg[level-1][mm*9+5] = val/ratio;
      val = matmg[level][m2*9+7]+matmg[level][m3*9+7];
      val += matmg[level][m2*9+8]+matmg[level][m3*9+6];
      matmg[level-1][mm*9+7] = val/ratio;
      matmg[level-1][mm*9] = matmg[level][m0*9]/ratio;
      matmg[level-1][mm*9+2] = matmg[level][m1*9+2]/ratio;
      matmg[level-1][mm*9+6] = matmg[level][m2*9+6]/ratio;
      matmg[level-1][mm*9+8] = matmg[level][m3*9+8]/ratio;      
    }
  }

}

void MultigridAlg::communications(BoutReal* x, int level) {
 
  MPI_Status  status[4];
  int stag,rtag,ierr;

  if(zNP > 1) {
    MPI_Datatype xvector;
    //    output<<"Start Z-comm"<<endl;
    ierr = MPI_Type_vector(lnx[level], 1, lnz[level]+2, MPI_DOUBLE, &xvector);
    ASSERT1(ierr == MPI_SUCCESS);
    
    ierr = MPI_Type_commit(&xvector);
    ASSERT1(ierr == MPI_SUCCESS);
    
    // Send to z+ and recieve from z-
    stag = rProcI;
    rtag = zProcM;
    // output<<"before to z+:"<<stag<<":"<<rtag<<endl;
    ierr = MPI_Sendrecv(&x[2*(lnz[level]+2)-2],1,xvector,zProcP,stag,
                        &x[lnz[level]+2],1,xvector,zProcM,rtag,commMG,status);
    ASSERT1(ierr == MPI_SUCCESS);
    // Send to z- and recieve from z+
    stag = rProcI+numP;
    rtag = zProcP+numP;
    // output<<"before to z-:"<<stag<<":"<<rtag<<endl;
    ierr = MPI_Sendrecv(&x[lnz[level]+3],1,xvector,zProcM,stag,
                        &x[2*(lnz[level]+2)-1],1,xvector,zProcP,rtag,
                        commMG,status);
    ASSERT1(ierr == MPI_SUCCESS);
    
    ierr = MPI_Type_free(&xvector);
    ASSERT1(ierr == MPI_SUCCESS);
  } else {
    for (int i=1;i<lnx[level]+1;i++) {
      x[i*(lnz[level]+2)] = x[(i+1)*(lnz[level]+2)-2];
      x[(i+1)*(lnz[level]+2)-1] = x[i*(lnz[level]+2)+1];
    }
  }
  if (xNP > 1) {
    // Send to x+ and recieve from x-
    stag = rProcI; 
    rtag = xProcM;
    ierr = MPI_Sendrecv(&x[lnx[level]*(lnz[level]+2)],lnz[level]+2,
                MPI_DOUBLE,xProcP,stag,&x[0],lnz[level]+2,MPI_DOUBLE,
                xProcM,rtag,commMG,status);
    ASSERT1(ierr == MPI_SUCCESS);
    
    // Send to x- and recieve from x+
    stag = rProcI+xNP;
    rtag = xProcP+xNP;;
    ierr = MPI_Sendrecv(&x[lnz[level]+2],lnz[level]+2,MPI_DOUBLE,xProcM,stag,
                        &x[(lnx[level]+1)*(lnz[level]+2)],lnz[level]+2,
                        MPI_DOUBLE,xProcP,rtag,commMG,status);
    ASSERT1(ierr == MPI_SUCCESS);
  }  else {
    for (int i=0;i<lnz[level]+2;i++) {
      x[i] = x[lnx[level]*(lnz[level]+2)+i];
      x[(lnx[level]+1)*(lnz[level]+2)+i] = x[(lnz[level]+2)+i];
    }
  }
}


void MultigridAlg::solveMG(BoutReal *sol,BoutReal *rhs,int level) {
  int m,MAXIT = 150;
  BoutReal ini_e,perror,error,rederr;
  BoutReal *y,*r;
  
  int ldim = (lnx[level]+2)*(lnz[level]+2);

#pragma omp parallel default(shared)
#pragma omp for
  for(int i = 0;i<ldim;i++) sol[i] = 0.0;

  communications(rhs,level);
  ini_e = vectorProd(level,rhs,rhs);
  if(ini_e < 0.0) {
    throw BoutException("In MG Initial Error %10.4e \n",ini_e);
  }
  ini_e = sqrt(ini_e);
  if((pcheck == 1) && (rProcI == 0)) 
    printf("%d \n  In MGsolve ini = %24.18f\n",numP,ini_e);
  y = new BoutReal[ldim];
  r = new BoutReal[ldim];
#pragma omp parallel default(shared)
#pragma omp for
  for(int i = 0;i<ldim;i++) r[i] = rhs[i];

  perror = ini_e;
  for(m=0;m<MAXIT;m++) {
#pragma omp parallel default(shared)
#pragma omp for
    for(int i = 0;i<ldim;i++) y[i] = 0.0;
    cycleMG(level,y,r);
#pragma omp parallel default(shared)
#pragma omp for
    for(int i = 0;i<ldim;i++) sol[i] = sol[i]+y[i];
    residualVec(level,sol,rhs,r);
    error = sqrt(vectorProd(level,r,r));
    if((pcheck == 1) && (rProcI == 0)) 
      printf("%d \n  In MGsolve error = %24.18f\n",m,error);
    if(error < rtol*ini_e+atol) break;
    if((fabs(perror-error)/error <rtol) || (error > dtol)) {
      throw BoutException("In MG Limited Error %10.4e \n",error);
      m-= 1;
      break;
    }
    perror = error;
  }

  if((rProcI == 0) && (pcheck == 1)) {
    rederr = log(error/ini_e)/((double)m+1.0);
    rederr = exp(rederr); 
    if(m == MAXIT) 
      printf("Reached maximum iteration: %14.8f\n",error);
    printf("The average error reduction of MG %d: %14.8f(%18.12f)\n",m+1,rederr,error);
    fflush(stdout);
  }

  delete [] r;
  delete [] y;
  return;
}
