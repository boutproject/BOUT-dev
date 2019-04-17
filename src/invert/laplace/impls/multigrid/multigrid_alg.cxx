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
#include <bout/openmpwrap.hxx>
#include "unused.hxx"

// Define basic multigrid algorithm

MultigridAlg::MultigridAlg(int level, int lx, int lz, int gx, int gz, MPI_Comm comm,
                           int check)
    : mglevel(level), pcheck(check), commMG(comm) {

  if(pcheck > 0) output<<"Construct MG "<<level<<endl; 

  // Memory allocate for Multigrid
  gnx.reallocate(mglevel);
  gnz.reallocate(mglevel);
  lnx.reallocate(mglevel);
  lnz.reallocate(mglevel);

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

  // Could be replaced with a Matrix
  matmg = new BoutReal *[mglevel];
  for(int i = 0;i<mglevel;i++) {
    matmg[i] = new BoutReal[(lnx[i]+2)*(lnz[i]+2)*9];
  }
}

MultigridAlg::~MultigridAlg() {
  output<<"End deconstruction Malg AAAA "<<numP<<endl;
  for(int i = 0;i<mglevel;i++) delete [] matmg[i];
  delete [] matmg;
}

void MultigridAlg::getSolution(BoutReal *x,BoutReal *b,int flag) {

  // swap ghost cells of initial guess
  communications(x, mglevel-1);
  // don't think the ghost cells of the rhs are used, so don't need to be communicated
  // /JTO 17-4-2019

  if(flag == 0) {
    //Solve exaclty
    if(mglevel == 1) pGMRES(x,b,mglevel-1,1);
    else if(mgplag == 1) pGMRES(x,b,mglevel-1,1);
    else solveMG(x,b,mglevel-1);
  }
  else {
    cycleMG(mglevel-1,x,b);
    if(flag > 1) {
      int level = mglevel-1;
      int ldim = (lnx[level]+2)*(lnz[level]+2);
      Array<BoutReal> y(ldim);
      Array<BoutReal> r(ldim);
      for(int n = 1;n<flag;n++) {
        residualVec(level, x, b, std::begin(r));
        BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
        for(int i = 0;i<ldim;i++) y[i] = 0.0;
        cycleMG(level, std::begin(y), std::begin(r));
        BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
        for(int i = 0;i<ldim;i++) x[i] = x[i]+y[i];
      }
    }
  }
}


void MultigridAlg::cycleMG(int level,BoutReal *sol,BoutReal *rhs)
{
  if(level == 0) {
    lowestSolver(sol,rhs,0);
  }
  else {
    Array<BoutReal> r((lnx[level] + 2) * (lnz[level] + 2));
    Array<BoutReal> pr((lnx[level - 1] + 2) * (lnz[level - 1] + 2));
    Array<BoutReal> y((lnx[level - 1] + 2) * (lnz[level - 1] + 2));
    Array<BoutReal> iy((lnx[level] + 2) * (lnz[level] + 2));

    smoothings(level,sol,rhs);

    residualVec(level, sol, rhs, std::begin(r));

    projection(level, std::begin(r), std::begin(pr));

    BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
    for(int i=0;i<(lnx[level-1]+2)*(lnz[level-1]+2);i++) y[i] = 0.0;

    cycleMG(level - 1, std::begin(y), std::begin(pr));

    prolongation(level - 1, std::begin(y), std::begin(iy));
    BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
    for(int i=0;i<(lnx[level]+2)*(lnz[level]+2);i++) 
       sol[i] += iy[i];

    smoothings(level,sol,rhs);
  }
  return;
}

void MultigridAlg::projection(int level,BoutReal *r,BoutReal *pr) 
{

BOUT_OMP(parallel default(shared))
  {
BOUT_OMP(for)
    for(int i=0;i<(lnx[level-1]+2)*(lnz[level-1]+2);i++) pr[i] = 0.;
    int xend = lnx[level-1]+1;
    int zend = lnz[level-1]+1;
BOUT_OMP(for collapse(2))
    for (int i=1; i< xend; i++) {
      for (int k=1; k< zend; k++) {
        int i2 = 2*i-1;
        int k2 = 2*k-1;
        int nn = i*(lnz[level-1]+2)+k;
        int n0 = i2*(lnz[level]+2)+k2;
        int n1 = n0 + 1;
        int n2 = n0 + lnz[level]+2;
        int n3 = n2 + 1;
        pr[nn] = (r[n0]+r[n1]+r[n2]+r[n3])/4.0;
      }
    }
  }
  communications(pr,level-1);
  return;
}

void MultigridAlg::prolongation(int level,BoutReal *x,BoutReal *ix) {

BOUT_OMP(parallel default(shared))
  {
BOUT_OMP(for)
    for(int i=0;i<(lnx[level+1]+2)*(lnz[level+1]+2);i++) ix[i] = 0.;

    int xend = lnx[level]+1;
    int zend = lnz[level]+1;
BOUT_OMP(for collapse(2))
    for (int i=1; i< xend; i++) {
      for (int k=1; k< zend; k++) {
        int i2 = 2*i-1;
        int k2 = 2*k-1;
        int nn = i*(lnz[level]+2)+k;
        int n0 = i2*(lnz[level+1]+2)+k2;
        int n1 = n0 + 1;
        int n2 = n0 + lnz[level+1]+2;
        int n3 = n2 +1;
        ix[n0] = x[nn];
        ix[n1] = x[nn];
        ix[n2] = x[nn];
        ix[n3] = x[nn];
      }
    }
  }
  communications(ix,level+1);
  return;
}

void MultigridAlg::smoothings(int level, BoutReal *x, BoutReal *b) {

  int dim;
  int mm = lnz[level]+2;
  dim = mm*(lnx[level]+2);
  if(mgsm == 0) {
    Array<BoutReal> x0(dim);
BOUT_OMP(parallel default(shared))
    for(int num =0;num < 2;num++) {
BOUT_OMP(for)
      for(int i = 0;i<dim;i++) x0[i] = x[i];    

      int xend = lnx[level]+1;
      int zend = lnz[level]+1;
BOUT_OMP(for collapse(2))
      for(int i=1;i<xend;i++)
        for(int k=1;k<zend;k++) {
          int nn = i*mm+k;
          BoutReal val = b[nn] - matmg[level][nn*9+3]*x0[nn-1]
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
  }
  else {
    for(int i = 1;i<lnx[level]+1;i++)
      for(int k=1;k<lnz[level]+1;k++) {
        int nn = i*mm+k;
        BoutReal val = b[nn] - matmg[level][nn*9+3]*x[nn-1]
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
        int nn = i*mm+k;
        BoutReal val = b[nn] - matmg[level][nn*9+3]*x[nn-1]
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
  int it,etest = 1,MAXIT;
  BoutReal ini_e,error,a0,a1,rederr,perror;
  BoutReal **v;
  BoutReal c[MAXGM+1],s[MAXGM+1],y[MAXGM+1],g[MAXGM+1],h[MAXGM+1][MAXGM+1];

  if((level == 0) || (iplag == 0)) MAXIT = 40000;
  else MAXIT = 500;

  int ldim = (lnx[level]+2)*(lnz[level]+2);
  // Could we use a Matrix here?
  v = new BoutReal *[MAXGM+1];
  for(int i=0;i<MAXGM+1;i++) v[i] = new BoutReal[ldim];

  Array<BoutReal> p(ldim);
  Array<BoutReal> q(ldim);
  Array<BoutReal> r(ldim);

  BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
  for(int i = 0;i<ldim;i++) sol[i] = 0.0;
  int num = 0;

  ini_e = sqrt(vectorProd(level,rhs,rhs));
  if((pcheck == 1) && (rProcI == 0)) {
    output<<numP<<"--In GMRES ini "<<ini_e<<endl;
  }
  if(ini_e < atol*rtol) {
    if((pcheck == 1) && (rProcI == 0)) {
      output<<numP<<"Don't need to solve. E= "<<ini_e<<endl;
    }
    // Clean up memory before returning from method
    for(int i=0;i<MAXGM+1;i++) {
      delete [] v[i];
    }
    delete [] v;
    return;
  }
BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
  for(int i = 0;i<ldim;i++) r[i] = 0.0;
  if (iplag == 0)
    smoothings(level, std::begin(r), rhs);
  else
    cycleMG(level, std::begin(r), rhs);
  BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
  for(int i = 0;i < ldim;i++) v[0][i] = r[i];
  perror = ini_e;
  do{
    a1 = vectorProd(level,v[0],v[0]);
    a1 = sqrt(a1);
    if(fabs(a1) < atol*rtol) {
      output<<num<<" First a1 in GMRES is wrong at level "<<level<<": "<<a1<<endl;
    }
    a0 = 1.0/a1;
    g[0] = a1;
BOUT_OMP(parallel default(shared))
    {
BOUT_OMP(for)
      for(int i=0;i<ldim;i++) v[0][i] *= a0;
BOUT_OMP(for)
      for(int i=1;i<MAXGM+1;i++) g[i] = 0.0;
    }
    for(it = 0;it<MAXGM;it++) {
      multiAVec(level, v[it], std::begin(q));
      BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
      for(int i=0;i<ldim;i++) v[it+1][i] = 0.0;

      if (iplag == 0)
        smoothings(level, v[it + 1], std::begin(q));
      else
        cycleMG(level, v[it + 1], std::begin(q));

      for(int i=0;i<it+1;i++) h[i][it] = vectorProd(level,v[it+1],v[i]);
      for(int i=0;i<it+1;i++) {
        a0 = -h[i][it];
        for(int k=0;k<ldim;k++) v[it+1][k] += a0*v[i][k];
      }
      a1 = vectorProd(level,v[it+1],v[it+1]);
      a1 = sqrt(a1);

      // if ldim==9 then there is only one grid point at this level, so the
      // solution will be exact, the residual will vanish and we will exit this
      // loop on the first iteration, so the value of a0=1/a1=infinity will
      // never be used. Therefore this check is not needed in that case
      if(fabs(a1) < atol*rtol && ldim > 9) {
        output<<num<<" Second a1 in GMRES is wrong at level "<<level<<": "<<a1<<endl;
      }
      a0 = 1.0/a1;
      h[it+1][it] = a1;
BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
      for(int i=0;i<ldim;i++) v[it+1][i] *= a0;

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
BOUT_OMP(parallel default(shared))
      {
BOUT_OMP(for)
        for(int i=0;i<ldim;i++) p[i] = sol[i];
BOUT_OMP(for)
        for(int k=0;k<ldim;k++)
          for(int i=0;i<=it;i++)
            p[k] += y[i]*v[i][k]; 
      }

      /* Get r_m and test convergence.*/
      residualVec(level, std::begin(p), rhs, std::begin(r));
      error = sqrt(vectorProd(level, std::begin(r), std::begin(r)));
      num += 1;
      if(error > dtol)
        throw BoutException("GMRES reached dtol with error %16.10f at iteration %d\n",error,num);
      if(num > MAXIT)
        throw BoutException("GMRES reached MAXIT with error %16.10f at iteration %d\n",error,num);
      if(error <= rtol*ini_e+atol) {
        etest = 0;
        break;
      }
      // J. Omotani, 27/2/2018: I think this test is intended to check for slow
      // convergence of the GMRES solve, and 'abort' if it is converging
      // slowly. This is OK on a coarse level solver, because at worst it means
      // the top-level iteration will have to continue but the top level
      // iteration should only be stopped by the previous test against the
      // tolerance.
      // Therefore, check that this is not the top-level solver before applying
      // this test.
      if( (level < mglevel-1) && (fabs(perror-error)/error < rtol) ) {
        if(it == 0) etest = 0;
        num -= 1;
        break;
      }
      perror = error;
    }
    /* Restart with new initial */
BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
    for(int i = 0;i<ldim;i++) v[0][i] = 0.0;
    if (iplag == 0)
      smoothings(level, v[0], std::begin(r));
    else
      cycleMG(level, v[0], std::begin(r));

    BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
    for(int i = 0;i<ldim;i++) sol[i] = p[i];
    if(num>MAXIT)
      throw BoutException(" GMRES Iteration limit.\n");
    //    if((etest == 1) & (xProcI == 0)) 
    //  printf("Restart GMRES  %d | %20.14f\n",num,error/ini_e);    
  } while(etest == 1); 

  if(level == 0) {
    if((rProcI == 0) && (pcheck == 1)) {
      rederr = log(error / ini_e) / static_cast<BoutReal>(num);
      rederr = exp(rederr);      
      printf("The average error reduction of GMRES %d: %14.8f(%18.10f)\n",num,rederr,error);
      fflush(stdout);
    }
  }
  else {
    if((rProcI == 0) && (pcheck == 1)) {
      rederr = log(error / ini_e) / static_cast<BoutReal>(num);
      rederr = exp(rederr);      
      printf("The average error reduction of PGMRES %d: %14.8f(%18.10f)\n",num,rederr,error);
      fflush(stdout);
    }
  }
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
BOUT_OMP(parallel default(shared) )
  {
    int xend = lnx[level]+1;
    int zend = lnz[level]+1;
BOUT_OMP(for reduction(+:ini_e) collapse(2))
    for(int i= 1;i<xend;i++){
      for(int k=1;k<zend;k++) {
        int ii = i*(lnz[level]+2)+k;
        ini_e += x[ii]*y[ii];
      }
    }
  }
  if(numP > 1) 
    MPI_Allreduce(&ini_e,&val,1,MPI_DOUBLE,MPI_SUM,commMG);
  else val = ini_e;

  return(val);  
}

void MultigridAlg::multiAVec(int level, BoutReal *x, BoutReal *b) {

  int mm = lnz[level]+2;
BOUT_OMP(parallel default(shared))
  {
BOUT_OMP(for)
    for(int i = 0;i<mm*(lnx[level]+2);i++) b[i] = 0.0;

    int xend = lnx[level]+1;
    int zend = lnz[level]+1;
BOUT_OMP(for collapse(2))
    for(int i=1;i<xend;i++) {
      for(int k=1;k<zend;k++) {
        int nn = i*mm+k;
        b[nn] = matmg[level][nn*9+4]*x[nn] + matmg[level][nn*9+3]*x[nn-1]
          +matmg[level][nn*9+5]*x[nn+1] + matmg[level][nn*9+1]*x[nn-mm]
          +matmg[level][nn*9+7]*x[nn+mm] +matmg[level][nn*9]*x[nn-mm-1]
          +matmg[level][nn*9+2]*x[nn-mm+1] + matmg[level][nn*9+6]*x[nn+mm-1]
          +matmg[level][nn*9+8]*x[nn+mm+1];
      } 
    }
  }
  communications(b,level);
}

void MultigridAlg::residualVec(int level, BoutReal *x, BoutReal *b,
BoutReal *r) {

  int mm;
  mm = lnz[level]+2;
BOUT_OMP(parallel default(shared))
  {
BOUT_OMP(for)
    for(int i = 0;i<mm*(lnx[level]+2);i++) r[i] = 0.0;

    int xend = lnx[level]+1;
    int zend = lnz[level]+1;
BOUT_OMP(for collapse(2))
    for(int i=1;i<xend;i++) {
      for(int k=1;k<zend;k++) {
        int nn = i*mm+k;
        BoutReal val = matmg[level][nn*9+4]*x[nn] + matmg[level][nn*9+3]*x[nn-1]
          +matmg[level][nn*9+5]*x[nn+1] + matmg[level][nn*9+1]*x[nn-mm]
          +matmg[level][nn*9+7]*x[nn+mm] +matmg[level][nn*9]*x[nn-mm-1]
          +matmg[level][nn*9+2]*x[nn-mm+1] + matmg[level][nn*9+6]*x[nn+mm-1]
          +matmg[level][nn*9+8]*x[nn+mm+1];
        r[nn] = b[nn]-val;
      } 
    }
  }
  communications(r,level);

}

void MultigridAlg::setMatrixC(int level) {

  BoutReal ratio = 8.0; 

BOUT_OMP(parallel default(shared))
  {
BOUT_OMP(for)
    for(int i=0;i<(lnx[level-1]+2)*(lnz[level-1]+2)*9;i++)
      matmg[level-1][i] = 0.0;

    int xend = lnx[level-1]+1;
    int zend = lnz[level-1]+1;
BOUT_OMP(for collapse(2))
    for(int i=1;i<xend;i++) {
      for(int k = 1;k<zend;k++) {
        int i2 = 2*i-1;
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

}

void MultigridAlg::communications(BoutReal* x, int level) {
  // Note: currently z-direction communications and x-direction communications are done
  // independently. They should really be done together if zNP>1 and xNP>1 (with a single
  // MPI_Waitall after all the MPI_Isend and MPI_Irecv calls, instead of the two
  // MPI_Waitalls there are now. As there are never any z-communications at the moment, it
  // is not worth implementing for now.

  MPI_Status status[4];
  int stag, rtag;
  MAYBE_UNUSED(int ierr);

  if(zNP > 1) {
    MPI_Request requests[] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};
    MPI_Datatype xvector;

    // Create datatype to hold z-direction guard cells, which are not continuous in memory
    ierr = MPI_Type_vector(lnx[level], 1, lnz[level]+2, MPI_DOUBLE, &xvector);
    ASSERT1(ierr == MPI_SUCCESS);
    
    ierr = MPI_Type_commit(&xvector);
    ASSERT1(ierr == MPI_SUCCESS);
    
    // Receive from z-
    rtag = zProcM;
    ierr = MPI_Irecv(&x[lnz[level]+2], 1, xvector, zProcM, rtag, commMG, &requests[2]);
    ASSERT1(ierr == MPI_SUCCESS);

    // Receive from z+
    rtag = zProcP+numP;
    ierr = MPI_Irecv(&x[2*(lnz[level]+2)-1], 1, xvector, zProcP, rtag, commMG,
        &requests[3]);
    ASSERT1(ierr == MPI_SUCCESS);

    // Send to z+
    stag = rProcI;
    ierr = MPI_Isend(&x[2*(lnz[level]+2)-2], 1, xvector, zProcP, stag, commMG,
        &requests[0]);
    ASSERT1(ierr == MPI_SUCCESS);

    // Send to z-
    stag = rProcI+numP;
    ierr = MPI_Isend(&x[lnz[level]+3], 1, xvector, zProcM, stag, commMG, &requests[1]);
    ASSERT1(ierr == MPI_SUCCESS);

    // Wait for communications to complete
    ierr = MPI_Waitall(4, requests, status);
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
    // Note: periodic x-direction not handled here

    MPI_Request requests[] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};

    if (xProcI > 0) {
      // Receive from x-
      rtag = xProcM;
      ierr = MPI_Irecv(&x[0], lnz[level]+2, MPI_DOUBLE, xProcM, rtag, commMG, &requests[2]);
      ASSERT1(ierr == MPI_SUCCESS);
    }
    
    if (xProcI < xNP - 1) {
      // Receive from x+
      rtag = xProcP+xNP;;
      ierr = MPI_Irecv(&x[(lnx[level]+1)*(lnz[level]+2)], lnz[level]+2, MPI_DOUBLE, xProcP,
          rtag, commMG, &requests[3]);
      ASSERT1(ierr == MPI_SUCCESS);

      // Send to x+
      stag = rProcI;
      ierr = MPI_Isend(&x[lnx[level]*(lnz[level]+2)], lnz[level]+2, MPI_DOUBLE, xProcP,
          stag, commMG, &requests[0]);
      ASSERT1(ierr == MPI_SUCCESS);
    }

    if (xProcI > 0) {
      // Send to x-
      stag = rProcI+xNP;
      ierr = MPI_Isend(&x[lnz[level]+2], lnz[level]+2, MPI_DOUBLE, xProcM, stag, commMG,
          &requests[1]);
      ASSERT1(ierr == MPI_SUCCESS);
    }

    // Wait for communications to complete
    ierr = MPI_Waitall(4, requests, status);
    ASSERT1(ierr == MPI_SUCCESS);
  } else {
    for (int i=0;i<lnz[level]+2;i++) {
      x[i] = x[lnx[level]*(lnz[level]+2)+i];
      x[(lnx[level]+1)*(lnz[level]+2)+i] = x[(lnz[level]+2)+i];
    }
  }
}


void MultigridAlg::solveMG(BoutReal *sol,BoutReal *rhs,int level) {
  int m,MAXIT = 150;
  BoutReal ini_e,perror,error,rederr;
  int ldim = (lnx[level]+2)*(lnz[level]+2);

BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
  for(int i = 0;i<ldim;i++) sol[i] = 0.0;

  ini_e = vectorProd(level,rhs,rhs);
  if(ini_e < 0.0)
    throw BoutException("In MG Initial Error %10.4e \n",ini_e);
  ini_e = sqrt(ini_e);
  if((pcheck == 1) && (rProcI == 0)) 
    printf("%d \n  In MGsolve ini = %24.18f\n",numP,ini_e);
  Array<BoutReal> y(ldim);
  Array<BoutReal> r(ldim);
  BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
  for(int i = 0;i<ldim;i++) r[i] = rhs[i];

  perror = ini_e;
  for(m=0;m<MAXIT;m++) {
BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
    for(int i = 0;i<ldim;i++) y[i] = 0.0;
    cycleMG(level, std::begin(y), std::begin(r));
    BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
    for(int i = 0;i<ldim;i++) sol[i] = sol[i]+y[i];
    residualVec(level, sol, rhs, std::begin(r));
    error = sqrt(vectorProd(level, std::begin(r), std::begin(r)));
    if((pcheck == 1) && (rProcI == 0)) 
      printf("%d \n  In MGsolve error = %24.18f\n",m,error);
    if(error < rtol*ini_e+atol) break;
    if((fabs(perror-error)/error <rtol) || (error > dtol))
      throw BoutException("In MG Limited Error %10.4e \n",error);
    perror = error;
  }

  if((rProcI == 0) && (pcheck == 1)) {
    rederr = log(error / ini_e) / (static_cast<BoutReal>(m) + 1.0);
    rederr = exp(rederr); 
    if(m == MAXIT) 
      printf("Reached maximum iteration: %14.8f\n",error);
    printf("The average error reduction of MG %d: %14.8f(%18.12f)\n",m+1,rederr,error);
    fflush(stdout);
  }

  return;
}
