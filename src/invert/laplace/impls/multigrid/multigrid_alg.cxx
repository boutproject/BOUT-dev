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

void MultigridAlg::initializeVectors() {
  // Set up working vectors
  r_array.reallocate(mglevel);
  pr_array.reallocate(mglevel);
  y_array.reallocate(mglevel);
  iy_array.reallocate(mglevel);
  for (int level=0; level<mglevel; ++level) {
    r_array[level] = bout::utils::make_unique<MultigridVector>(*this, level);
    pr_array[level] = bout::utils::make_unique<MultigridVector>(*this, level);
    y_array[level] = bout::utils::make_unique<MultigridVector>(*this, level);
    iy_array[level] = bout::utils::make_unique<MultigridVector>(*this, level);
  }

  int gmres_max_mglevel = 0;
  if (mgplag) {
    gmres_max_mglevel = mglevel;
  }
  // otherwise pGMRES is not used, so don't need to initialize these arrays
  v_gmres.reallocate(gmres_max_mglevel);
}

MultigridAlg::~MultigridAlg() {
  output<<"End deconstruction Malg AAAA "<<numP<<endl;
  for(int i = 0;i<mglevel;i++) delete [] matmg[i];
  delete [] matmg;
}

// Get array of MAXGM MultigridVector objects to be used by pGMRES.
// Don't need these for all levels, so only create on demand.
Array<std::unique_ptr<MultigridVector>>& MultigridAlg::get_v_gmres(int level) {
  if (v_gmres[level].empty()) {
    v_gmres[level].reallocate(MAXGM + 1);
    for (int i = 0; i < MAXGM + 1; ++i) {
      v_gmres[level][i] = bout::utils::make_unique<MultigridVector>(*this, level);
    }
  }

  return v_gmres[level];
}

void MultigridAlg::getSolution(MultigridVector& x, MultigridVector& b, int flag) {

  // swap ghost cells of initial guess
  x.communicate();
  // don't think the ghost cells of the rhs are used, so don't need to be communicated
  // /JTO 17-4-2019

  if(flag == 0) {
    //Solve exactly
    if(mglevel == 1) pGMRES(x,b,mglevel-1,1);
    else if(mgplag == 1) pGMRES(x,b,mglevel-1,1);
    else solveMG(x,b,mglevel-1);
  }
  else {
    cycleMG(mglevel-1,x,b);
    if(flag > 1) {
      int level = mglevel-1;
      int ldim = (lnx[level]+2)*(lnz[level]+2);
      // This is not very efficient, because we don't re-use the setup for y and r, but
      // this branch should not generally be used since it just does two V-cycles, rather
      // than checking a solution has been found
      MultigridVector y(*this, level);
      MultigridVector r(*this, level);
      for(int n = 1;n<flag;n++) {
        residualVec(level, x, b, r);
        BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
        for(int i = 0;i<ldim;i++) y[i] = 0.0;
        cycleMG(level, y, r);
        BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
        for(int i = 0;i<ldim;i++) x[i] = x[i]+y[i];
      }
    }
  }
}


void MultigridAlg::cycleMG(int level, MultigridVector& sol, MultigridVector& rhs)
{
  if(level == 0) {
    lowestSolver(sol,rhs,0);
  }
  else {
    smoothings(level,sol,rhs);

    MultigridVector& r = *r_array[level];
    MultigridVector& pr = *pr_array[level];
    MultigridVector& y = *y_array[level];
    MultigridVector& iy = *iy_array[level];

    residualVec(level, sol, rhs, r);

    projection(level, r, pr);

    BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
    for(int i=0;i<(lnx[level-1]+2)*(lnz[level-1]+2);i++) y[i] = 0.0;

    cycleMG(level - 1, y, pr);

    prolongation(level - 1, y, iy);
    BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
    for(int i=0;i<(lnx[level]+2)*(lnz[level]+2);i++) 
       sol[i] += iy[i];

    smoothings(level, sol, rhs);
  }
  return;
}

void MultigridAlg::projection(int level, MultigridVector& r, MultigridVector& pr)
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
  pr.communicate();
  return;
}

void MultigridAlg::prolongation(int level, MultigridVector& x,MultigridVector& ix) {

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
  ix.communicate();
  return;
}

void MultigridAlg::smoothings(int level, MultigridVector& x, MultigridVector& b) {

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
      x.communicate();
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
    x.communicate();
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
    x.communicate();
  }
  return;
}

void MultigridAlg::pGMRES(MultigridVector& sol,MultigridVector& rhs, int level, int iplag) {
  int it,etest = 1,MAXIT;
  BoutReal ini_e,error,a0,a1,rederr,perror;
  BoutReal c[MAXGM+1],s[MAXGM+1],y[MAXGM+1],g[MAXGM+1],h[MAXGM+1][MAXGM+1];
  MultigridVector& p = *pr_array[level];
  MultigridVector& q = *y_array[level];
  MultigridVector& r = *r_array[level];
  auto& v = get_v_gmres(level);

  if((level == 0) || (iplag == 0)) MAXIT = 40000;
  else MAXIT = 500;

  int ldim = (lnx[level]+2)*(lnz[level]+2);

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
    return;
  }
BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
  for(int i = 0;i<ldim;i++) r[i] = 0.0;
  if (iplag == 0)
    smoothings(level, r, rhs);
  else
    cycleMG(level, r, rhs);
  BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
  for(int i = 0;i < ldim;i++) (*v[0])[i] = r[i];
  perror = ini_e;
  do{
    a1 = vectorProd(level,*v[0],*v[0]);
    a1 = sqrt(a1);
    if(fabs(a1) < atol*rtol) {
      output<<num<<" First a1 in GMRES is wrong at level "<<level<<": "<<a1<<endl;
    }
    a0 = 1.0/a1;
    g[0] = a1;
BOUT_OMP(parallel default(shared))
    {
BOUT_OMP(for)
      for(int i=0;i<ldim;i++) (*v[0])[i] *= a0;
BOUT_OMP(for)
      for(int i=1;i<MAXGM+1;i++) g[i] = 0.0;
    }
    for(it = 0;it<MAXGM;it++) {
      multiAVec(level, *v[it], q);
      BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
      for(int i=0;i<ldim;i++) (*v[it+1])[i] = 0.0;

      if (iplag == 0)
        smoothings(level, *v[it + 1], q);
      else
        cycleMG(level, *v[it + 1], q);

      for(int i=0;i<it+1;i++) h[i][it] = vectorProd(level,*v[it+1],*v[i]);
      for(int i=0;i<it+1;i++) {
        a0 = -h[i][it];
        for(int k=0;k<ldim;k++) (*v[it+1])[k] += a0*(*v[i])[k];
      }
      a1 = vectorProd(level,*v[it+1],*v[it+1]);
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
      for(int i=0;i<ldim;i++) (*v[it+1])[i] *= a0;

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
            p[k] += y[i]*(*v[i])[k];
      }

      /* Get r_m and test convergence.*/
      residualVec(level, p, rhs, r);
      error = sqrt(vectorProd(level, r, r));
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
    for(int i = 0;i<ldim;i++) (*v[0])[i] = 0.0;
    if (iplag == 0)
      smoothings(level, *v[0], r);
    else
      cycleMG(level, *v[0], r);

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

void MultigridAlg::lowestSolver(MultigridVector& x, MultigridVector& b, int UNUSED(plag)) {
  pGMRES(x, b, 0, 0);
}

BoutReal MultigridAlg::vectorProd(int level, MultigridVector& x, MultigridVector& y) {
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

void MultigridAlg::multiAVec(int level, MultigridVector& x, MultigridVector& b) {

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
  b.communicate();
}

void MultigridAlg::residualVec(int level, MultigridVector& x, MultigridVector& b,
    MultigridVector& r) {

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
  r.communicate();
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

void MultigridAlg::solveMG(MultigridVector& sol, MultigridVector& rhs, int level) {
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
  MultigridVector& y = *y_array[level];
  MultigridVector& r = *r_array[level];
  BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
  for(int i = 0;i<ldim;i++) r[i] = rhs[i];

  perror = ini_e;
  for(m=0;m<MAXIT;m++) {
BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
    for(int i = 0;i<ldim;i++) y[i] = 0.0;
    cycleMG(level, y, r);
    BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
    for(int i = 0;i<ldim;i++) sol[i] = sol[i]+y[i];
    residualVec(level, sol, rhs, r);
    error = sqrt(vectorProd(level, r, r));
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

MultigridVector::MultigridVector(MultigridAlg& mgAlg_in, int level)
  : lnx(mgAlg_in.lnx[level]), lnz(mgAlg_in.lnz[level]), mgAlg(mgAlg_in) {

  // Create data
  data.reallocate((lnx + 2)*(lnz + 2));

  // Set up persistent communications
  MAYBE_UNUSED(int ierr);
  int rtag, stag;
  int numP = mgAlg.xNP*mgAlg.zNP;
  BoutReal* x = std::begin(data);

  // Count instances of this class, used to make sure tags don't conflict
  static int counter = -1;
  ++counter;

  if (mgAlg.zNP > 1) {
    MPI_Datatype xvector;

    ierr = MPI_Type_vector(lnx, 1, lnz+2, MPI_DOUBLE, &xvector);
    ASSERT0(ierr == MPI_SUCCESS);

    ierr = MPI_Type_commit(&xvector);
    ASSERT1(ierr == MPI_SUCCESS);

    // different tags for each receive on each level
    rtag = 4*counter*numP + mgAlg.zProcM;
    ierr = MPI_Recv_init(&x[lnz+2], 1, xvector, mgAlg.zProcM, rtag, mgAlg.commMG,
        &(zRequests[0]));
    ASSERT1(ierr == MPI_SUCCESS);

    // Receive from z+
    rtag = (4*counter+1)*numP + mgAlg.zProcP;
    ierr = MPI_Recv_init(&x[2*(lnz+2)-1], 1, xvector, mgAlg.zProcP, rtag, mgAlg.commMG,
        &(zRequests[1]));
    ASSERT1(ierr == MPI_SUCCESS);

    // Send to z+
    stag = 4*counter*numP + mgAlg.rProcI;
    ierr = MPI_Send_init(&x[2*(lnz+2)-2], 1, xvector, mgAlg.zProcP, stag, mgAlg.commMG,
        &(zRequests[2]));
    ASSERT1(ierr == MPI_SUCCESS);

    // Send to z-
    stag = (4*counter+1)*numP + mgAlg.rProcI;
    ierr = MPI_Send_init(&x[lnz+3], 1, xvector, mgAlg.zProcM, stag, mgAlg.commMG,
        &(zRequests[3]));
    ASSERT1(ierr == MPI_SUCCESS);

    ierr = MPI_Type_free(&xvector);
    ASSERT1(ierr == MPI_SUCCESS);
  }
  if (mgAlg.xNP > 1) {
    // Note: periodic x-direction not handled here
    if (mgAlg.xProcI > 0) {
      // Receive from x-
      rtag = (4*counter+2)*numP + mgAlg.xProcM;
      ierr = MPI_Irecv(&x[0], lnz+2, MPI_DOUBLE, mgAlg.xProcM, rtag, mgAlg.commMG,
          &(xRequests[0]));
      ASSERT1(ierr == MPI_SUCCESS);

      // Send to x-
      stag = (4*counter+3)*numP + mgAlg.rProcI;
      ierr = MPI_Isend(&x[lnz+2], lnz+2, MPI_DOUBLE, mgAlg.xProcM, stag, mgAlg.commMG,
          &(xRequests[2]));
      ASSERT1(ierr == MPI_SUCCESS);
    } else {
      xRequests[0] = MPI_REQUEST_NULL;
      xRequests[2] = MPI_REQUEST_NULL;
    }

    if (mgAlg.xProcI < mgAlg.xNP - 1) {
      // Receive from x+
      rtag = (4*counter+3)*numP + mgAlg.xProcP;
      ierr = MPI_Irecv(&x[(lnx+1)*(lnz+2)], lnz+2, MPI_DOUBLE, mgAlg.xProcP, rtag,
          mgAlg.commMG, &(xRequests[1]));
      ASSERT1(ierr == MPI_SUCCESS);

      // Send to x+
      stag = (4*counter+2) + mgAlg.rProcI;
      ierr = MPI_Isend(&x[lnx*(lnz+2)], lnz+2, MPI_DOUBLE, mgAlg.xProcP, stag,
          mgAlg.commMG, &(xRequests[3]));
      ASSERT1(ierr == MPI_SUCCESS);
    } else {
      xRequests[1] = MPI_REQUEST_NULL;
      xRequests[3] = MPI_REQUEST_NULL;
    }
  }
}

MultigridVector::~MultigridVector() {
  // Free persistent communications
  if (mgAlg.zNP > 1) {
    MPI_Request_free(&zRequests[0]);
    MPI_Request_free(&zRequests[1]);
    MPI_Request_free(&zRequests[2]);
    MPI_Request_free(&zRequests[3]);
  }
  if (mgAlg.xNP > 1) {
    MPI_Request_free(&xRequests[0]);
    MPI_Request_free(&xRequests[1]);
    MPI_Request_free(&xRequests[2]);
    MPI_Request_free(&xRequests[3]);
  }
}

void MultigridVector::communicate() {
  // Note: currently z-direction communications and x-direction communications are done
  // independently. They should really be done together if zNP>1 and xNP>1 (with a single
  // MPI_Waitall after all the MPI_Isend and MPI_Irecv calls, instead of the two
  // MPI_Waitalls there are now. As there are never any z-communications at the moment, it
  // is not worth implementing for now.

  MPI_Status status[4];
  MAYBE_UNUSED(int ierr);

  if(mgAlg.zNP > 1) {
    // post receives first
    ierr = MPI_Startall(2, zRequests);
    ASSERT1(ierr == MPI_SUCCESS);

    // start sends
    ierr = MPI_Startall(2, &(zRequests[2]));
    ASSERT1(ierr == MPI_SUCCESS);

    // Wait for communications to complete
    ierr = MPI_Waitall(4, zRequests, status);
    ASSERT1(ierr == MPI_SUCCESS);
  } else {
    for (int i=1;i<lnx+1;i++) {
      data[i*(lnz+2)] = data[(i+1)*(lnz+2)-2];
      data[(i+1)*(lnz+2)-1] = data[i*(lnz+2)+1];
    }
  }
  if (mgAlg.xNP > 1) {
    // post receives first
    ierr = MPI_Startall(2, xRequests);
    ASSERT1(ierr == MPI_SUCCESS);

    // start sends
    ierr = MPI_Startall(2, &(xRequests[2]));
    ASSERT1(ierr == MPI_SUCCESS);

    // Wait for communications to complete
    ierr = MPI_Waitall(4, xRequests, status);
    ASSERT1(ierr == MPI_SUCCESS);
  } else {
    for (int i=0;i<lnz+2;i++) {
      data[i] = data[lnx*(lnz+2)+i];
      data[(lnx+1)*(lnz+2)+i] = data[(lnz+2)+i];
    }
  }
}
