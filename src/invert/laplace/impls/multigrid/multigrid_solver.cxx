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

#include "bout/build_defines.hxx"

#if not BOUT_USE_METRIC_3D

#include "multigrid_laplace.hxx"
#include "bout/unused.hxx"
#include <bout/openmpwrap.hxx>

Multigrid1DP::Multigrid1DP(int level, int lx, int lz, int gx, int dl, int merge,
                           MPI_Comm comm, int check)
    : MultigridAlg(level, lx, lz, gx, lz, comm, check) {

  mglevel = level;

  if (pcheck > 0) {
    output << "Construct MG1DP " << endl;
  }
  commMG = comm;
  MPI_Comm_size(commMG, &xNP);
  MPI_Comm_rank(commMG, &rProcI);
  xProcI = rProcI;
  if (xNP > 1) {
    if (xProcI == 0) {
      xProcM = xNP - 1;
    } else {
      xProcM = xProcI - 1;
    }
    if (xProcI == xNP - 1) {
      xProcP = 0;
    } else {
      xProcP = xProcI + 1;
    }
  } else {
    xProcI = 0;
    xProcM = 0;
    xProcP = 0;
  }
  zNP = 1;
  numP = xNP;
  zProcI = xProcI;
  zProcP = xProcI;
  zProcM = xProcI;

  if (pcheck == 1) {
    output << "In MG1DP level " << mglevel << " xNP=" << xNP << "(" << xProcI << ")"
           << endl;
    output << "lest level is " << dl << "(" << numP << ")" << endl;
  }

  int nz, kk, nx;
  if (dl > 0) {
    // Find levels for more coarser spaces
    if (numP > merge) {
      int nn = numP;
      int mm = static_cast<int>(sqrt(numP));
      kk = 1;
      for (int n = nn; n > 1; n--) {
        if (nn % 2 != 0) {
          n = 1;
        } else {
          kk += 1;
          nn = nn / 2;
        }
      }
      nn = lnz[0];
      int kz = 1;
      for (int n = nn; n > 1; n--) {
        if (nn % 2 != 0) {
          n = 1;
        } else {
          kz += 1;
          nn = nn / 2;
        }
      }

      if (kz < kk) {
        kk = kz;
      }
      nz = 1;
      nx = xNP;
      for (int n = 0; n < kk; n++) {
        if (nz * 2 <= mm) {
          nz = 2 * nz;
          nx = nx / 2;
        } else {
          n = kk;
        }
      }

      lx = gnx[0] / nx;
      lz = lnz[0] / nz;
      // 2D parallelization
      kk = 1;
      int llx = lx;
      int llz = lz;
      for (int n = dl; n > 0; n--) {
        if ((llx % 2 == 0) && (llz % 2 == 0)) {
          kk += 1;
          llx = llx / 2;
          llz = llz / 2;
        } else {
          n = 1;
        }
      }
      if (kk > 0) {
        kflag = 1;
      } else {
        kflag = 2;
      }
    } else {
      kflag = 2;
    }
    if (kflag == 1) {
      if (pcheck == 1) {
        output << "To MG2DP " << kk << "xNP=" << nx << "(" << nz << ")" << endl;
        output << "lest level is " << dl - kk + 1 << "(" << lx << ", " << lz << ")"
               << endl;
      }
      int colors = rProcI / nz;
      int keys = rProcI / nz;
      MPI_Comm_split(commMG, colors, keys, &comm2D);
      rMG = bout::utils::make_unique<Multigrid2DPf1D>(
          kk, lx, lz, gnx[0], lnz[0], dl - kk + 1, nx, nz, commMG, pcheck);
    } else {
      int nn = gnx[0];
      int mm = gnz[0];
      int kk = 1;
      for (int n = dl; n > 0; n--) {
        if ((nn % 2 == 0) && (mm % 2 == 0)) {
          kk += 1;
          nn = nn / 2;
          mm = mm / 2;
        } else {
          n = 1;
        }
      }
      if (pcheck == 1) {
        output << "To Ser " << kk << " xNP=" << xNP << "(" << zNP << ")" << endl;
        output << kflag << " total dim " << gnx[0] << "(" << lnz[0] << ")" << endl;
      }
      sMG = bout::utils::make_unique<MultigridSerial>(kk, gnx[0], lnz[0], commMG, pcheck);
    }
  } else {
    kflag = 0;
  }
}

void Multigrid1DP::setMultigridC(int UNUSED(plag)) {

  int level = mglevel - 1;
  for (int n = level; n > 0; n--) {
    if (pcheck == 2) {
      output << n << "matrix in 1DP = " << lnx[n - 1] << "," << lnz[n - 1] << endl;
      output << gnx[n - 1] << "," << gnz[n - 1] << endl;
    }
    setMatrixC(n);
  }

  if (kflag == 1) {
    level = rMG->mglevel - 1;
    convertMatrixF2D(level);

    if (level > 0) {
      rMG->setMultigridC(0);
    }

    if (pcheck == 2) {
      for (int i = level; i >= 0; i--) {
        FILE* outf;
        char outfile[256];
        sprintf(outfile, "2DP_matC%1d_%d.mat", i, rMG->rProcI);
        output << "Out file= " << outfile << endl;
        outf = fopen(outfile, "w");
        int dim = (rMG->lnx[i] + 2) * (rMG->lnz[i] + 2);
        fprintf(outf, "dim = %d (%d, %d)\n", dim, rMG->lnx[i], rMG->lnz[i]);

        for (int ii = 0; ii < dim; ii++) {
          fprintf(outf, "%d ==", ii);
          for (int j = 0; j < 9; j++) {
            fprintf(outf, "%12.6f,", rMG->matmg[i][ii * 9 + j]);
          }
          fprintf(outf, "\n");
        }
        fclose(outf);
      }
    }
  } else if (kflag == 2) {
    level = sMG->mglevel - 1;
    convertMatrixFS(level);

    if (level > 0) {
      sMG->setMultigridC(0);
    }
    if (pcheck == 3) {
      for (int i = level; i >= 0; i--) {
        FILE* outf;
        char outfile[256];
        sprintf(outfile, "S1D_matC%1d_%d.mat", i, sMG->rProcI);
        output << "Out file= " << outfile << endl;
        outf = fopen(outfile, "w");
        int dim = (sMG->lnx[i] + 2) * (sMG->lnz[i] + 2);
        fprintf(outf, "dim = %d\n", dim);

        for (int ii = 0; ii < dim; ii++) {
          fprintf(outf, "%d ==", ii);
          for (int j = 0; j < 9; j++) {
            fprintf(outf, "%12.6f,", sMG->matmg[i][ii * 9 + j]);
          }
          fprintf(outf, "\n");
        }
        fclose(outf);
      }
    }
  }
}

void Multigrid1DP::setValueS() {

  if (kflag == 1) {
    rMG->mgplag = mgplag;
    rMG->mgsm = mgsm;
    rMG->cftype = cftype;
    rMG->rtol = rtol;
    rMG->atol = atol;
    rMG->dtol = dtol;
    rMG->omega = omega;
    rMG->setValueS();
  } else if (kflag == 2) {
    sMG->mgplag = mgplag;
    sMG->mgsm = mgsm;
    sMG->cftype = cftype;
    sMG->rtol = rtol;
    sMG->atol = atol;
    sMG->dtol = dtol;
    sMG->omega = omega;
  }
}

void Multigrid1DP::setPcheck(int check) {

  pcheck = check;
  if (kflag == 1) {
    rMG->setPcheck(check);
  } else if (kflag == 2) {
    sMG->pcheck = check;
  }
}

void Multigrid1DP::lowestSolver(BoutReal* x, BoutReal* b, int UNUSED(plag)) {

  if (kflag == 1) {
    int level = rMG->mglevel - 1;
    int dim = (rMG->lnx[level] + 2) * (rMG->lnz[level] + 2);
    Array<BoutReal> y(dim);
    Array<BoutReal> r(dim);

    int ggx = rMG->lnx[level];
    int dimg = (ggx + 2) * (gnz[0] + 2);
    Array<BoutReal> yl(dimg);
    Array<BoutReal> yg(dimg);

    int nx = (xProcI % rMG->zNP) * lnx[0];

    BOUT_OMP_PERF(parallel default(shared))
    {
      BOUT_OMP_PERF(for)
      for (int i = 0; i < dim; i++) {
        y[i] = 0.0;
        r[i] = 0.0;
      }

      BOUT_OMP_PERF(for)
      for (int i = 0; i < dimg; i++) {
        yl[i] = 0.0;
        yg[i] = 0.0;
      }

      int xend = lnx[0] + 1;
      int zend = lnz[0] + 1;
      BOUT_OMP_PERF(for collapse(2))
      for (int ix = 1; ix < xend; ix++) {
        for (int iz = 1; iz < zend; iz++) {
          int nn = (nx + ix) * (lnz[0] + 2) + iz;
          int mm = ix * (lnz[0] + 2) + iz;
          yl[nn] = b[mm];
        }
      }
    }
    bout::globals::mpi->MPI_Allreduce(std::begin(yl), std::begin(yg), dimg, MPI_DOUBLE,
                                      MPI_SUM, comm2D);

    int nz = (xProcI % rMG->zNP) * (rMG->lnz[level]);
    BOUT_OMP_PERF(parallel default(shared))
    {
      int xend = rMG->lnx[level] + 1;
      int zend = rMG->lnz[level] + 1;
      BOUT_OMP_PERF(for collapse(2))
      for (int ix = 1; ix < xend; ix++) {
        for (int iz = 1; iz < zend; iz++) {
          int nn = ix * (lnz[0] + 2) + nz + iz;
          int mm = ix * (rMG->lnz[level] + 2) + iz;
          r[mm] = yg[nn];
        }
      }
    }

    rMG->getSolution(std::begin(y), std::begin(r), 1);

    BOUT_OMP_PERF(parallel default(shared))
    {
    BOUT_OMP_PERF(for)
    for (int i = 0; i < dimg; i++) {
      yl[i] = 0.0;
      yg[i] = 0.0;
    }

    int xend = rMG->lnx[level] + 1;
    int zend = rMG->lnz[level] + 1;
      BOUT_OMP_PERF(for collapse(2))
      for (int ix = 1; ix < xend; ix++) {
        for (int iz = 1; iz < zend; iz++) {
          int nn = ix * (lnz[0] + 2) + nz + iz;
          int mm = ix * (rMG->lnz[level] + 2) + iz;
          yl[nn] = y[mm];
        }
      }
    }
    bout::globals::mpi->MPI_Allreduce(std::begin(yl), std::begin(yg), dimg, MPI_DOUBLE,
                                      MPI_SUM, comm2D);

    BOUT_OMP_PERF(parallel default(shared))
    {
      int xend = lnx[0] + 1;
      int zend = lnz[0] + 1;
      BOUT_OMP_PERF(for collapse(2))
      for (int ix = 1; ix < xend; ix++) {
        for (int iz = 1; iz < zend; iz++) {
          int nn = (nx + ix) * (lnz[0] + 2) + iz;
          int mm = ix * (lnz[0] + 2) + iz;
          x[mm] = yg[nn];
        }
      }
    }
    communications(x, 0);
  } else if (kflag == 2) {
    int level = sMG->mglevel - 1;
    int dim = (sMG->lnx[level] + 2) * (sMG->lnz[level] + 2);
    Array<BoutReal> y(dim);
    Array<BoutReal> r(dim);
    int nx = xProcI * lnx[0];
    BOUT_OMP_PERF(parallel default(shared))
    {
      BOUT_OMP_PERF(for)
      for (int i = 0; i < dim; i++) {
        y[i] = 0.0;
        r[i] = 0.0;
      }
      int xend = lnx[0] + 1;
      int zend = lnz[0] + 1;
      BOUT_OMP_PERF(for collapse(2))
      for (int ix = 1; ix < xend; ix++) {
        for (int iz = 1; iz < zend; iz++) {
          int nn = (nx + ix) * (lnz[0] + 2) + iz;
          int mm = ix * (lnz[0] + 2) + iz;
          y[nn] = b[mm];
        }
      }
    }
    bout::globals::mpi->MPI_Allreduce(std::begin(y), std::begin(r), dim, MPI_DOUBLE,
                                      MPI_SUM, commMG);
    BOUT_OMP_PERF(parallel default(shared))
    BOUT_OMP_PERF(for)
    for (int i = 0; i < dim; i++) {
      y[i] = 0.0;
    }
    sMG->getSolution(std::begin(y), std::begin(r), 1);

    BOUT_OMP_PERF(parallel default(shared))
    {
      int xend = lnx[0] + 1;
      int zend = lnz[0] + 1;
      BOUT_OMP_PERF(for collapse(2))
      for (int ix = 1; ix < xend; ix++) {
        for (int iz = 1; iz < zend; iz++) {
          int nn = (nx + ix) * (lnz[0] + 2) + iz;
          int mm = ix * (lnz[0] + 2) + iz;
          x[mm] = y[nn];
        }
      }
    }
    communications(x, 0);
  } else {
    pGMRES(x, b, 0, 0);
  }
}

void Multigrid1DP::convertMatrixF2D(int level) {

  int ggx = rMG->lnx[level];
  int dim = (ggx + 2) * (gnz[0] + 2);
  Array<BoutReal> yl(dim * 9);
  Array<BoutReal> yg(dim * 9);
  int nx = (xProcI % rMG->zNP) * lnx[0];
  BOUT_OMP_PERF(parallel default(shared))
  {
    BOUT_OMP_PERF(for)
    for (int i = 0; i < dim * 9; i++) {
      yl[i] = 0.0;
      yg[i] = 0.0;
    }
    BOUT_OMP_PERF(for)
    for (int i = 0; i < (rMG->lnx[level] + 2) * (rMG->lnz[level] + 2) * 9; i++) {
      rMG->matmg[level][i] = 0.0;
    }

    int xend = lnx[0] + 1;
    int zend = lnz[0] + 1;
    BOUT_OMP_PERF(for collapse(2))
    for (int ix = 1; ix < xend; ix++) {
      for (int iz = 1; iz < zend; iz++) {
        int nn = (nx + ix) * (lnz[0] + 2) + iz;
        int mm = ix * (lnz[0] + 2) + iz;
        for (int k = 0; k < 9; k++) {
          yl[nn * 9 + k] = matmg[0][mm * 9 + k];
        }
      }
    }
  }
  if (pcheck == 3) {
    FILE* outf;
    char outfile[256];
    sprintf(outfile, "2DP_CP_%d.mat", rProcI);
    output << "Out file= " << outfile << endl;
    outf = fopen(outfile, "w");
    fprintf(outf, "dim = (%d, %d)\n", ggx, gnz[0]);

    for (int ii = 0; ii < dim; ii++) {
      fprintf(outf, "%d ==", ii);
      for (int j = 0; j < 9; j++) {
        fprintf(outf, "%12.6f,", yl[ii * 9 + j]);
      }
      fprintf(outf, "\n");
    }
    fclose(outf);
  }
  bout::globals::mpi->MPI_Allreduce(std::begin(yl), std::begin(yg), dim * 9, MPI_DOUBLE,
                                    MPI_SUM, comm2D);

  if (pcheck == 3) {
    FILE* outf;
    char outfile[256];
    sprintf(outfile, "2DP_Conv_%d.mat", rProcI);
    output << "Out file= " << outfile << endl;
    outf = fopen(outfile, "w");
    fprintf(outf, "dim = (%d, %d)\n", ggx, gnz[0]);

    for (int ii = 0; ii < dim; ii++) {
      fprintf(outf, "%d ==", ii);
      for (int j = 0; j < 9; j++) {
        fprintf(outf, "%12.6f,", yg[ii * 9 + j]);
      }
      fprintf(outf, "\n");
    }
    fclose(outf);
  }
  int nz = (xProcI % rMG->zNP) * (rMG->lnz[level]);

  BOUT_OMP_PERF(parallel default(shared))
  {
    int xend = rMG->lnx[level] + 1;
    int zend = rMG->lnz[level] + 1;
    BOUT_OMP_PERF(for collapse(2))
    for (int ix = 1; ix < xend; ix++) {
      for (int iz = 1; iz < zend; iz++) {
        int nn = ix * (lnz[0] + 2) + nz + iz;
        int mm = ix * (rMG->lnz[level] + 2) + iz;
        for (int k = 0; k < 9; k++) {
          rMG->matmg[level][mm * 9 + k] = yg[nn * 9 + k];
        }
      }
    }
  }
}

void Multigrid1DP::convertMatrixFS(int level) {

  int dim = (gnx[0] + 2) * (gnz[0] + 2);
  Array<BoutReal> yl(dim * 9);
  BoutReal* yg = sMG->matmg[level];
  int nx = xProcI * lnx[0];
  BOUT_OMP_PERF(parallel default(shared))
  {
    BOUT_OMP_PERF(for)
    for (int i = 0; i < dim * 9; i++) {
      yl[i] = 0.0;
      yg[i] = 0.0;
    }
    int xend = lnx[0] + 1;
    int zend = lnz[0] + 1;
    BOUT_OMP_PERF(for collapse(2))
    for (int ix = 1; ix < xend; ix++) {
      for (int iz = 1; iz < zend; iz++) {
        int nn = (nx + ix) * (lnz[0] + 2) + iz;
        int mm = ix * (lnz[0] + 2) + iz;
        for (int k = 0; k < 9; k++) {
          yl[nn * 9 + k] = matmg[0][mm * 9 + k];
        }
      }
    }
  }
  bout::globals::mpi->MPI_Allreduce(std::begin(yl), yg, dim * 9, MPI_DOUBLE, MPI_SUM,
                                    commMG);
}

Multigrid2DPf1D::Multigrid2DPf1D(int level, int lx, int lz, int gx, int gz, int dl,
                                 int px, int pz, MPI_Comm comm, int check)
    : MultigridAlg(level, lx, lz, gx, gz, comm, check) {

  mglevel = level;

  /* Momory allocate for Multigrid */
  xNP = px;
  zNP = pz;
  numP = px * pz;
  commMG = comm;
  MPI_Comm_rank(comm, &rProcI);
  xProcI = rProcI / zNP;
  zProcI = rProcI % zNP;
  if (xProcI == 0) {
    xProcM = numP - zNP + zProcI;
  } else {
    xProcM = rProcI - zNP;
  }
  if (xProcI == xNP - 1) {
    xProcP = zProcI;
  } else {
    xProcP = rProcI + zNP;
  }
  if (zProcI == 0) {
    zProcM = rProcI + zNP - 1;
  } else {
    zProcM = rProcI - 1;
  }
  if (zProcI == zNP - 1) {
    zProcP = xProcI * zNP;
  } else {
    zProcP = rProcI + 1;
  }
  if (pcheck == 2) {
    output << "In 2DP " << mglevel << "xNP=" << xNP << "(" << zNP << ")" << dl << endl;
    for (int i = mglevel - 1; i >= 0; i--) {
      output << i << " loc dim " << lnx[i] << "," << lnz[i] << endl;
      output << i << " glo dim " << gnx[i] << "," << gnz[i] << endl;
    }
  }

  if (dl > 0) {
    int nn = gnx[0];
    int mm = gnz[0];
    int kk = 1;
    for (int n = dl; n > 0; n--) {
      if ((nn % 2 == 0) && (mm % 2 == 0)) {
        kk += 1;
        nn = nn / 2;
        mm = mm / 2;
      } else {
        n = 0;
      }
    }
    if (pcheck == 2) {
      output << "In 2DP To Ser" << kk << "xNP=" << xNP << "(" << zNP << ")" << endl;
      output << "total dim" << gnx[0] << "(" << gnz[0] << ")" << endl;
    }
    kflag = 2;
    sMG = bout::utils::make_unique<MultigridSerial>(kk, gnx[0], gnz[0], commMG, pcheck);
  } else {
    kflag = 0;
  }
}

void Multigrid2DPf1D::setMultigridC(int UNUSED(plag)) {

  int level = mglevel - 1;
  for (int n = level; n > 0; n--) {
    setMatrixC(n);
    if (pcheck == 2) {
      output << n << "Cmatrix in 2DP = " << lnx[n - 1] << "," << lnz[n - 1] << endl;
      output << gnx[n - 1] << "," << gnz[n - 1] << endl;
    }
  }
  if (kflag == 2) {
    level = sMG->mglevel - 1;
    convertMatrixFS(level);
    if (level > 0) {
      sMG->setMultigridC(0);
    }
    if (pcheck == 2) {
      for (int i = level; i >= 0; i--) {
        FILE* outf;
        char outfile[256];
        sprintf(outfile, "S2D_matC%1d_%d.mat", i, sMG->rProcI);
        outf = fopen(outfile, "w");
        output << "Out file= " << outfile << endl;
        int dim = (sMG->lnx[i] + 2) * (sMG->lnz[i] + 2);
        fprintf(outf, "dim = %d\n", dim);

        for (int ii = 0; ii < dim; ii++) {
          fprintf(outf, "%d ==", ii);
          for (int j = 0; j < 9; j++) {
            fprintf(outf, "%12.6f,", sMG->matmg[i][ii * 9 + j]);
          }
          fprintf(outf, "\n");
        }
        fclose(outf);
      }
    }
  }
}

void Multigrid2DPf1D::setValueS() {
  if (kflag == 2) {
    sMG->mgplag = mgplag;
    sMG->mgsm = mgsm;
    sMG->cftype = cftype;
    sMG->rtol = rtol;
    sMG->atol = atol;
    sMG->dtol = dtol;
    sMG->omega = omega;
  }
}

void Multigrid2DPf1D::setPcheck(int check) {

  pcheck = check;
  if (kflag == 2) {
    sMG->pcheck = check;
  }
}

void Multigrid2DPf1D::lowestSolver(BoutReal* x, BoutReal* b, int UNUSED(plag)) {

  if (kflag == 2) {
    int level = sMG->mglevel - 1;
    int dim = (sMG->lnx[level] + 2) * (sMG->lnz[level] + 2);
    Array<BoutReal> y(dim);
    Array<BoutReal> r(dim);
    int nx = xProcI * lnx[0];
    int nz = zProcI * lnz[0];
    BOUT_OMP_PERF(parallel default(shared))
    {
      BOUT_OMP_PERF(for)
      for (int i = 0; i < dim; i++) {
        y[i] = 0.0;
        r[i] = 0.0;
      }

      int xend = lnx[0] + 1;
      int zend = lnz[0] + 1;
      BOUT_OMP_PERF(for collapse(2))
      for (int ix = 1; ix < xend; ix++) {
        for (int iz = 1; iz < zend; iz++) {
          int nn = (nx + ix) * (gnz[0] + 2) + nz + iz;
          int mm = ix * (lnz[0] + 2) + iz;
          y[nn] = b[mm];
        }
      }
    }
    bout::globals::mpi->MPI_Allreduce(std::begin(y), std::begin(r), dim, MPI_DOUBLE,
                                      MPI_SUM, commMG);
    BOUT_OMP_PERF(parallel default(shared))
    BOUT_OMP_PERF(for)
    for (int i = 0; i < dim; i++) {
      y[i] = 0.0;
    }
    sMG->getSolution(std::begin(y), std::begin(r), 1);
    BOUT_OMP_PERF(parallel default(shared))
    {
      int xend = lnx[0] + 1;
      int zend = lnz[0] + 1;
      BOUT_OMP_PERF(for collapse(2))
      for (int ix = 1; ix < xend; ix++) {
        for (int iz = 1; iz < zend; iz++) {
          int nn = (nx + ix) * (gnz[0] + 2) + nz + iz;
          int mm = ix * (lnz[0] + 2) + iz;
          x[mm] = y[nn];
        }
      }
    }
    communications(x, 0);
  } else {
    pGMRES(x, b, 0, 0);
  }
}

void Multigrid2DPf1D::convertMatrixFS(int level) {

  int dim = (gnx[0] + 2) * (gnz[0] + 2);
  Array<BoutReal> yl(dim * 9);
  BoutReal* yg = sMG->matmg[level];
  int nx = xProcI * lnx[0];
  int nz = zProcI * lnz[0];
  BOUT_OMP_PERF(parallel default(shared))
  {
    BOUT_OMP_PERF(for)
    for (int i = 0; i < dim * 9; i++) {
      yl[i] = 0.0;
      yg[i] = 0.0;
    }
    int xend = lnx[0] + 1;
    int zend = lnz[0] + 1;
    BOUT_OMP_PERF(for collapse(2))
    for (int ix = 1; ix < xend; ix++) {
      for (int iz = 1; iz < zend; iz++) {
        int nn = (nx + ix) * (gnz[0] + 2) + nz + iz;
        int mm = ix * (lnz[0] + 2) + iz;
        for (int k = 0; k < 9; k++) {
          yl[nn * 9 + k] = matmg[0][mm * 9 + k];
        }
      }
    }
  }
  bout::globals::mpi->MPI_Allreduce(std::begin(yl), yg, dim * 9, MPI_DOUBLE, MPI_SUM,
                                    commMG);
}

MultigridSerial::MultigridSerial(int level, int gx, int gz, MPI_Comm comm, int check)
    : MultigridAlg(level, gx, gz, gx, gz, comm, check) {

  xNP = 1;
  zNP = 1;
  numP = 1;
  MPI_Comm_rank(commMG, &rProcI);
  xProcI = rProcI;
  zProcI = rProcI;
  xProcM = rProcI;
  xProcP = rProcI;
  zProcM = rProcI;
  zProcP = rProcI;
  if (pcheck == 2) {
    output << "In SerMG " << mglevel << "xNP=" << xNP << "(" << zNP << ")" << endl;
    for (int i = mglevel - 1; i >= 0; i--) {
      output << i << " Ser loc dim " << lnx[i] << "," << lnz[i] << endl;
      output << i << " Ser glo dim " << gnx[i] << "," << gnz[i] << endl;
    }
  }
}

#endif
