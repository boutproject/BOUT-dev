#include "laplacexz-cyclic.hxx"

#include <utils.hxx>
#include <fft.hxx>
#include <bout/constants.hxx>
#include <bout/sys/timer.hxx>

#include <output.hxx>

LaplaceXZcyclic::LaplaceXZcyclic(Mesh *m, Options *options) : LaplaceXZ(m, options), mesh(m) {

  // Number of Z Fourier modes, including DC
  nmode = (mesh->ngz-1)/2 + 1;

  // Number of independent systems of
  // equations to solve
  nsys = nmode * (mesh->yend - mesh->ystart+1);

  // Start and end X index
  xstart = mesh->xstart; // Starting X index
  if(mesh->firstX()) {
    xstart -= 1;
  }
  xend = mesh->xend;
  if(mesh->lastX()) {
    xend += 1;
  }

  // Number of points in X on this processor
  // including boundaries but not guard cells
  nloc = xend - xstart + 1;

  acoef  = matrix<dcomplex>(nsys, nloc);
  bcoef  = matrix<dcomplex>(nsys, nloc);
  ccoef  = matrix<dcomplex>(nsys, nloc);
  xcmplx = matrix<dcomplex>(nsys, nloc);
  rhscmplx = matrix<dcomplex>(nsys, nloc);

  k1d = new dcomplex[(mesh->ngz-1)/2 + 1];

  // Create a cyclic reduction object, operating on dcomplex values
  cr = new CyclicReduce<dcomplex>(mesh->getXcomm(), nloc);

  // Set default coefficients
  setCoefs(1.0, 0.0);
}

LaplaceXZcyclic::~LaplaceXZcyclic() {
  // Free coefficient arrays
  free_matrix(acoef);
  free_matrix(bcoef);
  free_matrix(ccoef);
  free_matrix(xcmplx);
  free_matrix(rhscmplx);

  delete[] k1d;

  // Delete tridiagonal solver
  delete cr;
}

void LaplaceXZcyclic::setCoefs(const Field2D &A2D, const Field2D &B2D) {
  Timer timer("invert");

  // Set coefficients

  int ind = 0;
  for(int y=mesh->ystart; y <= mesh->yend; y++) {
    for(int kz = 0; kz < nmode; kz++) {
      BoutReal kwave=kz*2.0*PI/(mesh->zlength());

      if(mesh->firstX()) {
        // Inner X boundary

        acoef[ind][0] =  0.0;
        bcoef[ind][0] =  1.0;
        ccoef[ind][0] =  -1.0;
      }

      // Bulk of the domain
      for(int x=mesh->xstart; x <= mesh->xend; x++) {
        acoef[ind][x-xstart] = 0.0; // X-1
        bcoef[ind][x-xstart] = 0.0; // Diagonal
        ccoef[ind][x-xstart] = 0.0; // X+1

        //////////////////////////////
        // B coefficient
        bcoef[ind][x-xstart] += B2D(x,y);

        //////////////////////////////
        // A coefficient

        // XX component

        // Metrics on x+1/2 boundary
        BoutReal J = 0.5*(mesh->J(x,y) + mesh->J(x+1,y));
        BoutReal g11 = 0.5*(mesh->g11(x,y) + mesh->g11(x+1,y));
        BoutReal dx = 0.5*(mesh->dx(x,y) + mesh->dx(x+1,y));
        BoutReal A = 0.5*(A2D(x,y) + A2D(x+1,y));

        BoutReal val = A * J * g11 / (mesh->J(x,y) * dx * mesh->dx(x,y));

        ccoef[ind][x-xstart] += val;
        bcoef[ind][x-xstart] -= val;

        // Metrics on x-1/2 boundary
        J = 0.5*(mesh->J(x,y) + mesh->J(x-1,y));
        g11 = 0.5*(mesh->g11(x,y) + mesh->g11(x-1,y));
        dx = 0.5*(mesh->dx(x,y) + mesh->dx(x-1,y));
        A = 0.5*(A2D(x,y) + A2D(x-1,y));

        val = A * J * g11 / (mesh->J(x,y) * dx * mesh->dx(x,y));
        acoef[ind][x-xstart] += val;
        bcoef[ind][x-xstart] -= val;

        // ZZ component
        bcoef[ind][x-xstart] -= A2D(x,y) * SQ(kwave) * mesh->g33(x,y);

      }

      // Outer X boundary
      if(mesh->lastX()) {
        // Outer X boundary

        acoef[ind][nloc-1] =  1.0;
        bcoef[ind][nloc-1] =  1.0;
        ccoef[ind][nloc-1] =  0.0;
      }

      /*
      if(y == mesh->ystart) {
        for(int i=0;i<nloc;i++) {
          output << i << ": " <<  acoef[ind][i] << ", " << bcoef[ind][i] << ", " << ccoef[ind][i] << endl;
        }
      }
      */

      ind++;
    }
  }
  // Set coefficients in tridiagonal solver
  cr->setCoefs(nsys, acoef, bcoef, ccoef);
}

Field3D LaplaceXZcyclic::solve(const Field3D &rhs, const Field3D &x0) {
  Timer timer("invert");
  
  // Create the rhs array
  int ind = 0;
  for(int y=mesh->ystart; y <= mesh->yend; y++) {

    if(mesh->firstX()) {
      // Inner X boundary

      for(int kz = 0; kz < nmode; kz++) {
        rhscmplx[ind + kz][0] = 0.0;
      }
    }

    // Bulk of the domain
    for(int x=mesh->xstart; x <= mesh->xend; x++) {
      // Fourier transform RHS, shifting into X-Z orthogonal coordinates
      ZFFT(&rhs(x,y,0), mesh->zShift(x, y), k1d);
      for(int kz = 0; kz < nmode; kz++) {
        rhscmplx[ind + kz][x-xstart] = k1d[kz];
      }
    }

    // Outer X boundary
    if(mesh->lastX()) {
      // Outer X boundary

      for(int kz = 0; kz < nmode; kz++) {
        rhscmplx[ind + kz][nloc-1] = 0.0;
      }
    }
    ind += nmode;
  }

  // Solve tridiagonal systems
  cr->solve(nsys, rhscmplx, xcmplx);

  // FFT back to real space

  Field3D result;
  result.allocate();

  ind = 0;
  for(int y=mesh->ystart; y <= mesh->yend; y++) {
    for(int x=xstart;x<=xend;x++) {
      for(int kz = 0; kz < nmode; kz++) {
        k1d[kz] = xcmplx[ind + kz][x-xstart];
      }

      // This shifts back to field-aligned coordinates
      ZFFT_rev(k1d, mesh->zShift(x, y), &result(x,y,0));
    }
    ind += nmode;
  }
  
  return result;
}
