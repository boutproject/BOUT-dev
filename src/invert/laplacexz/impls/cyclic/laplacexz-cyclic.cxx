#include "laplacexz-cyclic.hxx"

#include <utils.hxx>
#include <fft.hxx>
#include <bout/constants.hxx>
#include <bout/sys/timer.hxx>
#include <msg_stack.hxx>

#include <output.hxx>

LaplaceXZcyclic::LaplaceXZcyclic(Mesh *m, Options *options) : LaplaceXZ(m, options), mesh(m) {

  // Number of Z Fourier modes, including DC
  nmode = (m->LocalNz) / 2 + 1;

  // Number of independent systems of
  // equations to solve
  nsys = nmode * (m->yend - m->ystart + 1);

  // Start and end X index
  xstart = m->xstart; // Starting X index
  if (m->firstX()) {
    xstart -= 1;
  }
  xend = m->xend;
  if (m->lastX()) {
    xend += 1;
  }

  // Number of points in X on this processor
  // including boundaries but not guard cells
  nloc = xend - xstart + 1;

  acoef  = Matrix<dcomplex>(nsys, nloc);
  bcoef  = Matrix<dcomplex>(nsys, nloc);
  ccoef  = Matrix<dcomplex>(nsys, nloc);
  xcmplx = Matrix<dcomplex>(nsys, nloc);
  rhscmplx = Matrix<dcomplex>(nsys, nloc);

  k1d = Array<dcomplex>((m->LocalNz) / 2 + 1);
  k1d_2 = Array<dcomplex>((m->LocalNz) / 2 + 1);

  // Create a cyclic reduction object, operating on dcomplex values
  cr = new CyclicReduce<dcomplex>(mesh->getXcomm(), nloc);

  // Getting the boundary flags
  OPTION(options, inner_boundary_flags, 0);
  OPTION(options, outer_boundary_flags, 0);

  // Set default coefficients
  setCoefs(Field2D(1.0), Field2D(0.0));
}

LaplaceXZcyclic::~LaplaceXZcyclic() {
  // Delete tridiagonal solver
  delete cr;
}

void LaplaceXZcyclic::setCoefs(const Field2D &A2D, const Field2D &B2D) {
  TRACE("LaplaceXZcyclic::setCoefs");
  Timer timer("invert");
  
  // Set coefficients

  Coordinates *coord = mesh->coordinates();

  // NOTE: For now the X-Z terms are omitted, so check that they are small
  ASSERT2(max(abs(coord->g13)) < 1e-5);
  
  int ind = 0;
  for(int y=mesh->ystart; y <= mesh->yend; y++) {
    for(int kz = 0; kz < nmode; kz++) {
      BoutReal kwave=kz*2.0*PI/(coord->zlength());

      if(mesh->firstX()) {
        // Inner X boundary
        
        if( ((kz == 0) && (inner_boundary_flags & INVERT_DC_GRAD)) ||
            ((kz != 0) && (inner_boundary_flags & INVERT_AC_GRAD)) ) {
          // Neumann
          acoef(ind, 0) = 0.0;
          bcoef(ind, 0) = 1.0;
          ccoef(ind, 0) = -1.0;
        }else {
          // Dirichlet
          // This includes cases where the boundary is set to a value
          acoef(ind, 0) = 0.0;
          bcoef(ind, 0) = 0.5;
          ccoef(ind, 0) = 0.5;
        }
      }

      // Bulk of the domain
      for(int x=mesh->xstart; x <= mesh->xend; x++) {
        acoef(ind, x - xstart) = 0.0; // X-1
        bcoef(ind, x - xstart) = 0.0; // Diagonal
        ccoef(ind, x - xstart) = 0.0; // X+1

        //////////////////////////////
        // B coefficient
        bcoef(ind, x - xstart) += B2D(x, y);

        //////////////////////////////
        // A coefficient

        // XX component

        // Metrics on x+1/2 boundary
        BoutReal J = 0.5*(coord->J(x,y) + coord->J(x+1,y));
        BoutReal g11 = 0.5*(coord->g11(x,y) + coord->g11(x+1,y));
        BoutReal dx = 0.5*(coord->dx(x,y) + coord->dx(x+1,y));
        BoutReal A = 0.5*(A2D(x,y) + A2D(x+1,y));

        BoutReal val = A * J * g11 / (coord->J(x,y) * dx * coord->dx(x,y));

        ccoef(ind, x - xstart) += val;
        bcoef(ind, x - xstart) -= val;

        // Metrics on x-1/2 boundary
        J = 0.5*(coord->J(x,y) + coord->J(x-1,y));
        g11 = 0.5*(coord->g11(x,y) + coord->g11(x-1,y));
        dx = 0.5*(coord->dx(x,y) + coord->dx(x-1,y));
        A = 0.5*(A2D(x,y) + A2D(x-1,y));

        val = A * J * g11 / (coord->J(x,y) * dx * coord->dx(x,y));
        acoef(ind, x - xstart) += val;
        bcoef(ind, x - xstart) -= val;

        // ZZ component
        bcoef(ind, x - xstart) -= A2D(x, y) * SQ(kwave) * coord->g33(x, y);
      }

      // Outer X boundary
      if(mesh->lastX()) {
        // Outer X boundary
        if( ((kz == 0) && (outer_boundary_flags & INVERT_DC_GRAD)) ||
            ((kz != 0) && (outer_boundary_flags & INVERT_AC_GRAD)) ) {
          // Neumann
          acoef(ind, nloc - 1) = -1.0;
          bcoef(ind, nloc - 1) = 1.0;
          ccoef(ind, nloc - 1) = 0.0;
        }else {
          // Dirichlet
          acoef(ind, nloc - 1) = 0.5;
          bcoef(ind, nloc - 1) = 0.5;
          ccoef(ind, nloc - 1) = 0.0;
        }
      }
      
      ind++;
    }
  }
  // Set coefficients in tridiagonal solver
  cr->setCoefs(acoef, bcoef, ccoef);
}

Field3D LaplaceXZcyclic::solve(const Field3D &rhs, const Field3D &x0) {
  Timer timer("invert");

  Mesh *mesh = rhs.getMesh();
  // Create the rhs array
  int ind = 0;
  for(int y=mesh->ystart; y <= mesh->yend; y++) {

    if(mesh->firstX()) {
      // Inner X boundary
      
      if(inner_boundary_flags & INVERT_SET) {
        // Fourier transform x0 in Z at xstart-1 and xstart
        rfft(&x0(mesh->xstart - 1, y, 0), mesh->LocalNz, std::begin(k1d));
        rfft(&x0(mesh->xstart, y, 0), mesh->LocalNz, std::begin(k1d_2));
        for(int kz = 0; kz < nmode; kz++) {
          // Use the same coefficients as applied to the solution
          // so can either set gradient or value
          rhscmplx(ind + kz, 0) =
              bcoef(ind + kz, 0) * k1d[kz] + ccoef(ind + kz, 0) * k1d_2[kz];
        }
      }else if(inner_boundary_flags & INVERT_RHS) {
        // Fourier transform rhs in Z at xstart-1 and xstart
        rfft(&rhs(mesh->xstart - 1, y, 0), mesh->LocalNz, std::begin(k1d));
        rfft(&rhs(mesh->xstart, y, 0), mesh->LocalNz, std::begin(k1d_2));
        for(int kz = 0; kz < nmode; kz++) {
          // Use the same coefficients as applied to the solution
          // so can either set gradient or value
          rhscmplx(ind + kz, 0) =
              bcoef(ind + kz, 0) * k1d[kz] + ccoef(ind + kz, 0) * k1d_2[kz];
        }
      }else {
        for(int kz = 0; kz < nmode; kz++) {
          rhscmplx(ind + kz, 0) = 0.0;
        }
      }
    }

    // Bulk of the domain
    for(int x=mesh->xstart; x <= mesh->xend; x++) {
      // Fourier transform RHS
      rfft(&rhs(x, y, 0), mesh->LocalNz, std::begin(k1d));
      for(int kz = 0; kz < nmode; kz++) {
        rhscmplx(ind + kz, x - xstart) = k1d[kz];
      }
    }

    // Outer X boundary
    if(mesh->lastX()) {
      // Outer X boundary
      if(outer_boundary_flags & INVERT_SET) {
        // Fourier transform x0 in Z at xend and xend+1
        rfft(&x0(mesh->xend, y, 0), mesh->LocalNz, std::begin(k1d));
        rfft(&x0(mesh->xend + 1, y, 0), mesh->LocalNz, std::begin(k1d_2));
        for(int kz = 0; kz < nmode; kz++) {
          // Use the same coefficients as applied to the solution
          // so can either set gradient or value
          rhscmplx(ind + kz, nloc - 1) =
              acoef(ind + kz, nloc - 1) * k1d[kz] + bcoef(ind + kz, nloc - 1) * k1d_2[kz];
        }
      }else if(outer_boundary_flags & INVERT_RHS) {
        // Fourier transform rhs in Z at xstart-1 and xstart
        rfft(&rhs(mesh->xend, y, 0), mesh->LocalNz, std::begin(k1d));
        rfft(&rhs(mesh->xend + 1, y, 0), mesh->LocalNz, std::begin(k1d_2));
        for(int kz = 0; kz < nmode; kz++) {
          // Use the same coefficients as applied to the solution
          // so can either set gradient or value
          rhscmplx(ind + kz, nloc - 1) =
              acoef(ind + kz, nloc - 1) * k1d[kz] + bcoef(ind + kz, nloc - 1) * k1d_2[kz];
        }
      }else {
        for(int kz = 0; kz < nmode; kz++) {
          rhscmplx(ind + kz, nloc - 1) = 0.0;
        }
      }
      
      for(int kz = 0; kz < nmode; kz++) {
        rhscmplx(ind + kz, nloc - 1) = 0.0;
      }
    }
    ind += nmode;
  }

  // Solve tridiagonal systems
  cr->solve(rhscmplx, xcmplx);

  // FFT back to real space

  Field3D result(mesh);
  result.allocate();

  ind = 0;
  for(int y=mesh->ystart; y <= mesh->yend; y++) {
    for(int x=xstart;x<=xend;x++) {
      for(int kz = 0; kz < nmode; kz++) {
        k1d[kz] = xcmplx(ind + kz, x - xstart);
      }

      // This shifts back to field-aligned coordinates
      irfft(std::begin(k1d), mesh->LocalNz, &result(x, y, 0));
    }
    ind += nmode;
  }
  
  return result;
}
