#include <bout/constants.hxx>
#include <bout/mesh.hxx>
#include <boundary_standard.hxx>
#include <boutexception.hxx>
#include <derivs.hxx>
#include <fft.hxx>
#include <globals.hxx>
#include <invert_laplace.hxx>
#include <msg_stack.hxx>
#include <output.hxx>
#include <utils.hxx>

using bout::generator::Context;

// #define BOUNDARY_CONDITIONS_UPGRADE_EXTRAPOLATE_FOR_2ND_ORDER

///////////////////////////////////////////////////////////////
// Helpers

/** \brief Check that there are sufficient non-boundary points for desired B.C.

    Checks both the size of the global grid (i.e. if this B.C. could be ok
    for some parallel setup or not) and the local grid.

    Note the local grid check is not strictly necessary as this would typically
    lead to an out of bounds access error later but we add it here to provide a
    more explanatory message.
 */
#if CHECK > 0
void verifyNumPoints(BoundaryRegion* region, int ptsRequired) {
  TRACE("Verifying number of points available for BC");

  int ptsAvailGlobal, ptsAvailLocal, ptsAvail;
  std::string side, gridType;
  Mesh* mesh = region->localmesh;

  // Initialise var in case of no match and CHECK<=2
  ptsAvail = ptsRequired; // Ensures test passes without exception

  switch (region->location) {
  case BNDRY_XIN:
  case BNDRY_XOUT: {
    side = "x";

    // Here 2*mesh->xstart is the total number of guard/boundary cells
    ptsAvailGlobal = mesh->GlobalNx - 2 * mesh->xstart;

    // Work out how many processor local points we have excluding boundaries
    // but including ghost/guard cells
    ptsAvailLocal = mesh->LocalNx;
    if (mesh->firstX())
      ptsAvailLocal -= mesh->xstart;
    if (mesh->lastX())
      ptsAvailLocal -= mesh->xstart;

    // Now decide if it's a local or global limit, prefer global if a tie
    if (ptsAvailGlobal <= ptsAvailLocal) {
      ptsAvail = ptsAvailGlobal;
      gridType = "global";
    } else {
      ptsAvail = ptsAvailLocal;
      gridType = "local";
    }

    break;
  }
  case BNDRY_YUP:
  case BNDRY_YDOWN: {
    side = "y";

    //Here mesh->numberOfYBoundaries()*mesh->ystart is the total number of guard/boundary
    //cells
    ptsAvailGlobal = mesh->GlobalNy - mesh->numberOfYBoundaries()*2*mesh->ystart;

    // Work out how many processor local points we have excluding boundaries
    // but including ghost/guard cells
    ptsAvailLocal = mesh->LocalNy;
    if (mesh->firstY())
      ptsAvailLocal -= mesh->ystart;
    if (mesh->lastY())
      ptsAvailLocal -= mesh->ystart;

    // Now decide if it's a local or global limit, prefer global if a tie
    if (ptsAvailGlobal <= ptsAvailLocal) {
      ptsAvail = ptsAvailGlobal;
      gridType = "global";
    } else {
      ptsAvail = ptsAvailLocal;
      gridType = "local";
    }

    break;
  }
  default: {
#if CHECK > 2 // Only fail on Unrecognised boundary for extreme checking
    // location is an enum, so cast to int for clarity
    throw BoutException("Unrecognised boundary region ({:d}) for verifyNumPoints.",
                        static_cast<int>(region->location));
#endif
  }
  }

  // Now check we have enough points and if not throw an exception
  if (ptsAvail < ptsRequired) {
    throw BoutException(
        "Too few {:s} grid points for {:s} boundary, have {:d} but need at least {:d}",
        gridType, side, ptsAvail, ptsRequired);
  }
}
#else
// No-op for no checking
void verifyNumPoints(BoundaryRegion*, int) {}
#endif

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryDirichlet::clone(BoundaryRegion* region,
                                     const std::list<std::string>& args) {
  verifyNumPoints(region, 1);

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryDirichlet(region, newgen);
}

void BoundaryDirichlet::apply(Field2D& f) { BoundaryDirichlet::apply(f, 0.); }

void BoundaryDirichlet::apply(Field2D& f, BoutReal t) {
  // Set (at 2nd order) the value at the mid-point between the guard cell and the grid
  // cell to be val N.B. Only first guard cells (closest to the grid) should ever be used

  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  BoutReal val = 0.0;

  // Check for staggered grids

  CELL_LOC loc = f.getLocation();
  if (mesh->StaggerGrids and (loc == CELL_XLOW or loc == CELL_YLOW)) {
    // Staggered. Need to apply slightly differently

    if (loc == CELL_XLOW) {
      // shifted in X

      if (bndry->bx > 0) {
        // Outer x boundary

        for (; !bndry->isDone(); bndry->next1d()) {
          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh));
          }

          f(bndry->x, bndry->y) = val;

          // Need to set second guard cell, as may be used for interpolation or upwinding
          // derivatives
          for (int i = 1; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y;

            f(xi, yi) = 2 * f(xi - bndry->bx, yi) - f(xi - 2 * bndry->bx, yi);
          }
        }
      }
      if (bndry->bx < 0) {
        // Inner x boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {
          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh));
          }

          f(bndry->x - bndry->bx, bndry->y) = val;

          // Need to set second guard cell, as may be used for interpolation or upwinding
          // derivatives
          for (int i = 0; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y;

            f(xi, yi) = 2 * f(xi - bndry->bx, yi) - f(xi - 2 * bndry->bx, yi);
          }
        }
      }
      if (bndry->by != 0) {
        // y boundaries
        for (; !bndry->isDone(); bndry->next1d()) {
          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh));
          }
          f(bndry->x, bndry->y) = 2 * val - f(bndry->x - bndry->bx, bndry->y - bndry->by);

          // Need to set second guard cell, as may be used for interpolation or upwinding
          // derivatives
          for (int i = 1; i < bndry->width; i++) {
            int xi = bndry->x;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 2 * f(xi, yi - bndry->by) - f(xi, yi - 2 * bndry->by);
          }
        }
      }
    } else if (loc == CELL_YLOW) {
      // Y boundary, and field is shifted in Y

      if (bndry->by > 0) {
        // Upper y boundary

        for (; !bndry->isDone(); bndry->next1d()) {
          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh));
          }

          f(bndry->x, bndry->y) = val;

          // Need to set second guard cell, as may be used for interpolation or upwinding
          // derivatives
          for (int i = 1; i < bndry->width; i++) {
            int xi = bndry->x;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 2 * f(xi, yi - bndry->by) - f(xi, yi - 2 * bndry->by);
          }
        }
      }
      if (bndry->by < 0) {
        // Lower y boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {
          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh));
          }

          f(bndry->x, bndry->y - bndry->by) = val;

          // Need to set second guard cell, as may be used for interpolation or upwinding
          // derivatives
          for (int i = 0; i < bndry->width; i++) {
            int xi = bndry->x;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 2 * f(xi, yi - bndry->by) - f(xi, yi - 2 * bndry->by);
          }
        }
      }
      if (bndry->bx != 0) {
        // x boundaries
        for (; !bndry->isDone(); bndry->next1d()) {

          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh));
          }
          f(bndry->x, bndry->y) = 2 * val - f(bndry->x - bndry->bx, bndry->y - bndry->by);

          // Need to set second guard cell, as may be used for interpolation or upwinding
          // derivatives
          for (int i = 1; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y;
            f(xi, yi) = 2 * f(xi - bndry->bx, yi) - f(xi - 2 * bndry->bx, yi);
          }
        }
      }
    }
  } else {
    // Non-staggered, standard case
    for(; !bndry->isDone(); bndry->next1d()) {
      
      if(fg) {
	val = fg->generate(Context(bndry, loc, t, mesh));
      }
      
      f(bndry->x,bndry->y) = 2*val - f(bndry->x-bndry->bx, bndry->y-bndry->by);
			
      // Need to set second guard cell, as may be used for interpolation or upwinding derivatives
      for(int i=1;i<bndry->width;i++) {
	int xi = bndry->x + i*bndry->bx;
	int yi = bndry->y + i*bndry->by;						
	f(xi, yi) = 2*f(xi - bndry->bx, yi - bndry->by) - f(xi - 2*bndry->bx, yi - 2*bndry->by);	
      }
    }
  }
}

void BoundaryDirichlet::apply(Field3D& f) { BoundaryDirichlet::apply(f, 0.); }

void BoundaryDirichlet::apply(Field3D& f, BoutReal t) {
  // Set (at 2nd order) the value at the mid-point between the guard cell and the grid
  // cell to be val N.B. Only first guard cells (closest to the grid) should ever be used

  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  BoutReal val = 0.0;

  // Check for staggered grids

  CELL_LOC loc = f.getLocation();
  if (mesh->StaggerGrids && loc != CELL_CENTRE) {
    // Staggered. Need to apply slightly differently

    if (loc == CELL_XLOW) {
      // X boundary, and field is shifted in X

      if (bndry->bx > 0) {
        // Outer x boundary

        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh));
            }
            f(bndry->x, bndry->y, zk) = val;

            // Need to set second guard cell, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 1; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y;

              f(xi, yi, zk) =
                  2 * f(xi - bndry->bx, yi, zk) - f(xi - 2 * bndry->bx, yi, zk);
            }
          }
        }
      }
      if (bndry->bx < 0) {
        // Inner x boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh));
            }
            f(bndry->x - bndry->bx, bndry->y, zk) = val;
            f(bndry->x, bndry->y, zk) = f(bndry->x - bndry->bx, bndry->y, zk);

            // Need to set second guard cell, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 0; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y;

              f(xi, yi, zk) =
                  2 * f(xi - bndry->bx, yi, zk) - f(xi - 2 * bndry->bx, yi, zk);
            }
          }
        }
      }
      if (bndry->by != 0) {
        // y boundaries
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh));
            }
            f(bndry->x, bndry->y, zk) =
                2 * val - f(bndry->x - bndry->bx, bndry->y - bndry->by, zk);

            // Need to set second guard cell, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 1; i < bndry->width; i++) {
              int xi = bndry->x;
              int yi = bndry->y + i * bndry->by;

              f(xi, yi, zk) =
                  2 * f(xi, yi - bndry->by, zk) - f(xi, yi - 2 * bndry->by, zk);
            }
          }
        }
      }
    } else if (loc == CELL_YLOW) {
      // Shifted in Y

      if (bndry->by > 0) {
        // Upper y boundary boundary

        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh));
            }
            f(bndry->x, bndry->y, zk) = val;

            // Need to set second guard cell, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 1; i < bndry->width; i++) {
              int xi = bndry->x;
              int yi = bndry->y + i * bndry->by;

              f(xi, yi, zk) =
                  2.0 * f(xi, yi - bndry->by, zk) - f(xi, yi - 2 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->by < 0) {
        // Lower y boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh));
            }
            f(bndry->x, bndry->y - bndry->by, zk) = val;

            // Need to set second guard cell, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 0; i < bndry->width; i++) {
              int xi = bndry->x;
              int yi = bndry->y + i * bndry->by;

              f(xi, yi, zk) =
                  2 * f(xi, yi - bndry->by, zk) - f(xi, yi - 2 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->bx != 0) {
        // x boundaries
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh));
            }

            f(bndry->x, bndry->y, zk) =
                2 * val - f(bndry->x - bndry->bx, bndry->y - bndry->by, zk);

            // Need to set second guard cell, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 1; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y;

              f(xi, yi, zk) =
                  2 * f(xi - bndry->bx, yi, zk) - f(xi - 2 * bndry->bx, yi, zk);
            }
          }
        }
      }
    } else if (loc == CELL_ZLOW) {
      // Shifted in Z

      for(; !bndry->isDone(); bndry->next1d()) {
        // Calculate the X and Y normalised values half-way between the guard cell and grid cell
        BoutReal xnorm = 0.5*(   mesh->GlobalX(bndry->x)  // In the guard cell
                                 + mesh->GlobalX(bndry->x - bndry->bx) ); // the grid cell

        BoutReal ynorm = 0.5*(   mesh->GlobalY(bndry->y)  // In the guard cell
                                 + mesh->GlobalY(bndry->y - bndry->by) ); // the grid cell

        for(int zk=0;zk<mesh->LocalNz;zk++) {
          if(fg){
            val = fg->generate(xnorm,TWOPI*ynorm,TWOPI*(zk - 0.5)/(mesh->LocalNz), t);
          }
          f(bndry->x,bndry->y,zk) = 2*val - f(bndry->x-bndry->bx, bndry->y-bndry->by, zk);

          // We've set the first boundary point using extrapolation in
          // the line above.  The below block of code is attempting to
          // set the rest of the boundary cells also using
          // extrapolation. Whilst this choice doesn't impact 2nd order
          // methods it has been observed that with higher order
          // methods, which actually use these points, the use of
          // extrapolation can be unstable. For this reason we have
          // commented out the below block and replaced it with the loop
          // several lines below, which just sets all the rest of the
          // boundary points to be the specified value.  We've not
          // removed the commented out code as we may wish to revisit
          // this in the future, however it may be that this is
          // eventually removed.  It can be noted that we *don't* apply
          // this treatment for other boundary treatments,
          // i.e. elsewhere we tend to extrapolate.

          // // Need to set second guard cell, as may be used for interpolation or upwinding derivatives
          // for(int i=1;i<bndry->width;i++) {
          //   int xi = bndry->x + i*bndry->bx;
          //   int yi = bndry->y + i*bndry->by;

          //   f(xi, yi, zk) = 2*f(xi - bndry->bx, yi - bndry->by, zk) - f(xi - 2*bndry->bx, yi - 2*bndry->by, zk);
          //   // f(xi, yi, zk) = 3.0*f(xi - bndry->bx, yi - bndry->by, zk) - 3.0*f(xi - 2*bndry->bx, yi - 2*bndry->by, zk) + f(xi - 3*bndry->bx, yi - 3*bndry->by, zk);

          // }
        }

        // This loop is our alternative approach to setting the rest of the boundary
        // points. Instead of extrapolating we just use the generated values. This
        // can help with the stability of higher order methods.
        for (int i = 1; i < bndry->width; i++) {
          // Set any other guard cells using the values on the cells
          int xi = bndry->x + i*bndry->bx;
          int yi = bndry->y + i*bndry->by;
          xnorm = mesh->GlobalX(xi);
          ynorm = mesh->GlobalY(yi);
          for(int zk=0;zk<mesh->LocalNz;zk++) {
            if(fg) {
              val = fg->generate(xnorm,TWOPI*ynorm,TWOPI*(zk - 0.5)/(mesh->LocalNz), t);
            }
            f(xi, yi, zk) = val;
          }
        }
      }
    } else {
      throw BoutException("Unrecognised location");
    }
  } else {
    // Standard (non-staggered) case
    for (; !bndry->isDone(); bndry->next1d()) {
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        if (fg) {
          val = fg->generate(Context(bndry, zk, loc, t, mesh));
        }
        f(bndry->x, bndry->y, zk) =
            2 * val - f(bndry->x - bndry->bx, bndry->y - bndry->by, zk);

        // We've set the first boundary point using extrapolation in
        // the line above.  The below block of code is attempting to
        // set the rest of the boundary cells also using
        // extrapolation. Whilst this choice doesn't impact 2nd order
        // methods it has been observed that with higher order
        // methods, which actually use these points, the use of
        // extrapolation can be unstable. For this reason we have
        // commented out the below block and replaced it with the loop
        // several lines below, which just sets all the rest of the
        // boundary points to be the specified value.  We've not
        // removed the commented out code as we may wish to revisit
        // this in the future, however it may be that this is
        // eventually removed.  It can be noted that we *don't* apply
        // this treatment for other boundary treatments,
        // i.e. elsewhere we tend to extrapolate.

        // // Need to set second guard cell, as may be used for interpolation or upwinding
        // derivatives for(int i=1;i<bndry->width;i++) {
        //   int xi = bndry->x + i*bndry->bx;
        //   int yi = bndry->y + i*bndry->by;

        //   f(xi, yi, zk) = 2*f(xi - bndry->bx, yi - bndry->by, zk) - f(xi - 2*bndry->bx,
        //   yi - 2*bndry->by, zk);
        //   // f(xi, yi, zk) = 3.0*f(xi - bndry->bx, yi - bndry->by, zk) - 3.0*f(xi -
        //   2*bndry->bx, yi - 2*bndry->by, zk) + f(xi - 3*bndry->bx, yi - 3*bndry->by,
        //   zk);

        // }
      }

      // This loop is our alternative approach to setting the rest of the boundary
      // points. Instead of extrapolating we just use the generated values. This
      // can help with the stability of higher order methods.
      for (int i = 1; i < bndry->width; i++) {
        // Set any other guard cells using the values on the cells
        int xi = bndry->x + i * bndry->bx;
        int yi = bndry->y + i * bndry->by;
        for (int zk = 0; zk < mesh->LocalNz; zk++) {
          if (fg) {
            val = fg->generate(Context(bndry, zk, loc, t, mesh));
          }
          f(xi, yi, zk) = val;
        }
      }
    }
  }
}

void BoundaryDirichlet::apply_ddt(Field2D& f) {
  Field2D* dt = f.timeDeriv();
  for (bndry->first(); !bndry->isDone(); bndry->next())
    (*dt)(bndry->x, bndry->y) = 0.; // Set time derivative to zero
}

void BoundaryDirichlet::apply_ddt(Field3D& f) {
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  Field3D* dt = f.timeDeriv();

  for (bndry->first(); !bndry->isDone(); bndry->next())
    for (int z = 0; z < mesh->LocalNz; z++)
      (*dt)(bndry->x, bndry->y, z) = 0.; // Set time derivative to zero
}

///////////////////////////////////////////////////////////////
// New implementation, accurate to higher order

BoundaryOp* BoundaryDirichlet_O3::clone(BoundaryRegion* region,
                                        const std::list<std::string>& args) {
  verifyNumPoints(region, 2);
  std::shared_ptr<FieldGenerator> newgen = nullptr;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryDirichlet_O3(region, newgen);
}

void BoundaryDirichlet_O3::apply(Field2D& f) { BoundaryDirichlet_O3::apply(f, 0.); }

void BoundaryDirichlet_O3::apply(Field2D& f, BoutReal t) {
  // Set (at 2nd order) the value at the mid-point between the guard cell and the grid
  // cell to be val N.B. Only first guard cells (closest to the grid) should ever be used

  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  BoutReal val = 0.0;

  // Check for staggered grids

  CELL_LOC loc = f.getLocation();
  if (mesh->StaggerGrids and (loc == CELL_XLOW or loc == CELL_YLOW)) {
    // Staggered. Need to apply slightly differently

    if (loc == CELL_XLOW) {
      // Field is shifted in X

      if (bndry->bx > 0) {
        // Outer x boundary
        for (; !bndry->isDone(); bndry->next1d()) {
          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh));
          }

          f(bndry->x, bndry->y) = val;

          // Need to set second guard cell, as may be used for interpolation or upwinding
          // derivatives
          for (int i = 1; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 3.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + f(xi - 3 * bndry->bx, yi - 3 * bndry->by);
          }
        }
      }
      if (bndry->bx < 0) {
        // Inner x boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {
          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh));
          }
          f(bndry->x - bndry->bx, bndry->y) = val;

          // Need to set second guard cell, as may be used for interpolation or upwinding
          // derivatives
          for (int i = 0; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 3.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + f(xi - 3 * bndry->bx, yi - 3 * bndry->by);
          }
        }
      }
      if (bndry->by != 0) {
        for (; !bndry->isDone(); bndry->next1d()) {
          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh));
          }

          f(bndry->x, bndry->y) =
              (8. / 3) * val - 2. * f(bndry->x - bndry->bx, bndry->y - bndry->by)
              + f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by) / 3.;

          // Need to set second guard cell, as may be used for interpolation or upwinding
          // derivatives
          for (int i = 1; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 3.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + f(xi - 3 * bndry->bx, yi - 3 * bndry->by);
          }
        }
      }
    } else if (loc == CELL_YLOW) {
      // Field is shifted in Y

      if (bndry->by > 0) {
        // Upper y boundary

        for (; !bndry->isDone(); bndry->next1d()) {
          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh));
          }

          f(bndry->x, bndry->y) = val;

          // Need to set second guard cell, as may be used for interpolation or upwinding
          // derivatives
          for (int i = 1; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 3.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + f(xi - 3 * bndry->bx, yi - 3 * bndry->by);
          }
        }
      }
      if (bndry->by < 0) {
        // Lower y boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {
          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh));
          }

          f(bndry->x, bndry->y - bndry->by) = val;
          // Need to set second guard cell, as may be used for interpolation or upwinding
          // derivatives
          for (int i = 0; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 3.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + f(xi - 3 * bndry->bx, yi - 3 * bndry->by);
          }
        }
      }
      if (bndry->bx != 0) {
        // x boundaries
        for (; !bndry->isDone(); bndry->next1d()) {

          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh));
          }

          f(bndry->x, bndry->y) =
              (8. / 3) * val - 2. * f(bndry->x - bndry->bx, bndry->y - bndry->by)
              + f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by) / 3.;

          // Need to set second guard cell, as may be used for interpolation or upwinding
          // derivatives
          for (int i = 1; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 3.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + f(xi - 3 * bndry->bx, yi - 3 * bndry->by);
          }
        }
      }
    }
  } else {
    // Non-staggered, standard case

    for (; !bndry->isDone(); bndry->next1d()) {

      if (fg) {
        val = fg->generate(Context(bndry, loc, t, mesh));
      }

      f(bndry->x, bndry->y) =
          (8. / 3) * val - 2. * f(bndry->x - bndry->bx, bndry->y - bndry->by)
          + f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by) / 3.;

      // Need to set second guard cell, as may be used for interpolation or upwinding
      // derivatives
      for (int i = 1; i < bndry->width; i++) {
        int xi = bndry->x + i * bndry->bx;
        int yi = bndry->y + i * bndry->by;
        f(xi, yi) = 3.0 * f(xi - bndry->bx, yi - bndry->by)
                    - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                    + f(xi - 3 * bndry->bx, yi - 3 * bndry->by);
      }
    }
  }
}

void BoundaryDirichlet_O3::apply(Field3D& f) { BoundaryDirichlet_O3::apply(f, 0.); }

void BoundaryDirichlet_O3::apply(Field3D& f, BoutReal t) {
  // Set (at 2nd order) the value at the mid-point between the guard cell and the grid
  // cell to be val N.B. Only first guard cells (closest to the grid) should ever be used

  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  BoutReal val = 0.0;

  // Check for staggered grids

  CELL_LOC loc = f.getLocation();
  if (mesh->StaggerGrids && loc != CELL_CENTRE) {
    // Staggered. Need to apply slightly differently

    if (loc == CELL_XLOW) {
      // Field is shifted in X

      if (bndry->bx > 0) {
        // Outer x boundary

        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh));
            }
            f(bndry->x, bndry->y, zk) = val;

            // Need to set second guard cell, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 1; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 3.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->bx < 0) {
        // Inner x boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh));
            }
            f(bndry->x - bndry->bx, bndry->y, zk) = val;

            // Need to set remaining guard cells, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 0; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 3.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->by != 0) {
        // y boundaries

        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh));
            }

            f(bndry->x, bndry->y, zk) =
                (8. / 3) * val - 2. * f(bndry->x - bndry->bx, bndry->y - bndry->by, zk)
                + f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, zk) / 3.;

            // Need to set remaining guard cells, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 1; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 3.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk);
            }
          }
        }
      }
    } else if (loc == CELL_YLOW) {
      // Field is shifted in Y

      if (bndry->by > 0) {
        // Upper y boundary

        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh));
            }
            f(bndry->x, bndry->y, zk) = val;

            // Need to set remaining guard cells, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 1; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 3.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->by < 0) {
        // Lower y boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh));
            }
            f(bndry->x, bndry->y - bndry->by, zk) = val;

            // Need to set remaining guard cells, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 0; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 3.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->bx != 0) {
        // x boundaries
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh));
            }

            f(bndry->x, bndry->y, zk) =
                (8. / 3) * val - 2. * f(bndry->x - bndry->bx, bndry->y - bndry->by, zk)
                + f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, zk) / 3.;

            // Need to set remaining guard cells, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 1; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 3.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk);
            }
          }
        }
      }
    } else if (loc == CELL_ZLOW) {
      // Shifted in Z

      for(; !bndry->isDone(); bndry->next1d()) {
        // Calculate the X and Y normalised values half-way between the guard cell and grid cell
        BoutReal xnorm = 0.5*(   mesh->GlobalX(bndry->x)  // In the guard cell
                                 + mesh->GlobalX(bndry->x - bndry->bx) ); // the grid cell

        BoutReal ynorm = 0.5*(   mesh->GlobalY(bndry->y)  // In the guard cell
                                 + mesh->GlobalY(bndry->y - bndry->by) ); // the grid cell

        for(int zk=0;zk<mesh->LocalNz;zk++) {
          if(fg)
            val = fg->generate(xnorm,TWOPI*ynorm,TWOPI*(zk - 0.5)/(mesh->LocalNz), t);

          f(bndry->x,bndry->y,zk) = (8./3)*val - 2.*f(bndry->x-bndry->bx, bndry->y-bndry->by,zk) + f(bndry->x-2*bndry->bx, bndry->y-2*bndry->by,zk)/3.;

          // Need to set remaining guard cells, as may be used for interpolation or upwinding derivatives
          for(int i=1;i<bndry->width;i++) {
            int xi = bndry->x + i*bndry->bx;
            int yi = bndry->y + i*bndry->by;
            f(xi, yi, zk) = 3.0*f(xi - bndry->bx, yi - bndry->by, zk) - 3.0*f(xi - 2*bndry->bx, yi - 2*bndry->by, zk)
              + f(xi - 3*bndry->bx, yi - 3*bndry->by, zk);
          }
        }
      }
    } else {
      throw BoutException("Unrecognized location");
    }
  } else {
    // Standard (non-staggered) case
    for (; !bndry->isDone(); bndry->next1d()) {
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        if (fg) {
          val = fg->generate(Context(bndry, zk, loc, t, mesh));
        }

        f(bndry->x, bndry->y, zk) =
            (8. / 3) * val - 2. * f(bndry->x - bndry->bx, bndry->y - bndry->by, zk)
            + f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, zk) / 3.;

        // Need to set remaining guard cells, as may be used for interpolation or
        // upwinding derivatives
        for (int i = 1; i < bndry->width; i++) {
          int xi = bndry->x + i * bndry->bx;
          int yi = bndry->y + i * bndry->by;
          f(xi, yi, zk) = 3.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                          - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                          + f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk);
        }
      }
    }
  }
}

void BoundaryDirichlet_O3::apply_ddt(Field2D& f) {
  Field2D* dt = f.timeDeriv();
  for (bndry->first(); !bndry->isDone(); bndry->next())
    (*dt)(bndry->x, bndry->y) = 0.; // Set time derivative to zero
}

void BoundaryDirichlet_O3::apply_ddt(Field3D& f) {
  Mesh* mesh = bndry->localmesh;

  ASSERT1(mesh == f.getMesh());
  Field3D* dt = f.timeDeriv();

  bndry->first();
  for (bndry->first(); !bndry->isDone(); bndry->next()) {
    for (int z = 0; z < mesh->LocalNz; z++) {
      (*dt)(bndry->x, bndry->y, z) = 0.; // Set time derivative to zero
    }
  }
}

///////////////////////////////////////////////////////////////
// Extrapolate to calculate boundary cell to 4th-order

BoundaryOp* BoundaryDirichlet_O4::clone(BoundaryRegion* region,
                                        const std::list<std::string>& args) {
  verifyNumPoints(region, 3);
  std::shared_ptr<FieldGenerator> newgen = nullptr;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryDirichlet_O4(region, newgen);
}

void BoundaryDirichlet_O4::apply(Field2D& f) { BoundaryDirichlet_O4::apply(f, 0.); }

void BoundaryDirichlet_O4::apply(Field2D& f, BoutReal t) {
  // Set (at 2nd order) the value at the mid-point between the guard cell and the grid
  // cell to be val N.B. Only first guard cells (closest to the grid) should ever be used

  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  BoutReal val = 0.0;

  // Check for staggered grids

  CELL_LOC loc = f.getLocation();
  if (mesh->StaggerGrids and (loc == CELL_XLOW or loc == CELL_YLOW)) {
    // Staggered. Need to apply slightly differently

    if (loc == CELL_XLOW) {
      // Field is shifted in X

      if (bndry->bx > 0) {
        // Outer x boundary

        for (; !bndry->isDone(); bndry->next1d()) {
          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh));
          }
          f(bndry->x, bndry->y) = val;

          // Need to set remaining guard cells, as may be used for interpolation or
          // upwinding derivatives
          for (int i = 1; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 4.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 6.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + 4.0 * f(xi - 3 * bndry->bx, yi - 3 * bndry->by)
                        - f(xi - 4 * bndry->bx, yi - 4 * bndry->by);
          }
        }
      }

      if (bndry->bx < 0) {
        // Inner boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {
          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh));
          }

          f(bndry->x - bndry->bx, bndry->y) = val;

          // Need to set remaining guard cells, as may be used for interpolation or
          // upwinding derivatives
          for (int i = 0; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 4.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 6.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + 4.0 * f(xi - 3 * bndry->bx, yi - 3 * bndry->by)
                        - f(xi - 4 * bndry->bx, yi - 4 * bndry->by);
          }
        }
      }
      if (bndry->by != 0) {
        // y boundaries
        for (; !bndry->isDone(); bndry->next1d()) {

          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh));
          }

          f(bndry->x, bndry->y) =
              (16. / 5) * val - 3. * f(bndry->x - bndry->bx, bndry->y - bndry->by)
              + f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by)
              - (1. / 5) * f(bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by);

          // Need to set remaining guard cells, as may be used for interpolation or
          // upwinding derivatives
          for (int i = 1; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 4.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 6.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + 4.0 * f(xi - 3 * bndry->bx, yi - 3 * bndry->by)
                        - f(xi - 4 * bndry->bx, yi - 4 * bndry->by);
          }
        }
      }
    } else if (loc == CELL_YLOW) {
      // Field is shifted in Y

      if (bndry->by > 0) {
        // Outer y boundary

        for (; !bndry->isDone(); bndry->next1d()) {
          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh));
          }
          f(bndry->x, bndry->y) = val;

          // Need to set remaining guard cells, as may be used for interpolation or
          // upwinding derivatives
          for (int i = 1; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 4.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 6.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + 4.0 * f(xi - 3 * bndry->bx, yi - 3 * bndry->by)
                        - f(xi - 4 * bndry->bx, yi - 4 * bndry->by);
          }
        }
      }
      if (bndry->by < 0) {
        // Inner y boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {
          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh));
          }

          f(bndry->x, bndry->y - bndry->by) = val;

          // Need to set remaining guard cells, as may be used for interpolation or
          // upwinding derivatives
          for (int i = 0; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 4.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 6.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + 4.0 * f(xi - 3 * bndry->bx, yi - 3 * bndry->by)
                        - f(xi - 4 * bndry->bx, yi - 4 * bndry->by);
          }
        }
      }
      if (bndry->bx != 0) {
        // x boundaries.

        for (; !bndry->isDone(); bndry->next1d()) {

          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh));
          }

          f(bndry->x, bndry->y) =
              (16. / 5) * val - 3. * f(bndry->x - bndry->bx, bndry->y - bndry->by)
              + f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by)
              - (1. / 5) * f(bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by);

          // Need to set remaining guard cells, as may be used for interpolation or
          // upwinding derivatives
          for (int i = 1; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 4.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 6.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + 4.0 * f(xi - 3 * bndry->bx, yi - 3 * bndry->by)
                        - f(xi - 4 * bndry->bx, yi - 4 * bndry->by);
          }
        }
      }
    }
  } else {
    // Non-staggered, standard case

    for (; !bndry->isDone(); bndry->next1d()) {

      if (fg) {
        val = fg->generate(Context(bndry, loc, t, mesh));
      }

      f(bndry->x, bndry->y) =
          (16. / 5) * val - 3. * f(bndry->x - bndry->bx, bndry->y - bndry->by)
          + f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by)
          - (1. / 5) * f(bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by);

      // Need to set remaining guard cells, as may be used for interpolation or upwinding
      // derivatives
      for (int i = 1; i < bndry->width; i++) {
        int xi = bndry->x + i * bndry->bx;
        int yi = bndry->y + i * bndry->by;
        f(xi, yi) = 4.0 * f(xi - bndry->bx, yi - bndry->by)
                    - 6.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                    + 4.0 * f(xi - 3 * bndry->bx, yi - 3 * bndry->by)
                    - f(xi - 4 * bndry->bx, yi - 4 * bndry->by);
      }
    }
  }
}

void BoundaryDirichlet_O4::apply(Field3D& f) { BoundaryDirichlet_O4::apply(f, 0.); }

void BoundaryDirichlet_O4::apply(Field3D& f, BoutReal t) {
  // Set (at 2nd order) the value at the mid-point between the guard cell and the grid
  // cell to be val N.B. Only first guard cells (closest to the grid) should ever be used

  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  BoutReal val = 0.0;

  // Check for staggered grids

  CELL_LOC loc = f.getLocation();
  if (mesh->StaggerGrids && loc != CELL_CENTRE) {
    // Staggered. Need to apply slightly differently

    if (loc == CELL_XLOW) {
      // Field is shifted in X

      if (bndry->bx > 0) {
        // Outer x boundary

        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh));
            }
            f(bndry->x, bndry->y, zk) = val;

            // Need to set remaining guard cells, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 1; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 4.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 6.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + 4.0 * f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk)
                              - f(xi - 4 * bndry->bx, yi - 4 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->bx < 0) {
        // Inner x boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh));
            }
            f(bndry->x - bndry->bx, bndry->y, zk) = val;

            // Need to set remaining guard cells, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 0; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 4.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 6.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + 4.0 * f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk)
                              - f(xi - 4 * bndry->bx, yi - 4 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->by != 0) {
        // y boundaries

        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh));
            }
            f(bndry->x, bndry->y, zk) =
                (16. / 5) * val - 3. * f(bndry->x - bndry->bx, bndry->y - bndry->by, zk)
                + f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, zk)
                - (1. / 5) * f(bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, zk);

            // Need to set remaining guard cells, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 1; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 4.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 6.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + 4.0 * f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk)
                              - f(xi - 4 * bndry->bx, yi - 4 * bndry->by, zk);
            }
          }
        }
      }
    } else if (loc == CELL_YLOW) {
      // Y boundary, and field is shifted in Y

      if (bndry->by > 0) {
        // Outer y boundary

        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh));
            }

            f(bndry->x, bndry->y, zk) = val;

            // Need to set remaining guard cells, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 1; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 4.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 6.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + 4.0 * f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk)
                              - f(xi - 4 * bndry->bx, yi - 4 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->by < 0) {
        // Inner y boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh));
            }

            f(bndry->x, bndry->y - bndry->by, zk) = val;

            // Need to set remaining guard cells, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 0; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 4.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 6.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + 4.0 * f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk)
                              - f(xi - 4 * bndry->bx, yi - 4 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->bx != 0) {
        // x boundaries

        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh));
            }

            f(bndry->x, bndry->y, zk) =
                (16. / 5) * val - 3. * f(bndry->x - bndry->bx, bndry->y - bndry->by, zk)
                + f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, zk)
                - (1. / 5) * f(bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, zk);

            // Need to set remaining guard cells, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 1; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 4.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 6.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + 4.0 * f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk)
                              - f(xi - 4 * bndry->bx, yi - 4 * bndry->by, zk);
            }
          }
        }
      }
    } else if (loc == CELL_ZLOW) {
      // Shifted in Z
      for(; !bndry->isDone(); bndry->next1d()) {
        // Calculate the X and Y normalised values half-way between the guard cell and grid cell
        BoutReal xnorm = 0.5*(   mesh->GlobalX(bndry->x)  // In the guard cell
                                 + mesh->GlobalX(bndry->x - bndry->bx) ); // the grid cell

        BoutReal ynorm = 0.5*(   mesh->GlobalY(bndry->y)  // In the guard cell
                                 + mesh->GlobalY(bndry->y - bndry->by) ); // the grid cell

        for(int zk=0;zk<mesh->LocalNz;zk++) {
          if(fg)
            val = fg->generate(xnorm,TWOPI*ynorm,TWOPI*(zk - 0.5)/(mesh->LocalNz), t);

          f(bndry->x,bndry->y,zk) = (16./5)*val - 3.*f(bndry->x-bndry->bx, bndry->y-bndry->by,zk) + f(bndry->x-2*bndry->bx, bndry->y-2*bndry->by,zk) - (1./5)*f(bndry->x-3*bndry->bx, bndry->y-3*bndry->by,zk);

          // Need to set remaining guard cells, as may be used for interpolation or upwinding derivatives
          for(int i=1;i<bndry->width;i++) {
            int xi = bndry->x + i*bndry->bx;
            int yi = bndry->y + i*bndry->by;
            f(xi, yi, zk) = 4.0*f(xi - bndry->bx, yi - bndry->by, zk) - 6.0*f(xi - 2*bndry->bx, yi - 2*bndry->by, zk)
              + 4.0*f(xi - 3*bndry->bx, yi - 3*bndry->by, zk) - f(xi - 4*bndry->bx, yi - 4*bndry->by, zk);
          }
        }
      }
    } else {
      throw BoutException("Unrecognized location");
    }
  } else {
    // Standard (non-staggered) case
    for (; !bndry->isDone(); bndry->next1d()) {
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        if (fg) {
          val = fg->generate(Context(bndry, zk, loc, t, mesh));
        }

        f(bndry->x, bndry->y, zk) =
            (16. / 5) * val - 3. * f(bndry->x - bndry->bx, bndry->y - bndry->by, zk)
            + f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, zk)
            - (1. / 5) * f(bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, zk);

        // Need to set remaining guard cells, as may be used for interpolation or
        // upwinding derivatives
        for (int i = 1; i < bndry->width; i++) {
          int xi = bndry->x + i * bndry->bx;
          int yi = bndry->y + i * bndry->by;
          f(xi, yi, zk) = 4.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                          - 6.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                          + 4.0 * f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk)
                          - f(xi - 4 * bndry->bx, yi - 4 * bndry->by, zk);
        }
      }
    }
  }
}

void BoundaryDirichlet_O4::apply_ddt(Field2D& f) {
  Field2D* dt = f.timeDeriv();
  for (bndry->first(); !bndry->isDone(); bndry->next())
    (*dt)(bndry->x, bndry->y) = 0.; // Set time derivative to zero
}

void BoundaryDirichlet_O4::apply_ddt(Field3D& f) {
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  Field3D* dt = f.timeDeriv();
  for (bndry->first(); !bndry->isDone(); bndry->next())
    for (int z = 0; z < mesh->LocalNz; z++)
      (*dt)(bndry->x, bndry->y, z) = 0.; // Set time derivative to zero
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryDirichlet_4thOrder::clone(BoundaryRegion* region,
                                              const std::list<std::string>& args) {
  verifyNumPoints(region, 4);
  if (!args.empty()) {
    // First argument should be a value
    val = stringToReal(args.front());
    return new BoundaryDirichlet_4thOrder(region, val);
  }
  return new BoundaryDirichlet_4thOrder(region);
}

void BoundaryDirichlet_4thOrder::apply(Field2D& f) {
  // Set (at 4th order) the value at the mid-point between the guard cell and the grid
  // cell to be val
  for (bndry->first(); !bndry->isDone(); bndry->next1d()) {
    f(bndry->x, bndry->y) =
        128. / 35. * val - 4. * f(bndry->x - bndry->bx, bndry->y - bndry->by)
        + 2. * f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by)
        - 4. / 3. * f(bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by)
        + 1. / 7. * f(bndry->x - 4 * bndry->bx, bndry->y - 4 * bndry->by);
    f(bndry->x + bndry->bx, bndry->y + bndry->by) =
        -128. / 5. * val + 9. * f(bndry->x, bndry->y)
        + 18. * f(bndry->x - bndry->bx, bndry->y - bndry->by)
        - 4. * f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by)
        + 3. / 5. * f(bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by);
  }
}

void BoundaryDirichlet_4thOrder::apply(Field3D& f) {
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  // Set (at 4th order) the value at the mid-point between the guard cell and the grid
  // cell to be val
  for (bndry->first(); !bndry->isDone(); bndry->next1d())
    for (int z = 0; z < mesh->LocalNz; z++) {
      f(bndry->x, bndry->y, z) =
          128. / 35. * val - 4. * f(bndry->x - bndry->bx, bndry->y - bndry->by, z)
          + 2. * f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, z)
          - 4. / 3. * f(bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, z)
          + 1. / 7. * f(bndry->x - 4 * bndry->bx, bndry->y - 4 * bndry->by, z);
      f(bndry->x + bndry->bx, bndry->y + bndry->by, z) =
          -128. / 5. * val + 9. * f(bndry->x, bndry->y, z)
          + 18. * f(bndry->x - bndry->bx, bndry->y - bndry->by, z)
          - 4. * f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, z)
          + 3. / 5. * f(bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, z);
    }
}

void BoundaryDirichlet_4thOrder::apply_ddt(Field2D& f) {
  Field2D* dt = f.timeDeriv();
  for (bndry->first(); !bndry->isDone(); bndry->next())
    (*dt)(bndry->x, bndry->y) = 0.; // Set time derivative to zero
}

void BoundaryDirichlet_4thOrder::apply_ddt(Field3D& f) {
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  Field3D* dt = f.timeDeriv();
  for (bndry->first(); !bndry->isDone(); bndry->next())
    for (int z = 0; z < mesh->LocalNz; z++)
      (*dt)(bndry->x, bndry->y, z) = 0.; // Set time derivative to zero
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryNeumann_NonOrthogonal::clone(BoundaryRegion* region,
                                                 const std::list<std::string>& args) {
  verifyNumPoints(region, 1);
  if (!args.empty()) {
    output << "WARNING: arguments is set to BoundaryNeumann None Zero Gradient\n";
    // First argument should be a value
    val = stringToReal(args.front());
    return new BoundaryNeumann_NonOrthogonal(region, val);
  }
  return new BoundaryNeumann_NonOrthogonal(region);
}

void BoundaryNeumann_NonOrthogonal::apply(Field2D& f) {
#ifndef COORDINATES_USE_3D
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  Coordinates* metric = f.getCoordinates();

  // Calculate derivatives for metric use
  mesh->communicate(f);
  Field2D dfdy = DDY(f);
  // Loop over all elements and set equal to the next point in
  for (bndry->first(); !bndry->isDone(); bndry->next1d()) {
    // Interpolate (linearly) metrics to halfway between last cell and boundary cell
    BoutReal g11shift =
        0.5
        * (metric->g11(bndry->x, bndry->y) + metric->g11(bndry->x - bndry->bx, bndry->y));
    BoutReal g12shift =
        0.5
        * (metric->g12(bndry->x, bndry->y) + metric->g12(bndry->x - bndry->bx, bndry->y));
    // Have to use derivatives at last gridpoint instead of derivatives on boundary layer
    //   because derivative values don't exist in boundary region
    // NOTE: should be fixed to interpolate to boundary line
    BoutReal xshift = g12shift * dfdy(bndry->x - bndry->bx, bndry->y);

    if (bndry->bx != 0 && bndry->by == 0) {
      // x boundaries only
      BoutReal delta = bndry->bx * metric->dx(bndry->x, bndry->y);
      f(bndry->x, bndry->y) =
          f(bndry->x - bndry->bx, bndry->y) + delta / g11shift * (val - xshift);
      if (bndry->bx == 2) {
        f(bndry->x + bndry->bx, bndry->y) = f(bndry->x - 2 * bndry->bx, bndry->y)
                                            + 3.0 * delta / g11shift * (val - xshift);
      }
    } else if (bndry->by != 0 && bndry->bx == 0) {
      // y boundaries only
      //   no need to shift this b/c we want parallel nuemann not theta
      BoutReal delta = bndry->by * metric->dy(bndry->x, bndry->y);
      f(bndry->x, bndry->y) = f(bndry->x, bndry->y - bndry->by) + delta * val;
      if (bndry->width == 2) {
        f(bndry->x, bndry->y + bndry->by) =
            f(bndry->x, bndry->y - 2 * bndry->by) + 3.0 * delta * val;
      }
    } else {
      // set corners to zero
      f(bndry->x, bndry->y) = 0.0;
      if (bndry->width == 2) {
        f(bndry->x + bndry->bx, bndry->y + bndry->by) = 0.0;
      }
    }
  }
#else
  throw BoutException(
      "Applying boundary to Field2D not compatible with 3D metrics in all cases.");
#endif
}

void BoundaryNeumann_NonOrthogonal::apply(Field3D& f) {
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  Coordinates* metric = f.getCoordinates();

  // Calculate derivatives for metric use
  mesh->communicate(f);
  Field3D dfdy = DDY(f);
  Field3D dfdz = DDZ(f);
  // Loop over all elements and set equal to the next point in
  for (bndry->first(); !bndry->isDone(); bndry->next1d()) {
#ifdef COORDINATES_USE_3D
    for (int z = 0; z < mesh->LocalNz; z++) {
#else
      int z=0;
#endif
    // Interpolate (linearly) metrics to halfway between last cell and boundary cell
    BoutReal g11shift =
        0.5
        * (metric->g11(bndry->x, bndry->y, z) + metric->g11(bndry->x - bndry->bx, bndry->y, z));
    BoutReal g12shift =
        0.5
        * (metric->g12(bndry->x, bndry->y, z) + metric->g12(bndry->x - bndry->bx, bndry->y, z));
    BoutReal g13shift =
        0.5
        * (metric->g13(bndry->x, bndry->y, z) + metric->g13(bndry->x - bndry->bx, bndry->y, z));
    // Have to use derivatives at last gridpoint instead of derivatives on boundary layer
    //   because derivative values don't exist in boundary region
    // NOTE: should be fixed to interpolate to boundary line
#ifndef COORDINATES_USE_3D
    for (int z = 0; z < mesh->LocalNz; z++) {
#endif
      BoutReal xshift = g12shift * dfdy(bndry->x - bndry->bx, bndry->y, z)
                        + g13shift * dfdz(bndry->x - bndry->bx, bndry->y, z);
      if (bndry->bx != 0 && bndry->by == 0) {
        // x boundaries only
        BoutReal delta = bndry->bx * metric->dx(bndry->x, bndry->y, z);
        f(bndry->x, bndry->y, z) =
            f(bndry->x - bndry->bx, bndry->y, z) + delta / g11shift * (val - xshift);
        if (bndry->width == 2) {
          f(bndry->x + bndry->bx, bndry->y, z) =
              f(bndry->x - 2 * bndry->bx, bndry->y, z)
              + 3.0 * delta / g11shift * (val - xshift);
        }
      } else if (bndry->by != 0 && bndry->bx == 0) {
        // y boundaries only
        //   no need to shift this b/c we want parallel nuemann not theta
        BoutReal delta = bndry->by * metric->dy(bndry->x, bndry->y, z);
        f(bndry->x, bndry->y, z) = f(bndry->x, bndry->y - bndry->by, z) + delta * val;
        if (bndry->width == 2) {
          f(bndry->x, bndry->y + bndry->by, z) =
              f(bndry->x, bndry->y - 2 * bndry->by, z) + 3.0 * delta * val;
        }
      } else {
        // set corners to zero
        f(bndry->x, bndry->y, z) = 0.0;
        if (bndry->width == 2) {
          f(bndry->x + bndry->bx, bndry->y + bndry->by, z) = 0.0;
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryNeumann::clone(BoundaryRegion* region,
                                   const std::list<std::string>& args) {
  verifyNumPoints(region, 1);
  std::shared_ptr<FieldGenerator> newgen = nullptr;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryNeumann(region, newgen);
}

void BoundaryNeumann::apply(Field2D& f) { BoundaryNeumann::apply(f, 0.); }

void BoundaryNeumann::apply(Field2D& f, BoutReal t) {
  // Set (at 2nd order / 3rd order) the value at the mid-point between
  // the guard cell and the grid cell to be val
  // N.B. First guard cells (closest to the grid) is 2nd order, while
  // 2nd is 3rd order

#ifndef COORDINATES_USE_3D
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  Coordinates* metric = f.getCoordinates();

  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  BoutReal val = 0.0;

  // Check for staggered grids

  CELL_LOC loc = f.getLocation();
  if (mesh->StaggerGrids and (loc == CELL_XLOW or loc == CELL_YLOW)) {
    // Staggered. Need to apply slightly differently
    // Use one-sided differencing. Cell is now on
    // the boundary, so use one-sided differencing

    if (loc == CELL_XLOW) {
      // Field is shifted in X

      if (bndry->bx > 0) {
        // Outer x boundary

        for (; !bndry->isDone(); bndry->next1d()) {

          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh)) * metric->dx(bndry->x, bndry->y);
          }

          f(bndry->x, bndry->y) = (4. * f(bndry->x - bndry->bx, bndry->y)
                                   - f(bndry->x - 2 * bndry->bx, bndry->y) + 2. * val)
                                  / 3.;
          // Need to set second guard cell, as may be used for interpolation or upwinding
          // derivatives
          for (int i = 1; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            // Use third order extrapolation because boundary point is set to third order,
            // and these points may be used be used by 2nd order upwinding type schemes,
            // which require 3rd order
            f(xi, yi) = 3.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + f(xi - 3 * bndry->bx, yi - 3 * bndry->by);
          }
        }
      }
      if (bndry->bx < 0) {
        // Inner x boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {

          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh)) * metric->dx(bndry->x, bndry->y);
          }

          f(bndry->x - bndry->bx, bndry->y) =
              (4. * f(bndry->x - 2 * bndry->bx, bndry->y)
               - f(bndry->x - 3 * bndry->bx, bndry->y) - 2. * val)
              / 3.;

          // Need to set second guard cell, as may be used for interpolation or upwinding
          // derivatives
          for (int i = 0; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            // Use third order extrapolation because boundary point is set to third order,
            // and these points may be used be used by 2nd order upwinding type schemes,
            // which require 3rd order
            f(xi, yi) = 3.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + f(xi - 3 * bndry->bx, yi - 3 * bndry->by);
          }
        }
      }
      if (bndry->by != 0) {
        // y boundaries

        for (bndry->first(); !bndry->isDone(); bndry->next1d()) {
          BoutReal delta = bndry->bx * metric->dx(bndry->x, bndry->y)
                           + bndry->by * metric->dy(bndry->x, bndry->y);

          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh));
          }

          f(bndry->x, bndry->y) =
              f(bndry->x - bndry->bx, bndry->y - bndry->by) + delta * val;
          if (bndry->width == 2) {
            f(bndry->x + bndry->bx, bndry->y + bndry->by) =
                f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by) + 3.0 * delta * val;
          }
        }
      }
    } else if (loc == CELL_YLOW) {
      // Y boundary, and field is shifted in Y

      if (bndry->by > 0) {
        // Outer y boundary

        for (; !bndry->isDone(); bndry->next1d()) {
          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh)) * metric->dy(bndry->x, bndry->y);
          }
          f(bndry->x, bndry->y) = (4. * f(bndry->x, bndry->y - bndry->by)
                                   - f(bndry->x, bndry->y - 2 * bndry->by) + 2. * val)
                                  / 3.;

          // Need to set second guard cell, as may be used for interpolation or upwinding
          // derivatives
          for (int i = 1; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            // Use third order extrapolation because boundary point is set to third order,
            // and these points may be used be used by 2nd order upwinding type schemes,
            // which require 3rd order
            f(xi, yi) = 3.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + f(xi - 3 * bndry->bx, yi - 3 * bndry->by);
          }
        }
      }
      if (bndry->by < 0) {
        // Inner y boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {

          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh)) * metric->dy(bndry->x, bndry->y - bndry->by);
          }
          f(bndry->x, bndry->y - bndry->by) =
              (4. * f(bndry->x, bndry->y - 2 * bndry->by)
               - f(bndry->x, bndry->y - 3 * bndry->by) - 2. * val)
              / 3.;

          // Need to set second guard cell, as may be used for interpolation or upwinding
          // derivatives
          for (int i = 0; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            // Use third order extrapolation because boundary point is set to third order,
            // and these points may be used be used by 2nd order upwinding type schemes,
            // which require 3rd order
            f(xi, yi) = 3.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + f(xi - 3 * bndry->bx, yi - 3 * bndry->by);
          }
        }
      }
      if (bndry->bx != 0) {
        // x boundaries
        for (bndry->first(); !bndry->isDone(); bndry->next1d()) {
          BoutReal delta = bndry->bx * metric->dx(bndry->x, bndry->y)
                           + bndry->by * metric->dy(bndry->x, bndry->y);

          if (fg) {
            val = fg->generate(Context(bndry, loc, t, mesh));
          }

          f(bndry->x, bndry->y) =
              f(bndry->x - bndry->bx, bndry->y - bndry->by) + delta * val;
          if (bndry->width == 2) {
            f(bndry->x + bndry->bx, bndry->y + bndry->by) =
                f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by) + 3.0 * delta * val;
          }
        }
      }
    }
  } else {
    // Non-staggered, standard case

    for (bndry->first(); !bndry->isDone(); bndry->next1d()) {
      BoutReal delta = bndry->bx * metric->dx(bndry->x, bndry->y)
                       + bndry->by * metric->dy(bndry->x, bndry->y);

      if (fg) {
        val = fg->generate(Context(bndry, loc, t, mesh));
      }

      f(bndry->x, bndry->y) = f(bndry->x - bndry->bx, bndry->y - bndry->by) + delta * val;
      if (bndry->width == 2) {
        f(bndry->x + bndry->bx, bndry->y + bndry->by) =
            f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by) + 3.0 * delta * val;
      }
    }
  }
#else
  throw BoutException(
      "Applying boundary to Field2D not compatible with 3D metrics in all cases.");
#endif
}

void BoundaryNeumann::apply(Field3D& f) { BoundaryNeumann::apply(f, 0.); }

void BoundaryNeumann::apply(Field3D& f, BoutReal t) {
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  Coordinates* metric = f.getCoordinates();

  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  BoutReal val = 0.0;

  // Check for staggered grids

  CELL_LOC loc = f.getLocation();
  if (mesh->StaggerGrids && loc != CELL_CENTRE) {
    // Staggered. Need to apply slightly differently
    // Use one-sided differencing. Cell is now on
    // the boundary, so use one-sided differencing

    if (loc == CELL_XLOW) {
      // Field is shifted in X

      if (bndry->bx > 0) {
        // Outer x boundary
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh)) * metric->dx(bndry->x, bndry->y, zk);
            }

            f(bndry->x, bndry->y, zk) =
                (4. * f(bndry->x - bndry->bx, bndry->y, zk)
                 - f(bndry->x - 2 * bndry->bx, bndry->y, zk) + 2. * val)
                / 3.;

            // Need to set second guard cell, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 1; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              // Use third order extrapolation because boundary point is set to third
              // order, and these points may be used be used by 2nd order upwinding type
              // schemes, which require 3rd order
              f(xi, yi, zk) = 3.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->bx < 0) {
        // Inner x boundary
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh)) * metric->dx(bndry->x - bndry->bx, bndry->y, zk);
            }

            f(bndry->x - bndry->bx, bndry->y, zk) =
                (4. * f(bndry->x - 2 * bndry->bx, bndry->y, zk)
                 - f(bndry->x - 3 * bndry->bx, bndry->y, zk) - 2. * val)
                / 3.;

            // Need to set second guard cell, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 0; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              // Use third order extrapolation because boundary point is set to third
              // order, and these points may be used be used by 2nd order upwinding type
              // schemes, which require 3rd order
              f(xi, yi, zk) = 3.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->by != 0) {
        for (; !bndry->isDone(); bndry->next1d()) {
#ifndef COORDINATES_USE_3D
          BoutReal delta = bndry->bx * metric->dx(bndry->x, bndry->y)
                           + bndry->by * metric->dy(bndry->x, bndry->y);
#endif
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
#ifdef COORDINATES_USE_3D
            BoutReal delta = bndry->bx * metric->dx(bndry->x, bndry->y, zk)
                             + bndry->by * metric->dy(bndry->x, bndry->y, zk);
#endif
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh));
            }
            f(bndry->x, bndry->y, zk) =
                f(bndry->x - bndry->bx, bndry->y - bndry->by, zk) + delta * val;
            if (bndry->width == 2) {
              f(bndry->x + bndry->bx, bndry->y + bndry->by, zk) =
                  f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, zk)
                  + 3.0 * delta * val;
            }
          }
        }
      }
    } else if (loc == CELL_YLOW) {
      // Field is shifted in Y

      if (bndry->by > 0) {
        // Outer y boundary
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh)) * metric->dy(bndry->x, bndry->y, zk);
            }
            f(bndry->x, bndry->y, zk) =
                (4. * f(bndry->x, bndry->y - bndry->by, zk)
                 - f(bndry->x, bndry->y - 2 * bndry->by, zk) + 2. * val)
                / 3.;

            // Need to set second guard cell, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 1; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              // Use third order extrapolation because boundary point is set to third
              // order, and these points may be used be used by 2nd order upwinding type
              // schemes, which require 3rd order
              f(xi, yi, zk) = 3.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->by < 0) {
        // Inner y boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh)) * metric->dy(bndry->x, bndry->y - bndry->by, zk);
            }

            f(bndry->x, bndry->y - bndry->by, zk) =
                (4. * f(bndry->x, bndry->y - 2 * bndry->by, zk)
                 - f(bndry->x, bndry->y - 3 * bndry->by, zk) - 2. * val)
                / 3.;

            // Need to set second guard cell, as may be used for interpolation or
            // upwinding derivatives
            for (int i = 0; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              // Use third order extrapolation because boundary point is set to third
              // order, and these points may be used be used by 2nd order upwinding type
              // schemes, which require 3rd order
              f(xi, yi, zk) = 3.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->bx != 0) {
        // x boundaries.
        for (; !bndry->isDone(); bndry->next1d()) {
#ifndef COORDINATES_USE_3D
	  int zk=0;
          BoutReal delta = bndry->bx * metric->dx(bndry->x, bndry->y, zk)
                         + bndry->by * metric->dy(bndry->x, bndry->y, zk);
#endif
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
#ifdef COORDINATES_USE_3D
          BoutReal delta = bndry->bx * metric->dx(bndry->x, bndry->y, zk)
                         + bndry->by * metric->dy(bndry->x, bndry->y, zk);
#endif
            if (fg) {
              val = fg->generate(Context(bndry, zk, loc, t, mesh));
            }
            f(bndry->x, bndry->y, zk) =
                f(bndry->x - bndry->bx, bndry->y - bndry->by, zk) + delta * val;
            if (bndry->width == 2) {
              f(bndry->x + bndry->bx, bndry->y + bndry->by, zk) =
                  f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, zk)
                  + 3.0 * delta * val;
            }
          }
        }
      }
    } else if (loc == CELL_ZLOW) {
      // Shifted in Z
      for(; !bndry->isDone(); bndry->next1d()) {
        // Calculate the X and Y normalised values half-way between the guard cell and grid cell
        BoutReal xnorm = 0.5*(   mesh->GlobalX(bndry->x)  // In the guard cell
                                 + mesh->GlobalX(bndry->x - bndry->bx) ); // the grid cell

        BoutReal ynorm = 0.5*(   mesh->GlobalY(bndry->y)  // In the guard cell
                                 + mesh->GlobalY(bndry->y - bndry->by) ); // the grid cell


        for(int zk=0;zk<mesh->LocalNz;zk++) {
	  BoutReal delta = bndry->bx*metric->dx(bndry->x,bndry->y, zk)+bndry->by*metric->dy(bndry->x,bndry->y, zk);
          if(fg){
            val = fg->generate(xnorm,TWOPI*ynorm,TWOPI*(zk - 0.5)/(mesh->LocalNz),t);
          }
          f(bndry->x,bndry->y, zk) = f(bndry->x-bndry->bx, bndry->y-bndry->by, zk) + delta*val;
          if (bndry->width == 2){
            f(bndry->x + bndry->bx, bndry->y + bndry->by, zk) = f(bndry->x - 2*bndry->bx, bndry->y - 2*bndry->by, zk) + 3.0*delta*val;
          }
        }
      }
    } else {
      throw BoutException("Unrecognized location");
    }
  } else {
    for (; !bndry->isDone(); bndry->next1d()) {
#ifdef COORDINATES_USE_3D
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
	BoutReal delta = bndry->bx * metric->dx(bndry->x, bndry->y, zk)
	  + bndry->by * metric->dy(bndry->x, bndry->y, zk);
#else
      BoutReal delta = bndry->bx * metric->dx(bndry->x, bndry->y)
                     + bndry->by * metric->dy(bndry->x, bndry->y);
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
#endif
        if (fg) {
          val = fg->generate(Context(bndry, zk, loc, t, mesh));
        }
        f(bndry->x, bndry->y, zk) =
            f(bndry->x - bndry->bx, bndry->y - bndry->by, zk) + delta * val;
        if (bndry->width == 2) {
          f(bndry->x + bndry->bx, bndry->y + bndry->by, zk) =
              f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, zk)
              + 3.0 * delta * val;
        }
      }
    }
  }
}

void BoundaryNeumann::apply_ddt(Field2D& f) {
  Field2D* dt = f.timeDeriv();
  for (bndry->first(); !bndry->isDone(); bndry->next())
    (*dt)(bndry->x, bndry->y) = 0.; // Set time derivative to zero
}

void BoundaryNeumann::apply_ddt(Field3D& f) {
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  Field3D* dt = f.timeDeriv();
  for (bndry->first(); !bndry->isDone(); bndry->next())
    for (int z = 0; z < mesh->LocalNz; z++)
      (*dt)(bndry->x, bndry->y, z) = 0.; // Set time derivative to zero
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryNeumann_O4::clone(BoundaryRegion* region,
                                      const std::list<std::string>& args) {
  std::shared_ptr<FieldGenerator> newgen = nullptr;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryNeumann_O4(region, newgen);
}

void BoundaryNeumann_O4::apply(Field2D& f) { BoundaryNeumann_O4::apply(f, 0.); }

void BoundaryNeumann_O4::apply(Field2D& f, BoutReal t) {
#ifndef COORDINATES_USE_3D
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());

  // Set (at 4th order) the value at the mid-point between the guard cell and the grid
  // cell to be val N.B. Only first guard cells (closest to the grid) should ever be used
  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  BoutReal val = 0.0;

  // Check for staggered grids
  CELL_LOC loc = f.getLocation();
  if (mesh->StaggerGrids && loc != CELL_CENTRE) {
    throw BoutException("neumann_o4 not implemented with staggered grid yet");
  } else {
    // Non-staggered, standard case

    Coordinates* coords = f.getCoordinates();

    for (bndry->first(); !bndry->isDone(); bndry->next1d()) {
      BoutReal delta = bndry->bx * coords->dx(bndry->x, bndry->y)
                       + bndry->by * coords->dy(bndry->x, bndry->y);

      if (fg) {
        val = fg->generate(Context(bndry, loc, t, mesh));
      }

      f(bndry->x, bndry->y) =
          12. * delta * val / 11.
          + (+17. * f(bndry->x - bndry->bx, bndry->y - bndry->by)
             + 9. * f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by)
             - 5. * f(bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by)
             + f(bndry->x - 4 * bndry->bx, bndry->y - 4 * bndry->by))
                / 22.;

      if (bndry->width == 2) {
        throw BoutException("neumann_o4 with a boundary width of 2 not implemented yet");
      }
    }
  }
#else
  throw BoutException(
      "Applying boundary to Field2D not compatible with 3D metrics in all cases.");
#endif
}

void BoundaryNeumann_O4::apply(Field3D& f) { BoundaryNeumann_O4::apply(f, 0.); }

void BoundaryNeumann_O4::apply(Field3D& f, BoutReal t) {
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  BoutReal val = 0.0;

  // Check for staggered grids
  CELL_LOC loc = f.getLocation();
  if (mesh->StaggerGrids && loc != CELL_CENTRE) {
    throw BoutException("neumann_o4 not implemented with staggered grid yet");
  } else {
    Coordinates* coords = f.getCoordinates();
    for (; !bndry->isDone(); bndry->next1d()) {
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
	BoutReal delta = bndry->bx * coords->dx(bndry->x, bndry->y, zk)
                       + bndry->by * coords->dy(bndry->x, bndry->y, zk);
        if (fg) {
          val = fg->generate(Context(bndry, zk, loc, t, mesh));
        }

        f(bndry->x, bndry->y, zk) =
            12. * delta * val / 11.
            + (+17. * f(bndry->x - bndry->bx, bndry->y - bndry->by, zk)
               + 9. * f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, zk)
               - 5. * f(bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, zk)
               + f(bndry->x - 4 * bndry->bx, bndry->y - 4 * bndry->by, zk))
                  / 22.;

        if (bndry->width == 2) {
          throw BoutException(
              "neumann_o4 with a boundary width of 2 not implemented yet");
        }
      }
    }
  }
}

void BoundaryNeumann_O4::apply_ddt(Field2D& f) {
  Field2D* dt = f.timeDeriv();
  for (bndry->first(); !bndry->isDone(); bndry->next())
    (*dt)(bndry->x, bndry->y) = 0.; // Set time derivative to zero
}

void BoundaryNeumann_O4::apply_ddt(Field3D& f) {
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  Field3D *dt = f.timeDeriv();
  for (bndry->first(); !bndry->isDone(); bndry->next())
    for (int z = 0; z < mesh->LocalNz; z++)
      (*dt)(bndry->x, bndry->y, z) = 0.; // Set time derivative to zero
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryNeumann_4thOrder::clone(BoundaryRegion* region,
                                            const std::list<std::string>& args) {
  verifyNumPoints(region, 4);
  if (!args.empty()) {
    // First argument should be a value
    val = stringToReal(args.front());
    return new BoundaryNeumann_4thOrder(region, val);
  }
  return new BoundaryNeumann_4thOrder(region);
}

void BoundaryNeumann_4thOrder::apply(Field2D& f) {
#ifndef COORDINATES_USE_3D
  Coordinates* metric = f.getCoordinates();
  // Set (at 4th order) the gradient at the mid-point between the guard cell and the grid
  // cell to be val This sets the value of the co-ordinate derivative, i.e. DDX/DDY not
  // Grad_par/Grad_perp.x
  for (bndry->first(); !bndry->isDone(); bndry->next1d()) {
    BoutReal delta = -(bndry->bx * metric->dx(bndry->x, bndry->y)
                       + bndry->by * metric->dy(bndry->x, bndry->y));
    f(bndry->x, bndry->y) =
        12. * delta / 11. * val
        + 17. / 22. * f(bndry->x - bndry->bx, bndry->y - bndry->by)
        + 9. / 22. * f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by)
        - 5. / 22. * f(bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by)
        + 1. / 22. * f(bndry->x - 4 * bndry->bx, bndry->y - 4 * bndry->by);
    f(bndry->x + bndry->bx, bndry->y + bndry->by) =
        -24. * delta * val + 27. * f(bndry->x, bndry->y)
        - 27. * f(bndry->x - bndry->bx, bndry->y - bndry->by)
        + f(bndry->x - 2 * bndry->bx,
            bndry->y
                - 2 * bndry->by); // The f(bndry->x-4*bndry->bx,bndry->y-4*bndry->by) term
                                  // vanishes, so that this sets to zero the 4th order
                                  // central difference first derivative at the point half
                                  // way between the guard cell and the grid cell
  }
#else
  throw BoutException("void BoundaryNeumann_4thOrder::apply(Field2D& f) not implemented with 3D metrics");
#endif
}

void BoundaryNeumann_4thOrder::apply(Field3D& f) {
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  Coordinates *metric = f.getCoordinates();
  // Set (at 4th order) the gradient at the mid-point between the guard cell and the grid cell to be val
  // This sets the value of the co-ordinate derivative, i.e. DDX/DDY not Grad_par/Grad_perp.x
  for (bndry->first(); !bndry->isDone(); bndry->next1d())
    for (int z = 0;z < mesh->LocalNz; z++) {
      BoutReal delta = -(bndry->bx*metric->dx(bndry->x,bndry->y, z)+bndry->by*metric->dy(bndry->x,bndry->y, z));
      f(bndry->x,bndry->y,z) = 12.*delta/11.*val + 17./22.*f(bndry->x-bndry->bx,bndry->y-bndry->by,z) + 9./22.*f(bndry->x-2*bndry->bx,bndry->y-2*bndry->by,z) - 5./22.*f(bndry->x-3*bndry->bx,bndry->y-3*bndry->by,z) + 1./22.*f(bndry->x-4*bndry->bx,bndry->y-4*bndry->by,z);
      f(bndry->x+bndry->bx,bndry->y+bndry->by,z) = -24.*delta*val + 27.*f(bndry->x,bndry->y,z) - 27.*f(bndry->x-bndry->bx,bndry->y-bndry->by,z) + f(bndry->x-2*bndry->bx,bndry->y-2*bndry->by,z); // The f(bndry->x-4*bndry->bx,bndry->y-4*bndry->by,z) term vanishes, so that this sets to zero the 4th order central difference first derivative at the point half way between the guard cell and the grid cell
    }
}

void BoundaryNeumann_4thOrder::apply_ddt(Field2D& f) {
  Field2D* dt = f.timeDeriv();
  for (bndry->first(); !bndry->isDone(); bndry->next())
    (*dt)(bndry->x, bndry->y) = 0.; // Set time derivative to zero
}

void BoundaryNeumann_4thOrder::apply_ddt(Field3D& f) {
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  Field3D *dt = f.timeDeriv();
  for (bndry->first(); !bndry->isDone(); bndry->next())
    for (int z = 0;z < mesh->LocalNz; z++)
      (*dt)(bndry->x,bndry->y,z) = 0.; // Set time derivative to zero
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryNeumannPar::clone(BoundaryRegion* region,
                                      const std::list<std::string>& args) {
  verifyNumPoints(region, 1);
  if (!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryNeumann2\n";
  }
  return new BoundaryNeumannPar(region);
}

void BoundaryNeumannPar::apply(Field2D& f) {
#ifndef COORDINATES_USE_3D
  Coordinates* metric = f.getCoordinates();
  // Loop over all elements and set equal to the next point in
  for (bndry->first(); !bndry->isDone(); bndry->next())
    f(bndry->x, bndry->y) =
        f(bndry->x - bndry->bx, bndry->y - bndry->by)
        * sqrt(metric->g_22(bndry->x, bndry->y)
               / metric->g_22(bndry->x - bndry->bx, bndry->y - bndry->by));
#else
  throw BoutException(
          "Applying boundary to Field2D not compatible with 3D metrics in all cases.");
#endif
}

void BoundaryNeumannPar::apply(Field3D& f) {
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  Coordinates *metric = f.getCoordinates();
  for (bndry->first(); !bndry->isDone(); bndry->next())
    for (int z = 0;z < mesh->LocalNz; z++)
      f(bndry->x,bndry->y,z) = f(bndry->x - bndry->bx,bndry->y - bndry->by,z)*sqrt(metric->g_22(bndry->x, bndry->y, z)/metric->g_22(bndry->x - bndry->bx, bndry->y - bndry->by,z));
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryRobin::clone(BoundaryRegion* region,
                                 const std::list<std::string>& args) {
  verifyNumPoints(region, 1);
  BoutReal a = 0.5, b = 1.0, g = 0.;

  auto it = args.begin();

  if (it != args.end()) {
    // First argument is 'a'
    a = stringToReal(*it);
    it++;

    if (it != args.end()) {
      // Second is 'b'
      b = stringToReal(*it);
      it++;

      if (it != args.end()) {
        // Third is 'g'
        g = stringToReal(*it);
        it++;
        if (it != args.end()) {
          output
              << "WARNING: BoundaryRobin takes maximum of 3 arguments. Ignoring extras\n";
        }
      }
    }
  }

  return new BoundaryRobin(region, a, b, g);
}

void BoundaryRobin::apply(Field2D& f) {
  if (fabs(bval) < 1.e-12) {
    // No derivative term so just constant value
    for (bndry->first(); !bndry->isDone(); bndry->next())
      f(bndry->x, bndry->y) = gval / aval;
  } else {
    BoutReal sign = 1.;
    if ((bndry->bx < 0) || (bndry->by < 0))
      sign = -1.;
    for (bndry->first(); !bndry->isDone(); bndry->next())
      f(bndry->x, bndry->y) =
          f(bndry->x - bndry->bx, bndry->y - bndry->by)
          + sign * (gval - aval * f(bndry->x - bndry->bx, bndry->y - bndry->by)) / bval;
  }
}

void BoundaryRobin::apply(Field3D& f) {
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  if (fabs(bval) < 1.e-12) {
    for (bndry->first(); !bndry->isDone(); bndry->next())
      for (int z = 0;z < mesh->LocalNz; z++)
	f(bndry->x, bndry->y, z) = gval / aval;
  } else {
    BoutReal sign = 1.;
    if ((bndry->bx < 0) || (bndry->by < 0))
      sign = -1.;
    for (bndry->first(); !bndry->isDone(); bndry->next())
      for (int z = 0; z < mesh->LocalNz; z++)
        f(bndry->x, bndry->y, z) =
            f(bndry->x - bndry->bx, bndry->y - bndry->by, z)
            + sign * (gval - aval * f(bndry->x - bndry->bx, bndry->y - bndry->by, z))
                  / bval;
  }
}

///////////////////////////////////////////////////////////////

void BoundaryConstGradient::apply(Field2D& f) {
  // Loop over all elements and set equal to the next point in
  for (bndry->first(); !bndry->isDone(); bndry->next())
    f(bndry->x, bndry->y) = 2. * f(bndry->x - bndry->bx, bndry->y - bndry->by)
                            - f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by);
}

void BoundaryConstGradient::apply(Field3D& f) {
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  for (bndry->first(); !bndry->isDone(); bndry->next())
    for (int z = 0;z < mesh->LocalNz; z++)
      f(bndry->x, bndry->y, z) = 2.*f(bndry->x - bndry->bx, bndry->y - bndry->by, z) - f(bndry->x - 2*bndry->bx,bndry->y - 2*bndry->by,z);
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryConstGradient::clone(BoundaryRegion* region,
                                         const std::list<std::string>& args) {
  verifyNumPoints(region, 2);
  if (!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryConstGradient\n";
  }
  return new BoundaryConstGradient(region);
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryZeroLaplace::clone(BoundaryRegion* region,
                                       const std::list<std::string>& args) {
  verifyNumPoints(region, 2);
  if (!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryZeroLaplace\n";
  }
  return new BoundaryZeroLaplace(region);
}

void BoundaryZeroLaplace::apply(Field2D& f) {
#ifndef COORDINATES_USE_3D
  Coordinates* metric = f.getCoordinates();
  if ((bndry->location != BNDRY_XIN) && (bndry->location != BNDRY_XOUT)) {
    // Can't apply this boundary condition to non-X boundaries
    throw BoutException(
        "ERROR: Can't apply Zero Laplace condition to non-X boundaries\n");
  }
  // Constant X derivative
  int bx = bndry->bx;
  // Loop over the Y dimension
  for (bndry->first(); !bndry->isDone(); bndry->nextY()) {
    int x = bndry->x;
    int y = bndry->y;
    BoutReal g = (f(x - bx, y) - f(x - 2 * bx, y)) / metric->dx(x - bx, y);
    // Loop in X towards edge of domain
    do {
      f(x, y) = f(x - bx, y) + g * metric->dx(x, y);
      bndry->nextX();
      x = bndry->x;
      y = bndry->y;
    } while (!bndry->isDone());
  }
#else
  throw BoutException(
      "Applying boundary to Field2D not compatible with 3D metrics in all cases.");
#endif
}

void BoundaryZeroLaplace::apply(Field3D& f) {
#ifndef COORDINATES_USE_3D
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  int ncz = mesh->LocalNz;

  Coordinates* metric = f.getCoordinates();

  Array<dcomplex> c0(ncz / 2 + 1);
  Array<dcomplex> c1(ncz / 2 + 1);

  if ((bndry->location != BNDRY_XIN) && (bndry->location != BNDRY_XOUT)) {
    // Can't apply this boundary condition to non-X boundaries
    throw BoutException(
        "ERROR: Can't apply Zero Laplace condition to non-X boundaries\n");
  }

  int bx = bndry->bx;
  // Loop over the Y dimension
  for (bndry->first(); !bndry->isDone(); bndry->nextY()) {
    // bndry->(x,y) is the first point in the boundary
    // bndry->(x-bx,y) is the last "real" point in the domain

    int x = bndry->x;
    int y = bndry->y;

    // Take FFT of last 2 points in domain
    rfft(f(x - bx, y), mesh->LocalNz, c0.begin());
    rfft(f(x - 2 * bx, y), mesh->LocalNz, c1.begin());
    c1[0] = c0[0] - c1[0]; // Only need gradient

    // Solve  metric->g11*d2f/dx2 - metric->g33*kz^2f = 0
    // Assume metric->g11, metric->g33 constant -> exponential growth or decay

    // Loop in X towards edge of domain
    do {
      // kz = 0 solution
      c0[0] += c1[0]; // Straight line

      // kz != 0 solution
      BoutReal coef =
          -1.0 * sqrt(metric->g33(x, y) / metric->g11(x, y)) * metric->dx(x, y);
      for (int jz = 1; jz <= ncz / 2; jz++) {
        BoutReal kwave = jz * 2.0 * PI / metric->zlength()(x,y); // wavenumber in [rad^-1]
        c0[jz] *= exp(coef * kwave);                        // The decaying solution only
      }
      // Reverse FFT
      irfft(c0.begin(), mesh->LocalNz, f(x, y));

      bndry->nextX();
      x = bndry->x;
      y = bndry->y;
    } while (!bndry->isDone());
  }
#else
  throw BoutException(
      "Applying boundary to Field3D not compatible with 3D metrics in ZeroLaplace case.");
#endif
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryZeroLaplace2::clone(BoundaryRegion* region,
                                        const std::list<std::string>& args) {
  verifyNumPoints(region, 3);
  if (!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryZeroLaplace2\n";
  }
  return new BoundaryZeroLaplace2(region);
}

void BoundaryZeroLaplace2::apply(Field2D& f) {
#ifndef COORDINATES_USE_3D
  if ((bndry->location != BNDRY_XIN) && (bndry->location != BNDRY_XOUT)) {
    // Can't apply this boundary condition to non-X boundaries
    throw BoutException(
        "ERROR: Can't apply Zero Laplace condition to non-X boundaries\n");
  }

  Coordinates* metric = f.getCoordinates();

  // Constant X derivative
  int bx = bndry->bx;
  // Loop over the Y dimension
  for (bndry->first(); !bndry->isDone(); bndry->nextY()) {
    int x = bndry->x;
    int y = bndry->y;
    BoutReal g = (f(x - bx, y) - f(x - 2 * bx, y)) / metric->dx(x - bx, y);
    // Loop in X towards edge of domain
    do {
      f(x, y) = f(x - bx, y) + g * metric->dx(x, y);
      bndry->nextX();
      x = bndry->x;
      y = bndry->y;
    } while (!bndry->isDone());
  }
#else
  throw BoutException(
      "Applying boundary to Field2D not compatible with 3D metrics in all cases.");
#endif
}

void BoundaryZeroLaplace2::apply(Field3D& f) {
#ifndef COORDINATES_USE_3D
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  int ncz = mesh->LocalNz;

  ASSERT0(ncz % 2 == 0); // Allocation assumes even number

  // allocate memory
  Array<dcomplex> c0(ncz / 2 + 1), c1(ncz / 2 + 1), c2(ncz / 2 + 1);

  if ((bndry->location != BNDRY_XIN) && (bndry->location != BNDRY_XOUT)) {
    // Can't apply this boundary condition to non-X boundaries
    throw BoutException(
        "ERROR: Can't apply Zero Laplace condition to non-X boundaries\n");
  }

  int bx = bndry->bx;
  // Loop over the Y dimension
  for (bndry->first(); !bndry->isDone(); bndry->nextY()) {
    // bndry->(x,y) is the first point in the boundary
    // bndry->(x-bx,y) is the last "real" point in the domain

    int x = bndry->x;
    int y = bndry->y;

    // Take FFT of last 2 points in domain
    rfft(f(x - bx, y), ncz, c1.begin());
    rfft(f(x - 2 * bx, y), ncz, c2.begin());

    // Loop in X towards edge of domain
    do {
      for (int jz = 0; jz <= ncz / 2; jz++) {
        dcomplex la, lb, lc;
        laplace_tridag_coefs(x - bx, y, jz, la, lb, lc);
        if (bx > 0) {
          // Outer boundary
          swap(la, lc);
        }
        c0[jz] = -(lb * c1[jz] + lc * c2[jz]) / la;
      }
      // Reverse FFT
      irfft(c0.begin(), ncz, f(x, y));
      // cycle c0 -> c1 -> c2 -> c0
      swap(c0, c2);
      swap(c2, c1);

      bndry->nextX();
      x = bndry->x;
      y = bndry->y;
    } while (!bndry->isDone());
  }
#else
  throw BoutException("Applying boundary to Field3D not compatible with 3D metrics in "
                      "ZeroLaplace2 case.");
#endif
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryConstLaplace::clone(BoundaryRegion* region,
                                        const std::list<std::string>& args) {
  verifyNumPoints(region, 2);
  if (!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryConstLaplace\n";
  }
  return new BoundaryConstLaplace(region);
}

void BoundaryConstLaplace::apply(Field2D& f) {
  if ((bndry->location != BNDRY_XIN) && (bndry->location != BNDRY_XOUT)) {
    // Can't apply this boundary condition to non-X boundaries
    throw BoutException(
        "ERROR: Can't apply Zero Laplace condition to non-X boundaries\n");
  }

  // Constant X second derivative
  int bx = bndry->bx;
  // Loop over the Y dimension
  for (bndry->first(); !bndry->isDone(); bndry->nextY()) {
    int x = bndry->x;
    int y = bndry->y;
    // Calculate the Laplacian on the last point
    dcomplex la, lb, lc;
    laplace_tridag_coefs(x - 2 * bx, y, 0, la, lb, lc);
    dcomplex val =
        la * f(x - bx - 1, y) + lb * f(x - 2 * bx, y) + lc * f(x - 2 * bx + 1, y);
    // Loop in X towards edge of domain
    do {
      laplace_tridag_coefs(x - bx, y, 0, la, lb, lc);
      if (bx < 0) { // Lower X
        f(x, y) = ((val - lb * f(x - bx, y) + lc * f(x - 2 * bx, y)) / la).real();
      } else // Upper X
        f(x, y) = ((val - lb * f(x - bx, y) + la * f(x - 2 * bx, y)) / lc).real();

      bndry->nextX();
      x = bndry->x;
      y = bndry->y;
    } while (!bndry->isDone());
  }
}

void BoundaryConstLaplace::apply(Field3D& f) {
#ifndef COORDINATES_USE_3D
  if ((bndry->location != BNDRY_XIN) && (bndry->location != BNDRY_XOUT)) {
    // Can't apply this boundary condition to non-X boundaries
    throw BoutException(
        "ERROR: Can't apply Zero Laplace condition to non-X boundaries\n");
  }

  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  Coordinates *metric = f.getCoordinates();
  
  int ncz = mesh->LocalNz;

  // Allocate memory
  Array<dcomplex> c0(ncz / 2 + 1), c1(ncz / 2 + 1), c2(ncz / 2 + 1);

  int bx = bndry->bx;
  // Loop over the Y dimension
  for (bndry->first(); !bndry->isDone(); bndry->nextY()) {
    int x = bndry->x;
    int y = bndry->y;

    // Take FFT of last 3 points in domain
    rfft(f(x - bx, y), ncz, c0.begin());
    rfft(f(x - 2 * bx, y), ncz, c1.begin());
    rfft(f(x - 3 * bx, y), ncz, c2.begin());
    dcomplex k0lin = (c1[0] - c0[0]) / metric->dx(x - bx, y); // for kz=0 solution

    // Calculate Delp2 on point MXG+1 (and put into c1)
    for (int jz = 0; jz <= ncz / 2; jz++) {
      dcomplex la, lb, lc;
      laplace_tridag_coefs(x - 2 * bx, y, jz, la, lb, lc);
      if (bx < 0) { // Inner X
        c1[jz] = la * c0[jz] + lb * c1[jz] + lc * c2[jz];
      } else { // Outer X
        c1[jz] = la * c2[jz] + lb * c1[jz] + lc * c0[jz];
      }
    }
    // Solve  metric->g11*d2f/dx2 - metric->g33*kz^2f = 0
    // Assume metric->g11, metric->g33 constant -> exponential growth or decay
    BoutReal xpos = 0.0;
    // Loop in X towards edge of domain
    do {
      // kz = 0 solution
      xpos -= metric->dx(x, y);
      c2[0] = c0[0] + k0lin * xpos + 0.5 * c1[0] * xpos * xpos / metric->g11(x - bx, y);
      // kz != 0 solution
      BoutReal coef = -1.0 * sqrt(metric->g33(x - bx, y) / metric->g11(x - bx, y))
                      * metric->dx(x - bx, y);
      for (int jz = 1; jz <= ncz / 2; jz++) {
#warning TODO: fix
        BoutReal kwave = jz * 2.0 * PI / metric->zlength()(0,0); // wavenumber in [rad^-1]
        c0[jz] *= exp(coef * kwave);                        // The decaying solution only
        // Add the particular solution
        c2[jz] = c0[jz] - c1[jz] / (metric->g33(x - bx, y) * kwave * kwave);
      }
      // Reverse FFT
      irfft(c2.begin(), ncz, f(x, y));

      bndry->nextX();
      x = bndry->x;
      y = bndry->y;
    } while (!bndry->isDone());
  }
#else
  throw BoutException("Applying boundary to Field3D not compatible with 3D metrics in "
                      "ConstLaplace case.");
#endif
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryDivCurl::clone(BoundaryRegion* region,
                                   const std::list<std::string>& args) {
  if (!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryDivCurl\n";
  }
  return new BoundaryDivCurl(region);
}

void BoundaryDivCurl::apply(Vector2D& UNUSED(f)) {
  throw BoutException("ERROR: DivCurl boundary not yet implemented for 2D vectors\n");
}

void BoundaryDivCurl::apply(Vector3D& var) {
#ifndef COORDINATES_USE_3D
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == var.x.getMesh());

  int jx, jy, jz, jzp, jzm;
  BoutReal tmp;

  Coordinates* metric = mesh->getCoordinates(var.getLocation());

  int ncz = mesh->LocalNz;

  if (bndry->location != BNDRY_XOUT) {
    throw BoutException("ERROR: DivCurl boundary only works for outer X currently\n");
  }
  var.toCovariant();

  if (mesh->xstart > 2) {
    throw BoutException(
        "Error: Div = Curl = 0 boundary condition doesn't work for MXG > 2. Sorry\n");
  }

  jx = mesh->xend + 1;
  for (jy = 1; jy < mesh->LocalNy - 1; jy++) {
    for (jz = 0; jz < ncz; jz++) {
      jzp = (jz + 1) % ncz;
      jzm = (jz - 1 + ncz) % ncz;

      // dB_y / dx = dB_x / dy

      // dB_x / dy
      tmp = (var.x(jx - 1, jy + 1, jz) - var.x(jx - 1, jy - 1, jz))
            / (metric->dy(jx - 1, jy - 1) + metric->dy(jx - 1, jy));

      var.y(jx, jy, jz) =
          var.y(jx - 2, jy, jz) + (metric->dx(jx - 2, jy) + metric->dx(jx - 1, jy)) * tmp;
      if (mesh->xstart == 2)
        // 4th order to get last point
        var.y(jx + 1, jy, jz) = var.y(jx - 3, jy, jz) + 4. * metric->dx(jx, jy) * tmp;

      // dB_z / dx = dB_x / dz

      tmp = (var.x(jx - 1, jy, jzp) - var.x(jx - 1, jy, jzm)) / (2. * metric->dz(jx - 1, jy));

      var.z(jx, jy, jz) =
          var.z(jx - 2, jy, jz) + (metric->dx(jx - 2, jy) + metric->dx(jx - 1, jy)) * tmp;
      if (mesh->xstart == 2)
        var.z(jx + 1, jy, jz) = var.z(jx - 3, jy, jz) + 4. * metric->dx(jx, jy) * tmp;

      // d/dx( Jmetric->g11 B_x ) = - d/dx( Jmetric->g12 B_y + Jmetric->g13 B_z)
      //                    - d/dy( JB^y ) - d/dz( JB^z )

      tmp = -(metric->J(jx, jy) * metric->g12(jx, jy) * var.y(jx, jy, jz)
              + metric->J(jx, jy) * metric->g13(jx, jy) * var.z(jx, jy, jz)
              - metric->J(jx - 2, jy) * metric->g12(jx - 2, jy) * var.y(jx - 2, jy, jz)
              + metric->J(jx - 2, jy) * metric->g13(jx - 2, jy) * var.z(jx - 2, jy, jz))
            / (metric->dx(jx - 2, jy)
               + metric->dx(jx - 1, jy)); // First term (d/dx) using vals calculated above
      tmp -= (metric->J(jx - 1, jy + 1) * metric->g12(jx - 1, jy + 1)
                  * var.x(jx - 1, jy + 1, jz)
              - metric->J(jx - 1, jy - 1) * metric->g12(jx - 1, jy - 1)
                    * var.x(jx - 1, jy - 1, jz)
              + metric->J(jx - 1, jy + 1) * metric->g22(jx - 1, jy + 1)
                    * var.y(jx - 1, jy + 1, jz)
              - metric->J(jx - 1, jy - 1) * metric->g22(jx - 1, jy - 1)
                    * var.y(jx - 1, jy - 1, jz)
              + metric->J(jx - 1, jy + 1) * metric->g23(jx - 1, jy + 1)
                    * var.z(jx - 1, jy + 1, jz)
              - metric->J(jx - 1, jy - 1) * metric->g23(jx - 1, jy - 1)
                    * var.z(jx - 1, jy - 1, jz))
             / (metric->dy(jx - 1, jy - 1) + metric->dy(jx - 1, jy)); // second (d/dy)
      tmp -= (metric->J(jx - 1, jy) * metric->g13(jx - 1, jy)
                  * (var.x(jx - 1, jy, jzp) - var.x(jx - 1, jy, jzm))
              + metric->J(jx - 1, jy) * metric->g23(jx - 1, jy)
                    * (var.y(jx - 1, jy, jzp) - var.y(jx - 1, jy, jzm))
              + metric->J(jx - 1, jy) * metric->g33(jx - 1, jy)
                    * (var.z(jx - 1, jy, jzp) - var.z(jx - 1, jy, jzm)))
             / (2. * metric->dz(jx - 1, jy));

      var.x(jx, jy, jz) =
          (metric->J(jx - 2, jy) * metric->g11(jx - 2, jy) * var.x(jx - 2, jy, jz)
           + (metric->dx(jx - 2, jy) + metric->dx(jx - 1, jy)) * tmp)
          / metric->J(jx, jy) * metric->g11(jx, jy);
      if (mesh->xstart == 2)
        var.x(jx + 1, jy, jz) =
            (metric->J(jx - 3, jy) * metric->g11(jx - 3, jy) * var.x(jx - 3, jy, jz)
             + 4. * metric->dx(jx, jy) * tmp)
            / metric->J(jx + 1, jy) * metric->g11(jx + 1, jy);
    }
  }
#else
  throw BoutException(
      "Applying boundary to Vector3D not compatible with 3D metrics in DivCurl.");
#endif
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryFree::clone(BoundaryRegion* region,
                                const std::list<std::string>& args) {
  if (!args.empty()) {
    // First argument should be a value
    val = stringToReal(args.front());
    return new BoundaryFree(region, val);
  }
  return new BoundaryFree(region);
}

void BoundaryFree::apply(Field2D& UNUSED(f)) {
  // Do nothing for free boundary
}

void BoundaryFree::apply(Field3D& UNUSED(f)) {
  // Do nothing for free boundary
}

void BoundaryFree::apply_ddt(Field2D& UNUSED(f)) {
  // Do nothing for free boundary
}

void BoundaryFree::apply_ddt(Field3D& UNUSED(f)) {
  // Do nothing for free boundary
}
///////////////////////////////////////////////////////////////
// New free boundary implementation.  Uses last grid points to extrapolate into the guard
// cells. Written by L. Easy.
///////////////////////////////////////////////////////////////

// 2nd order extrapolation:

BoundaryOp* BoundaryFree_O2::clone(BoundaryRegion* region,
                                   const std::list<std::string>& args) {
  verifyNumPoints(region, 2);
  if (!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryFree\n";
  }
  return new BoundaryFree_O2(region);
}

void BoundaryFree_O2::apply(Field2D& f) {
  // Set (at 2nd order) the value at the mid-point between the guard cell and the grid
  // cell to be val N.B. Only first guard cells (closest to the grid) should ever be used

  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  bndry->first();

  // Check for staggered grids

  CELL_LOC loc = f.getLocation();
  if (mesh->StaggerGrids and (loc == CELL_XLOW or loc == CELL_YLOW)) {
    // Staggered. Need to apply slightly differently

    if (loc == CELL_XLOW) {
      // Field is shifted in X

      if (bndry->bx > 0) {
        // Outer x boundary
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int i = 0; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 2 * f(xi - bndry->bx, yi - bndry->by)
                        - f(xi - 2 * bndry->bx, yi - 2 * bndry->by);
          }
        }
      }
      if (bndry->bx < 0) {
        // Inner x boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int i = -1; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 2 * f(xi - bndry->bx, yi - bndry->by)
                        - f(xi - 2 * bndry->bx, yi - 2 * bndry->by);
          }
        }
      }
      if (bndry->by != 0) {
        // y boundaries
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int i = 0; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 2 * f(xi - bndry->bx, yi - bndry->by)
                        - f(xi - 2 * bndry->bx, yi - 2 * bndry->by);
          }
        }
      }
    } else if (loc == CELL_YLOW) {
      // Field is shifted in Y

      if (bndry->by > 0) {
        // Upper y boundary

        for (; !bndry->isDone(); bndry->next1d()) {
          for (int i = 0; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 2 * f(xi - bndry->bx, yi - bndry->by)
                        - f(xi - 2 * bndry->bx, yi - 2 * bndry->by);
          }
        }
      }
      if (bndry->by < 0) {
        // Lower y boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int i = -1; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 2 * f(xi - bndry->bx, yi - bndry->by)
                        - f(xi - 2 * bndry->bx, yi - 2 * bndry->by);
          }
        }
      }
      if (bndry->bx != 0) {
        // x boundaries
        for (; !bndry->isDone(); bndry->next1d()) {

          for (int i = 0; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 2 * f(xi - bndry->bx, yi - bndry->by)
                        - f(xi - 2 * bndry->bx, yi - 2 * bndry->by);
          }
        }
      }
    }
  } else {
    // Non-staggered, standard case

    for (; !bndry->isDone(); bndry->next1d()) {

      for (int i = 0; i < bndry->width; i++) {
        int xi = bndry->x + i * bndry->bx;
        int yi = bndry->y + i * bndry->by;
        f(xi, yi) = 2 * f(xi - bndry->bx, yi - bndry->by)
                    - f(xi - 2 * bndry->bx, yi - 2 * bndry->by);
      }
    }
  }
}

void BoundaryFree_O2::apply(Field3D& f) {
  // Extrapolate from the last evolved simulation cells into the guard cells at 3rd order.

  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  bndry->first();

  // Check for staggered grids

  CELL_LOC loc = f.getLocation();
  if (mesh->StaggerGrids and (loc == CELL_XLOW or loc == CELL_YLOW)) {
    // Staggered. Need to apply slightly differently

    if (loc == CELL_XLOW) {
      // Field is shifted in X

      if (bndry->bx > 0) {
        // Outer x boundary

        for (; !bndry->isDone(); bndry->next1d()) {

          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            for (int i = 0; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 2 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->bx < 0) {
        // Inner x boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {

          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            for (int i = -1; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 2 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->by != 0) {
        // y boundaries

        for (; !bndry->isDone(); bndry->next1d()) {

          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            for (int i = 0; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 2 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk);
            }
          }
        }
      }
    } else if (loc == CELL_YLOW) {
      // Field is shifted in Y

      if (bndry->by > 0) {
        // Upper y boundary
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            for (int i = 0; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 2 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->by < 0) {
        // Lower y boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {

          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            for (int i = -1; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 2 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->bx != 0) {
        // x boundaries
        for (; !bndry->isDone(); bndry->next1d()) {

          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            for (int i = 0; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 2 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk);
            }
          }
        }
      }
    }
  } else {
    // Standard (non-staggered) case
    for (; !bndry->isDone(); bndry->next1d()) {

      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        for (int i = 0; i < bndry->width; i++) {
          int xi = bndry->x + i * bndry->bx;
          int yi = bndry->y + i * bndry->by;
          f(xi, yi, zk) = 2 * f(xi - bndry->bx, yi - bndry->by, zk)
                          - f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk);
        }
      }
    }
  }
}

void BoundaryFree_O2::apply_ddt(Field2D& f) {
  Field2D* dt = f.timeDeriv();
  for (bndry->first(); !bndry->isDone(); bndry->next())
    (*dt)(bndry->x, bndry->y) = 0.; // Set time derivative to zero
}

void BoundaryFree_O2::apply_ddt(Field3D& f) {
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  Field3D *dt = f.timeDeriv();
  for (bndry->first(); !bndry->isDone(); bndry->next())
    for (int z = 0; z < mesh->LocalNz; z++)
      (*dt)(bndry->x,bndry->y,z) = 0.; // Set time derivative to zero
}

//////////////////////////////////
// Third order extrapolation:
//////////////////////////////////
BoundaryOp* BoundaryFree_O3::clone(BoundaryRegion* region,
                                   const std::list<std::string>& args) {
  verifyNumPoints(region, 3);

  if (!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryConstLaplace\n";
  }
  return new BoundaryFree_O3(region);
}

void BoundaryFree_O3::apply(Field2D& f) {

  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  bndry->first();

  // Check for staggered grids

  CELL_LOC loc = f.getLocation();
  if (mesh->StaggerGrids and (loc == CELL_XLOW or loc == CELL_YLOW)) {
    // Staggered. Need to apply slightly differently

    if (loc == CELL_XLOW) {
      // Field is shifted in X

      if (bndry->bx > 0) {
        // Outer x boundary
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int i = 0; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 3.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + f(xi - 3 * bndry->bx, yi - 3 * bndry->by);
          }
        }
      }
      if (bndry->bx < 0) {
        // Inner x boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int i = -1; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 3.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + f(xi - 3 * bndry->bx, yi - 3 * bndry->by);
          }
        }
      }
      if (bndry->by != 0) {
        // y boundaries
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int i = 0; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 3.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + f(xi - 3 * bndry->bx, yi - 3 * bndry->by);
          }
        }
      }
    } else if (loc == CELL_YLOW) {
      // Field is shifted in Y

      if (bndry->by > 0) {
        // Upper y boundary

        for (; !bndry->isDone(); bndry->next1d()) {
          for (int i = 0; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 3.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + f(xi - 3 * bndry->bx, yi - 3 * bndry->by);
          }
        }
      }
      if (bndry->by < 0) {
        // Lower y boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int i = -1; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 3.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + f(xi - 3 * bndry->bx, yi - 3 * bndry->by);
          }
        }
      }
      if (bndry->bx != 0) {
        // x boundaries
        for (; !bndry->isDone(); bndry->next1d()) {

          for (int i = 0; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            f(xi, yi) = 3.0 * f(xi - bndry->bx, yi - bndry->by)
                        - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                        + f(xi - 3 * bndry->bx, yi - 3 * bndry->by);
          }
        }
      }
    }
  } else {
    // Non-staggered, standard case

    for (; !bndry->isDone(); bndry->next1d()) {

      for (int i = 0; i < bndry->width; i++) {
        int xi = bndry->x + i * bndry->bx;
        int yi = bndry->y + i * bndry->by;
        f(xi, yi) = 3.0 * f(xi - bndry->bx, yi - bndry->by)
                    - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by)
                    + f(xi - 3 * bndry->bx, yi - 3 * bndry->by);
      }
    }
  }
}

void BoundaryFree_O3::apply(Field3D& f) {
  // Extrapolate from the last evolved simulation cells into the guard cells at 3rd order.

  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  bndry->first();

  // Check for staggered grids

  CELL_LOC loc = f.getLocation();
  if (mesh->StaggerGrids and (loc == CELL_XLOW or loc == CELL_YLOW)) {
    // Staggered. Need to apply slightly differently

    if (loc == CELL_XLOW) {
      // Field is shifted in X

      if (bndry->bx > 0) {
        // Outer x boundary

        for (; !bndry->isDone(); bndry->next1d()) {

          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            for (int i = 0; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 3.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->bx < 0) {
        // Inner x boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {

          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            for (int i = -1; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 3.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->by != 0) {
        // y boundaries

        for (; !bndry->isDone(); bndry->next1d()) {

          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            for (int i = 0; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 3.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk);
            }
          }
        }
      }
    } else if (loc == CELL_YLOW) {
      // Field is shifted in Y

      if (bndry->by > 0) {
        // Upper y boundary
        for (; !bndry->isDone(); bndry->next1d()) {
          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            for (int i = 0; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 3.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->by < 0) {
        // Lower y boundary. Set one point inwards
        for (; !bndry->isDone(); bndry->next1d()) {

          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            for (int i = -1; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 3.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk);
            }
          }
        }
      }
      if (bndry->bx != 0) {
        // x boundaries
        for (; !bndry->isDone(); bndry->next1d()) {

          for (int zk = 0; zk < mesh->LocalNz; zk++) {
            for (int i = 0; i < bndry->width; i++) {
              int xi = bndry->x + i * bndry->bx;
              int yi = bndry->y + i * bndry->by;
              f(xi, yi, zk) = 3.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                              - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                              + f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk);
            }
          }
        }
      }
    }
  } else {
    // Standard (non-staggered) case
    for (; !bndry->isDone(); bndry->next1d()) {

      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        for (int i = 0; i < bndry->width; i++) {
          int xi = bndry->x + i * bndry->bx;
          int yi = bndry->y + i * bndry->by;
          f(xi, yi, zk) = 3.0 * f(xi - bndry->bx, yi - bndry->by, zk)
                          - 3.0 * f(xi - 2 * bndry->bx, yi - 2 * bndry->by, zk)
                          + f(xi - 3 * bndry->bx, yi - 3 * bndry->by, zk);
        }
      }
    }
  }
}

void BoundaryFree_O3::apply_ddt(Field2D& f) {
  Field2D* dt = f.timeDeriv();
  for (bndry->first(); !bndry->isDone(); bndry->next())
    (*dt)(bndry->x, bndry->y) = 0.; // Set time derivative to zero
}

void BoundaryFree_O3::apply_ddt(Field3D& f) {
  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());
  Field3D *dt = f.timeDeriv();
  for (bndry->first(); !bndry->isDone(); bndry->next())
    for (int z = 0; z < mesh->LocalNz; z++)
      (*dt)(bndry->x,bndry->y,z) = 0.; // Set time derivative to zero
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryRelax::cloneMod(BoundaryOp* operation,
                                    const std::list<std::string>& args) {
  auto* result = new BoundaryRelax(operation, r);

  if (!args.empty()) {
    // First argument should be the rate
    BoutReal val = stringToReal(args.front());
    val = fabs(val); // Should always be positive
    result->r = val;
  }

  return result;
}

void BoundaryRelax::apply(Field2D& f, BoutReal t) {
  // Just apply the original boundary condition to f
  op->apply(f, t);
}

void BoundaryRelax::apply(Field3D& f, BoutReal UNUSED(t)) {
  // Just apply the original boundary condition to f
  op->apply(f);
}

void BoundaryRelax::apply_ddt(Field2D& f) {
  TRACE("BoundaryRelax::apply_ddt(Field2D)");

  // Make a copy of f
  Field2D g = f;
  // Apply the boundary to g
  op->apply(g);

  bndry->first();

  // Set time-derivatives
  for (bndry->first(); !bndry->isDone(); bndry->next()) {
    ddt(f)(bndry->x, bndry->y) = r * (g(bndry->x, bndry->y) - f(bndry->x, bndry->y));
  }
}

void BoundaryRelax::apply_ddt(Field3D& f) {
  TRACE("BoundaryRelax::apply_ddt(Field3D)");

  Mesh* mesh = bndry->localmesh;
  ASSERT1(mesh == f.getMesh());

  // Make a copy of f
  Field3D g = f; // NOTE: This is not very efficient... copying entire field
  // Apply the boundary to g
  op->apply(g);
  // Set time-derivatives
  for (bndry->first(); !bndry->isDone(); bndry->next())
    for (int z = 0; z < mesh->LocalNz; z++) {
      ddt(f)(bndry->x, bndry->y, z) =
          r * (g(bndry->x, bndry->y, z) - f(bndry->x, bndry->y, z));
    }
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryWidth::cloneMod(BoundaryOp* operation,
                                    const std::list<std::string>& args) {
  auto* result = new BoundaryWidth(operation, width);

  if (args.empty()) {
    output << "WARNING: BoundaryWidth expected 1 argument\n";
  } else {
    // First argument should be the rate
    int val = stringToInt(args.front());
    result->width = val;
  }

  return result;
}

void BoundaryWidth::apply(Field2D& f, BoutReal t) {
  // Pointer to boundary region shared between all BoundaryOp, BoundaryModifiers
  int oldwid = bndry->width;
  bndry->width = width;
  op->apply(f, t);
  bndry->width = oldwid;
}

void BoundaryWidth::apply(Field3D& f, BoutReal t) {
  int oldwid = bndry->width;
  bndry->width = width;
  op->apply(f, t);
  bndry->width = oldwid;
}

void BoundaryWidth::apply_ddt(Field2D& f) {
  int oldwid = bndry->width;
  bndry->width = width;
  op->apply_ddt(f);
  bndry->width = oldwid;
}

void BoundaryWidth::apply_ddt(Field3D& f) {
  int oldwid = bndry->width;
  bndry->width = width;
  op->apply_ddt(f);
  bndry->width = oldwid;
}

///////////////////////////////////////////////////////////////
BoundaryOp* BoundaryToFieldAligned::cloneMod(BoundaryOp* operation,
                                             const std::list<std::string>& args) {
  auto* result = new BoundaryToFieldAligned(operation);

  if (!args.empty()) {
    output << "WARNING: BoundaryToFieldAligned expected no argument\n";
    // Shouldn't we throw ?
  }

  return result;
}

void BoundaryToFieldAligned::apply(Field2D& f, BoutReal t) { op->apply(f, t); }

void BoundaryToFieldAligned::apply(Field3D &f, BoutReal t) {
  ASSERT1(bndry->localmesh == f.getMesh());

  // NOTE: This is not very efficient... updating entire field
  f = fromFieldAligned(f);

  // Apply the boundary to shifted field
  op->apply(f, t);

  // Shift back
  f = toFieldAligned(f);

  // This is inefficient -- could instead use the shiftZ just in the bndry
  // but this is not portable to other parallel transforms -- we could instead
  // have a flag to define the region in which we want to apply to/fromFieldAligned
}

void BoundaryToFieldAligned::apply_ddt(Field2D &f) {
  op->apply_ddt(f);
}

void BoundaryToFieldAligned::apply_ddt(Field3D &f) {
  ASSERT1(bndry->localmesh == f.getMesh());

  f = fromFieldAligned(f);
  ddt(f) = fromFieldAligned(ddt(f));
  op->apply_ddt(f);
  ddt(f) = toFieldAligned(ddt(f));
}

///////////////////////////////////////////////////////////////
BoundaryOp* BoundaryFromFieldAligned::cloneMod(BoundaryOp* operation,
                                               const std::list<std::string>& args) {
  auto* result = new BoundaryFromFieldAligned(operation);

  if (!args.empty()) {
    output << "WARNING: BoundaryFromFieldAligned expected no argument\n";
    // Shouldn't we throw ?
  }

  return result;
}

void BoundaryFromFieldAligned::apply(Field2D& f, BoutReal t) { op->apply(f, t); }

void BoundaryFromFieldAligned::apply(Field3D &f, BoutReal t) {
  ASSERT1(bndry->localmesh == f.getMesh());

  // NOTE: This is not very efficient... shifting entire field
  f = toFieldAligned(f);

  // Apply the boundary to shifted field
  op->apply(f, t);

  // Shift back
  f = fromFieldAligned(f);

  // This is inefficient -- could instead use the shiftZ just in the bndry
  // but this is not portable to other parallel transforms -- we could instead
  // have a flag to define the region in which we want to apply to/fromFieldAligned
}

void BoundaryFromFieldAligned::apply_ddt(Field2D &f) {
  op->apply_ddt(f);
}

void BoundaryFromFieldAligned::apply_ddt(Field3D &f) {
  ASSERT1(bndry->localmesh == f.getMesh());

  f = toFieldAligned(f);
  ddt(f) = toFieldAligned(ddt(f));
  op->apply_ddt(f);
  ddt(f) = fromFieldAligned(ddt(f));
}
