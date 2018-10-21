#include <globals.hxx>
#include <boundary_standard.hxx>
#include <invert_laplace.hxx>
#include <fft.hxx>
#include <boutexception.hxx>
#include <utils.hxx>
#include <output.hxx>
#include <msg_stack.hxx>
#include <bout/constants.hxx>
#include <derivs.hxx>

///////////////////////////////////////////////////////////////
// Helpers

/** \brief Check that there are sufficient non-boundary points for desired B.C.
    
    Checks both the size of the global grid (i.e. if this B.C. could be ok
    for some parallel setup or not) and the local grid.

    Note the local grid check is not strictly necessary as this would typically
    lead to an out of bounds access error later but we add it here to provide a
    more explanatory message.
 */
namespace {
  void verifyNumPoints(BoundaryRegion *region, int ptsRequired) {
    TRACE("Verifying number of points available for BC");

#ifndef CHECK
    return; //No checking so just return
#else

    Mesh* localmesh = region->localmesh;

    int ptsAvailGlobal, ptsAvailLocal, ptsAvail;
    string side, gridType;

    //Initialise var in case of no match and CHECK<=2
    ptsAvail = ptsRequired; //Ensures test passes without exception

    switch(region->location) {
    case BNDRY_XIN:
    case BNDRY_XOUT: {
      side = "x";

      //Here 2*localmesh->xstart is the total number of guard/boundary cells
      ptsAvailGlobal = localmesh->GlobalNx - 2*localmesh->xstart;

      //Work out how many processor local points we have excluding boundaries
      //but including ghost/guard cells
      ptsAvailLocal  = localmesh->LocalNx;
      if(localmesh->firstX()) ptsAvailLocal -= localmesh->xstart;
      if(localmesh->lastX())  ptsAvailLocal -= localmesh->xstart;

      //Now decide if it's a local or global limit, prefer global if a tie
      if(ptsAvailGlobal <= ptsAvailLocal){
        ptsAvail = ptsAvailGlobal;
        gridType = "global";
      }else{
        ptsAvail = ptsAvailLocal;
        gridType = "local";
      }

      break;
    }
    case BNDRY_YUP:
    case BNDRY_YDOWN: {
      side = "y";

      //Here 2*localmesh->ystart is the total number of guard/boundary cells
      ptsAvailGlobal = localmesh->GlobalNy - 2*localmesh->ystart;

      //Work out how many processor local points we have excluding boundaries
      //but including ghost/guard cells
      ptsAvailLocal  = localmesh->LocalNy;
      if(localmesh->firstY()) ptsAvailLocal -= localmesh->ystart;
      if(localmesh->lastY())  ptsAvailLocal -= localmesh->ystart;

      //Now decide if it's a local or global limit, prefer global if a tie
      if(ptsAvailGlobal <= ptsAvailLocal){
        ptsAvail = ptsAvailGlobal;
        gridType = "global";
      }else{
        ptsAvail = ptsAvailLocal;
        gridType = "local";
      }

      break;
    }
#if CHECK > 2 //Only fail on Unrecognised boundary for extreme checking
    default : {
      throw BoutException("Unrecognised boundary region (%s) for verifyNumPoints.",region->location);
    }
#endif
    }

    //Now check we have enough points and if not throw an exception
    if(ptsAvail < ptsRequired){
      throw BoutException("Too few %s grid points for %s boundary, have %d but need at least %d",
                          gridType.c_str(),side.c_str(),ptsAvail,ptsRequired);
    }

#endif
  }

  // 2nd order extrapolation to a point
  template<typename T>
  void extrap2nd(T &f, int x, int bx, int y, int by, int z) {
    f(x, y, z) = 2*f(x - bx, y - by, z) - f(x - 2*bx, y - 2*by, z);
  }

  // 3rd order extrapolation to a point
  template<typename T>
  void extrap3rd(T &f, int x, int bx, int y, int by, int z) {
    f(x, y, z) = 3.0*f(x - bx, y - by, z) - 3.0*f(x - 2*bx, y - 2*by, z) + f(x - 3*bx, y - 3*by, z);
  }

  // 4th order extrapolation to a point
  template<typename T>
  void extrap4th(T &f, int x, int bx, int y, int by, int z) {
    f(x, y, z) = 4.0*f(x - bx, y - by, z) - 6.0*f(x - 2*bx, y - 2*by, z)
      + 4.0*f(x - 3*bx, y - 3*by, z) - f(x - 4*bx, y - 4*by, z);
  }

  // 5th order extrapolation to a point
  template<typename T>
  void extrap5th(T &f, int x, int bx, int y, int by, int z) {
    f(x, y, z) = 5.0*f(x - bx, y - by, z) - 10.0*f(x - 2*bx, y - 2*by, z)
      + 10.0*f(x - 3*bx, y - 3*by, z) - 5.0*f(x - 4*bx, y - 4*by, z)
      + f(x - 5*bx, y - 5*by, z);
  }
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryDirichlet::clone(BoundaryRegion *region, const list<string> &args){
  verifyNumPoints(region, 1);

  std::shared_ptr<FieldGenerator> newgen = nullptr;
  if(!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryDirichlet(region, newgen);
}

// Override apply(), using this private method to provide both Field2D and
// Field3D versions, for BoundaryDirichlet because we apply a funny hack to the
// extra guard cells
template<typename T>
void BoundaryDirichlet::applyTemplate(T &f,BoutReal t) {
  // Set (at 2nd order) the value at the mid-point between the guard cell and the grid cell to be val
  // N.B. Only first guard cells (closest to the grid) should ever be used

  Mesh* localmesh = f.getMesh();

  // Check for staggered grids
  CELL_LOC loc = f.getLocation();
  ASSERT1(localmesh->StaggerGrids || loc == CELL_CENTRE);

  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator>  fg = gen;
  if(!fg) {
    fg = f.getBndryGenerator(bndry->location);
  }

  BoutReal val = 0.0;

  if (loc == CELL_CENTRE) {
    // Unstaggered case
    for(; !bndry->isDone(); bndry->next1d()) {
      // Calculate the X and Y normalised values half-way between the guard cell and grid cell 
      BoutReal xnorm = 0.5*(   localmesh->GlobalX(bndry->x)  // In the guard cell
          + localmesh->GlobalX(bndry->x - bndry->bx) ); // the grid cell

      BoutReal ynorm = 0.5*(   localmesh->GlobalY(bndry->y)  // In the guard cell
          + localmesh->GlobalY(bndry->y - bndry->by) ); // the grid cell

      for(int zk=0; zk<f.getNz(); zk++) {
        BoutReal znorm = localmesh->GlobalZ(zk);
        if(fg){
          val = fg->generate(xnorm,TWOPI*ynorm,TWOPI*znorm, t);
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
        xnorm = localmesh->GlobalX(xi);
        ynorm = localmesh->GlobalY(yi);
        for(int zk=0; zk<f.getNz() ; zk++) {
          BoutReal znorm = localmesh->GlobalZ(zk);
          if(fg) {
            val = fg->generate(xnorm,TWOPI*ynorm,TWOPI*znorm, t);
          }
          f(xi, yi, zk) = val;
        }
      }
    }
  } else if( loc == CELL_XLOW ) {
    // Field is shifted in X
    if(bndry->bx > 0) {
      // Outer x boundary
      for(; !bndry->isDone(); bndry->next1d()) {
        BoutReal xnorm = 0.5*(   localmesh->GlobalX(bndry->x)
            + localmesh->GlobalX(bndry->x - bndry->bx) );
        BoutReal ynorm = localmesh->GlobalY(bndry->y);

        for(int zk=0; zk<f.getNz(); zk++) {
          BoutReal znorm = localmesh->GlobalZ(zk);
          if(fg){
            val = fg->generate(xnorm,TWOPI*ynorm,TWOPI*znorm, t);
          }
          f(bndry->x,bndry->y, zk) = val;

          // Need to set second guard cell, as may be used for interpolation or upwinding derivatives
          for(int i=1;i<bndry->width;i++) {
            int xi = bndry->x + i*bndry->bx;
            int yi = bndry->y ;

            f(xi, yi, zk) = 2*f(xi - bndry->bx, yi , zk) - f(xi- 2*bndry->bx, yi , zk);
          }
        }
      }
    }
    if (bndry->bx < 0){
      // Inner x boundary. Set one point inwards
      for(; !bndry->isDone(); bndry->next1d()) {

        BoutReal xnorm = 0.5*(   localmesh->GlobalX(bndry->x)
            + localmesh->GlobalX(bndry->x - bndry->bx) );
        BoutReal ynorm = localmesh->GlobalY(bndry->y);

        for(int zk=0; zk<f.getNz(); zk++) {
          BoutReal znorm = localmesh->GlobalZ(zk);
          if(fg){
            val = fg->generate(xnorm,TWOPI*ynorm,TWOPI*znorm, t);
          }
          f(bndry->x - bndry->bx,bndry->y, zk) = val;

          // Need to set second guard cell, as may be used for interpolation or upwinding derivatives
          for(int i=0;i<bndry->width;i++) {
            int xi = bndry->x + i*bndry->bx;
            int yi = bndry->y ;

            f(xi, yi, zk) = 2*f(xi - bndry->bx, yi , zk) - f(xi- 2*bndry->bx, yi , zk);
          }
        }
      }
    }
    if(bndry->by !=0){
      // y boundaries
      for(; !bndry->isDone(); bndry->next1d()) {
        // x norm is shifted by half a grid point because it is staggered.
        // y norm is located half way between first grid cell and guard cell.
        BoutReal xnorm = 0.5*(   localmesh->GlobalX(bndry->x) + localmesh->GlobalX(bndry->x - 1) );
        BoutReal ynorm = 0.5*(   localmesh->GlobalY(bndry->y) + localmesh->GlobalY(bndry->y - bndry->by) );

        for(int zk=0; zk<f.getNz(); zk++) {
          BoutReal znorm = localmesh->GlobalZ(zk);
          if(fg){
            val = fg->generate(xnorm,TWOPI*ynorm,TWOPI*znorm, t);
          }
          f(bndry->x,bndry->y,zk) = 2*val - f(bndry->x-bndry->bx, bndry->y-bndry->by, zk);

          // Need to set second guard cell, as may be used for interpolation or upwinding derivatives
          for(int i=1;i<bndry->width;i++) {
            int xi = bndry->x ;
            int yi = bndry->y + i*bndry->by;

            f(xi, yi, zk) = 2*f(xi, yi - bndry->by, zk) - f(xi, yi - 2*bndry->by, zk);
          }
        }
      }
    }
  } else if( loc == CELL_YLOW ) {
    // Shifted in Y
    if(bndry->by > 0) {
      // Upper y boundary boundary
      for(; !bndry->isDone(); bndry->next1d()) {
        BoutReal xnorm = localmesh->GlobalX(bndry->x);
        BoutReal ynorm = 0.5*(   localmesh->GlobalY(bndry->y) + localmesh->GlobalY(bndry->y - bndry->by) );
        for(int zk=0; zk<f.getNz(); zk++) {
          BoutReal znorm = localmesh->GlobalZ(zk);
          if(fg){
            val = fg->generate(xnorm,TWOPI*ynorm,TWOPI*znorm, t);
          }
          f(bndry->x,bndry->y,zk) = val;

          // Need to set second guard cell, as may be used for interpolation or upwinding derivatives
          for(int i=1;i<bndry->width;i++) {
            int xi = bndry->x ;
            int yi = bndry->y + i*bndry->by;

            f(xi, yi, zk) = 2.0*f(xi, yi - bndry->by, zk) - f(xi, yi - 2*bndry->by, zk);
          }
        }
      }
    }
    if(bndry->by < 0){
      // Lower y boundary. Set one point inwards
      for(; !bndry->isDone(); bndry->next1d()) {

        BoutReal xnorm = localmesh->GlobalX(bndry->x);
        BoutReal ynorm = 0.5*(   localmesh->GlobalY(bndry->y) + localmesh->GlobalY(bndry->y - bndry->by) );

        for(int zk=0; zk<f.getNz(); zk++) {
          BoutReal znorm = localmesh->GlobalZ(zk);
          if(fg){
            val = fg->generate(xnorm,TWOPI*ynorm,TWOPI*znorm, t);
          }
          f(bndry->x,bndry->y - bndry->by, zk) = val;

          // Need to set second guard cell, as may be used for interpolation or upwinding derivatives
          for(int i=0;i<bndry->width;i++) {
            int xi = bndry->x ;
            int yi = bndry->y + i*bndry->by;

            f(xi, yi, zk) = 2*f(xi, yi - bndry->by, zk) - f(xi, yi - 2*bndry->by, zk);
          }
        }
      }
    }
    if(bndry->bx != 0){
      // x boundaries
      for(; !bndry->isDone(); bndry->next1d()) {
        // x norm is located half way between first grid cell and guard cell.
        // y norm is shifted by half a grid point because it is staggered.
        BoutReal xnorm = 0.5*( localmesh->GlobalX(bndry->x) + localmesh->GlobalX(bndry->x - bndry->bx) );
        BoutReal ynorm = 0.5*( localmesh->GlobalY(bndry->y) + localmesh->GlobalY(bndry->y - 1) );

        for(int zk=0; zk<f.getNz(); zk++) {
          BoutReal znorm = localmesh->GlobalZ(zk);
          if(fg) {
            val = fg->generate(xnorm,TWOPI*ynorm,TWOPI*znorm, t);
          }

          f(bndry->x,bndry->y,zk) = 2*val - f(bndry->x-bndry->bx, bndry->y-bndry->by, zk);

          // Need to set second guard cell, as may be used for interpolation or upwinding derivatives
          for(int i=1;i<bndry->width;i++) {
            int xi = bndry->x + i*bndry->bx;
            int yi = bndry->y ;

            f(xi, yi, zk) = 2*f(xi - bndry->bx, yi , zk) - f(xi - 2*bndry->bx, yi, zk);
          }
        }
      }
    }
  } else if (loc == CELL_ZLOW) {
    // Staggered in Z. Note there are no z-boundaries.
    for(; !bndry->isDone(); bndry->next1d()) {
      // Calculate the X and Y normalised values half-way between the guard cell and grid cell 
      BoutReal xnorm = 0.5*(   localmesh->GlobalX(bndry->x)  // In the guard cell
          + localmesh->GlobalX(bndry->x - bndry->bx) ); // the grid cell

      BoutReal ynorm = 0.5*(   localmesh->GlobalY(bndry->y)  // In the guard cell
          + localmesh->GlobalY(bndry->y - bndry->by) ); // the grid cell

      for(int zk=0; zk<f.getNz(); zk++) {
        // It shouldn't matter if znorm<0 because the expression in fg->generate should be periodic in z
        BoutReal znorm = 0.5*( localmesh->GlobalZ(zk) + localmesh->GlobalZ(zk - 1) ); // znorm is shifted by half a grid point because it is staggered
        if(fg){
          val = fg->generate(xnorm,TWOPI*ynorm,TWOPI*znorm, t);
        }
        f(bndry->x,bndry->y,zk) = 2*val - f(bndry->x-bndry->bx, bndry->y-bndry->by, zk);

        // Need to set second guard cell, as may be used for interpolation or upwinding derivatives
        for(int i=1;i<bndry->width;i++) {
          int xi = bndry->x + i*bndry->bx;
          int yi = bndry->y ;

          f(xi, yi, zk) = 2*f(xi - bndry->bx, yi , zk) - f(xi - 2*bndry->bx, yi, zk);
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////
// New implementation, accurate to higher order

BoundaryOp* BoundaryDirichlet_O3::clone(BoundaryRegion *region, const list<string> &args){
  verifyNumPoints(region, 2);

  std::shared_ptr<FieldGenerator> newgen = nullptr;
  if(!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryDirichlet_O3(region, newgen);
}

void BoundaryDirichlet_O3::applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = (8./3)*val - 2.*f(x - bx, y - by, z) + f(x - 2*bx, y - 2*by, z)/3.;
}
void BoundaryDirichlet_O3::applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = (8./3)*val - 2.*f(x - bx, y - by, z) + f(x - 2*bx, y - 2*by, z)/3.;
}

void BoundaryDirichlet_O3::applyAtPointStaggered(Field2D &f, BoutReal val, int x, int UNUSED(bx), int y, int UNUSED(by), int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = val;
}
void BoundaryDirichlet_O3::applyAtPointStaggered(Field3D &f, BoutReal val, int x, int UNUSED(bx), int y, int UNUSED(by), int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = val;
}

void BoundaryDirichlet_O3::extrapFurther(Field2D &f, int x, int bx, int y, int by, int z) {
  extrap3rd(f, x, bx, y, by, z);
}
void BoundaryDirichlet_O3::extrapFurther(Field3D &f, int x, int bx, int y, int by, int z) {
  extrap3rd(f, x, bx, y, by, z);
}

///////////////////////////////////////////////////////////////
// Extrapolate to calculate boundary cell to 4th-order

BoundaryOp* BoundaryDirichlet_O4::clone(BoundaryRegion *region, const list<string> &args){
  verifyNumPoints(region, 3);

  std::shared_ptr<FieldGenerator> newgen = nullptr;
  if(!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryDirichlet_O4(region, newgen);
}

void BoundaryDirichlet_O4::applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = (16./5)*val - 3.*f(x - bx, y - by, z) + f(x - 2*bx, y - 2*by, z) - (1./5)*f(x - 3*bx, y - 3*by, z);
}
void BoundaryDirichlet_O4::applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = (16./5)*val - 3.*f(x - bx, y - by, z) + f(x - 2*bx, y - 2*by, z) - (1./5)*f(x - 3*bx, y - 3*by, z);
}

void BoundaryDirichlet_O4::applyAtPointStaggered(Field2D &f, BoutReal val, int x, int UNUSED(bx), int y, int UNUSED(by), int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = val;
}
void BoundaryDirichlet_O4::applyAtPointStaggered(Field3D &f, BoutReal val, int x, int UNUSED(bx), int y, int UNUSED(by), int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = val;
}

void BoundaryDirichlet_O4::extrapFurther(Field2D &f, int x, int bx, int y, int by, int z) {
  extrap4th(f, x, bx, y, by, z);
}
void BoundaryDirichlet_O4::extrapFurther(Field3D &f, int x, int bx, int y, int by, int z) {
  extrap4th(f, x, bx, y, by, z);
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryDirichlet_smooth::clone(BoundaryRegion *region, const list<string> &args){
  verifyNumPoints(region, 2);

  std::shared_ptr<FieldGenerator> newgen = nullptr;
  if(!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryDirichlet_smooth(region, newgen);
}

void BoundaryDirichlet_smooth::applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = 5./3.*val - 0.5*f(x - bx, y - by, z) - 1./6.*f(x - 2*bx, y - 2*by, z);
}
void BoundaryDirichlet_smooth::applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  // Dirichlet bc using val and first grid point would be
  // fb = 2*val - f0
  // using val and second grid point would be
  // fb = 4/3*val - 1/3*f1
  // Here we apply the bc using the average of the two, to try and suppress
  // grid-scale oscillations at the boundary
  f(x, y, z) = 5./3.*val - 0.5*f(x - bx, y - by, z) - 1./6.*f(x - 2*bx, y - 2*by, z);
}

void BoundaryDirichlet_smooth::applyAtPointStaggered(Field2D &f, BoutReal val, int x, int UNUSED(bx), int y, int UNUSED(by), int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = val;
}
void BoundaryDirichlet_smooth::applyAtPointStaggered(Field3D &f, BoutReal val, int x, int UNUSED(bx), int y, int UNUSED(by), int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = val;
}
///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryDirichlet_2ndOrder::clone(BoundaryRegion *region, const list<string> &args) {
  output << "WARNING: Use of boundary condition \"dirichlet_2ndorder\" is deprecated!\n";
  output << "         Consider using \"dirichlet\" instead\n";
  verifyNumPoints(region, 2);

  std::shared_ptr<FieldGenerator> newgen = nullptr;
  if(!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryDirichlet_2ndOrder(region, newgen);
}

// Set (at 2nd order) the value at the mid-point between the guard cell and the grid cell to be val
// N.B. Only first guard cells (closest to the grid) should ever be used
void BoundaryDirichlet_2ndOrder::applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = 8./3.*val - 2.*f(x - bx, y - by, z) + 1./3.*f(x - 2*bx, y - 2*by, z);
}
void BoundaryDirichlet_2ndOrder::applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = 8./3.*val - 2.*f(x - bx, y - by, z) + 1./3.*f(x - 2*bx, y - 2*by, z);
}

void BoundaryDirichlet_2ndOrder::applyAtPointStaggered(Field2D &f, BoutReal val, int x, int UNUSED(bx), int y, int UNUSED(by), int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = val;
}
void BoundaryDirichlet_2ndOrder::applyAtPointStaggered(Field3D &f, BoutReal val, int x, int UNUSED(bx), int y, int UNUSED(by), int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = val;
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryDirichlet_O5::clone(BoundaryRegion *region, const list<string> &args) {
  verifyNumPoints(region, 4);

  std::shared_ptr<FieldGenerator> newgen = nullptr;
  if(!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryDirichlet_O5(region, newgen);
}

void BoundaryDirichlet_O5::applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = 128./35.*val - 4.*f(x - bx, y - by, z) + 2.*f(x - 2*bx, y - 2*by, z) - 4./5.*f(x - 3*bx, y - 3*by, z) + 1./7.*f(x - 4*bx, y - 4*by, z);
}
void BoundaryDirichlet_O5::applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = 128./35.*val - 4.*f(x - bx, y - by, z) + 2.*f(x - 2*bx, y - 2*by, z) - 4./5.*f(x - 3*bx, y - 3*by, z) + 1./7.*f(x - 4*bx, y - 4*by, z);
}

void BoundaryDirichlet_O5::applyAtPointStaggered(Field2D &f, BoutReal val, int x, int UNUSED(bx), int y, int UNUSED(by), int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = val;
}
void BoundaryDirichlet_O5::applyAtPointStaggered(Field3D &f, BoutReal val, int x, int UNUSED(bx), int y, int UNUSED(by), int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = val;
}

void BoundaryDirichlet_O5::extrapFurther(Field2D &f, int x, int bx, int y, int by, int z) {
  // Changing this extrapolation to not depend on val, so just using grid point
  // values. JTO 16/10/2018
  extrap5th(f, x, bx, y, by, z);
}
void BoundaryDirichlet_O5::extrapFurther(Field3D &f, int x, int bx, int y, int by, int z) {
  // Changing this extrapolation to not depend on val, so just using grid point
  // values. Not sure if this is the correct order... JTO 16/10/2018
  extrap5th(f, x, bx, y, by, z);
}


///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryNeumann_NonOrthogonal::clone(BoundaryRegion *region, const list<string> &args) {
  verifyNumPoints(region, 1);
  if(!args.empty()) {
    output << "WARNING: argument is set to BoundaryNeumann_NonOrthogonal\n";
    // First argument should be a value
    val = stringToReal(args.front());
    return new BoundaryNeumann_NonOrthogonal(region, val);
  }
  return new BoundaryNeumann_NonOrthogonal(region);
}

template<typename T>
void BoundaryNeumann_NonOrthogonal::applyTemplate(T &f, BoutReal UNUSED(t)) {
  Mesh* localmesh = f.getMesh();
  Coordinates *metric = f.getCoordinates();
  // Calculate derivatives for metric use
  localmesh->communicate(f);
  T dfdy = DDY(f);
  T dfdz = DDZ(f);
  // Loop over all elements and set equal to the next point in
  for(bndry->first(); !bndry->isDone(); bndry->next1d()) {
    // Interpolate (linearly) metrics to halfway between last cell and boundary cell
    BoutReal g11shift = 0.5*(metric->g11(bndry->x,bndry->y) + metric->g11(bndry->x-bndry->bx,bndry->y));
    BoutReal g12shift = 0.5*(metric->g12(bndry->x,bndry->y) + metric->g12(bndry->x-bndry->bx,bndry->y));
    BoutReal g13shift = 0.5*(metric->g13(bndry->x,bndry->y) + metric->g13(bndry->x-bndry->bx,bndry->y));
    // Have to use derivatives at last gridpoint instead of derivatives on boundary layer
    //   because derivative values don't exist in boundary region
    // NOTE: should be fixed to interpolate to boundary line
    for(int z=0;z<localmesh->LocalNz;z++) {
      BoutReal xshift = g12shift*dfdy(bndry->x-bndry->bx,bndry->y,z) 
        + g13shift*dfdz(bndry->x-bndry->bx,bndry->y,z);
      if(bndry->bx != 0 && bndry->by == 0) {
        // x boundaries only
        BoutReal delta = bndry->bx*metric->dx(bndry->x, bndry->y);
        f(bndry->x, bndry->y, z) = f(bndry->x - bndry->bx, bndry->y, z) + delta/g11shift*(val - xshift);
        if (bndry->width == 2){
          f(bndry->x + bndry->bx, bndry->y, z) = f(bndry->x - 2*bndry->bx, bndry->y, z) + 3.0*delta/g11shift*(val - xshift);
        }
      } else if(bndry->by != 0 && bndry->bx == 0) {
        // y boundaries only
        //   no need to shift this b/c we want parallel nuemann not theta
        BoutReal delta = bndry->by*metric->dy(bndry->x, bndry->y);
        f(bndry->x, bndry->y, z) = f(bndry->x, bndry->y - bndry->by, z) + delta*val;
        if (bndry->width == 2){
          f(bndry->x, bndry->y + bndry->by, z) = f(bndry->x, bndry->y - 2*bndry->by, z) + 3.0*delta*val;
        }
      } else {
        // set corners to zero
        f(bndry->x, bndry->y, z) = 0.0;
        if (bndry->width == 2){
          f(bndry->x + bndry->bx, bndry->y + bndry->by, z) = 0.0;
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryNeumann2::clone(BoundaryRegion *region, const list<string> &args) {
  output << "WARNING: Use of boundary condition \"neumann2\" is deprecated!\n";
  output << "         Consider using \"neumann\" instead\n";
  verifyNumPoints(region, 2);

  std::shared_ptr<FieldGenerator> newgen = nullptr;
  if(!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryNeumann2\n";
  }
  return new BoundaryNeumann2(region, newgen);
}

void BoundaryNeumann2::applyAtPoint(Field2D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = (4.*f(x - bx, y - by, z) - f(x - 2*bx, y - 2*by, z))/3.;
}
void BoundaryNeumann2::applyAtPoint(Field3D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = (4.*f(x - bx, y - by, z) - f(x - 2*bx, y - 2*by, z))/3.;
}

void BoundaryNeumann2::applyAtPointStaggered(Field2D &UNUSED(f), BoutReal UNUSED(val), int UNUSED(x), int UNUSED(bx), int UNUSED(y), int UNUSED(by), int UNUSED(z), Coordinates* UNUSED(metric)) {
  throw BoutException("BoundaryNeumann2 not implemented for staggered grids");
}
void BoundaryNeumann2::applyAtPointStaggered(Field3D &UNUSED(f), BoutReal UNUSED(val), int UNUSED(x), int UNUSED(bx), int UNUSED(y), int UNUSED(by), int UNUSED(z), Coordinates* UNUSED(metric)) {
  throw BoutException("BoundaryNeumann2 not implemented for staggered grids");
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryNeumann_2ndOrder::clone(BoundaryRegion *region, const list<string> &args) {
  output << "WARNING: Use of boundary condition \"neumann_2ndorder\" is deprecated!\n";
  output << "         Consider using \"neumann\" instead\n";
  verifyNumPoints(region, 1);

  std::shared_ptr<FieldGenerator> newgen = nullptr;
  if(!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryNeumann_2ndOrder(region, newgen);
}

void BoundaryNeumann_2ndOrder::applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) {
  BoutReal delta = bx*metric->dx(x, y) + by*metric->dy(x, y);
  f(x, y, z) = f(x - bx, y - by, z) + val*delta;
}
void BoundaryNeumann_2ndOrder::applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) {
  BoutReal delta = bx*metric->dx(x, y) + by*metric->dy(x, y);
  f(x, y, z) = f(x - bx, y - by, z) + val*delta;
}

void BoundaryNeumann_2ndOrder::applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) {
  BoutReal delta = bx*metric->dx(x, y) + by*metric->dy(x, y);
  f(x, y, z) = (4.*f(x - bx, y - by, z) - f(x - 2*bx, y - 2*by, z) + 2.*delta*val)/3.;
}
void BoundaryNeumann_2ndOrder::applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) {
  BoutReal delta = bx*metric->dx(x, y) + by*metric->dy(x, y);
  f(x, y, z) = (4.*f(x - bx, y - by, z) - f(x - 2*bx, y - 2*by, z) + 2.*delta*val)/3.;
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryNeumann::clone(BoundaryRegion *region, const list<string> &args){
  verifyNumPoints(region, 1);

  std::shared_ptr<FieldGenerator> newgen = nullptr;
  if(!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryNeumann(region, newgen);
}

void BoundaryNeumann::applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) {
  BoutReal delta = bx*metric->dx(x, y) + by*metric->dy(x, y);
  f(x, y, z) = f(x - bx, y - by, z) + delta*val;
}
void BoundaryNeumann::applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) {
  BoutReal delta = bx*metric->dx(x, y) + by*metric->dy(x, y);
  f(x, y, z) = f(x - bx, y - by, z) + delta*val;
}

// For staggered case need to apply slightly differently Use one-sided
// differencing. Cell is now on the boundary, so use one-sided differencing
void BoundaryNeumann::applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) {
  BoutReal delta = bx*metric->dx(x, y) + by*metric->dy(x, y);
  f(x, y, z) = (4.*f(x - bx, y - by, z) - f(x - 2*bx, y - 2*by, z) + 2.*delta*val)/3.;
}
void BoundaryNeumann::applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) {
  BoutReal delta = bx*metric->dx(x, y) + by*metric->dy(x, y);
  f(x, y, z) = (4.*f(x - bx, y - by, z) - f(x - 2*bx, y - 2*by, z) + 2.*delta*val)/3.;
}


///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryNeumann_O4::clone(BoundaryRegion *region, const list<string> &args){
  verifyNumPoints(region, 4);

  std::shared_ptr<FieldGenerator> newgen = nullptr;
  if(!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryNeumann_O4(region, newgen);
}

void BoundaryNeumann_O4::applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) {
  BoutReal delta = bx*metric->dx(x, y) + by*metric->dy(x, y);
  f(x, y, z) = 12.*delta*val/11. +
               ( 17.*f(x - bx, y - by, z) + 9.*f(x - 2*bx, y - 2*by, z)
                 - 5.*f(x - 3*bx, y - 3*by, z) + f(x - 4*bx, y - 4*by, z))/22.;
}
void BoundaryNeumann_O4::applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) {
  BoutReal delta = bx*metric->dx(x, y) + by*metric->dy(x, y);
  f(x, y, z) = 12.*delta*val/11. +
               ( 17.*f(x - bx, y - by, z) + 9.*f(x - 2*bx, y - 2*by, z)
                 - 5.*f(x - 3*bx, y - 3*by, z) + f(x - 4*bx, y - 4*by, z))/22.;
}

void BoundaryNeumann_O4::applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) {
  BoutReal delta = bx*metric->dx(x, y) + by*metric->dy(x, y);
  f(x, y, z) = 12./25.*(delta*val
                        + 4.*f(x - bx, y - by, z) - 3.*f(x - 2*bx, y - 2*by, z)
                        + 4./3.*f(x - 3*bx, y - 3*by, z) - 1./4.*f(x - 4*bx, y - 4*by, z));
}
void BoundaryNeumann_O4::applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) {
  BoutReal delta = bx*metric->dx(x, y) + by*metric->dy(x, y);
  f(x, y, z) = 12./25.*(delta*val
                        + 4.*f(x - bx, y - by, z) - 3.*f(x - 2*bx, y - 2*by, z)
                        + 4./3.*f(x - 3*bx, y - 3*by, z) - 1./4.*f(x - 4*bx, y - 4*by, z));
}

void BoundaryNeumann_O4::extrapFurther(Field2D &f, int x, int bx, int y, int by, int z) {
  extrap5th(f, x, bx, y, by, z);
}
void BoundaryNeumann_O4::extrapFurther(Field3D &f, int x, int bx, int y, int by, int z) {
  extrap5th(f, x, bx, y, by, z);
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryNeumann_4thOrder::clone(BoundaryRegion *region, const list<string> &args) {
  verifyNumPoints(region, 4);

  std::shared_ptr<FieldGenerator> newgen = nullptr;
  if(!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryNeumann_4thOrder(region, newgen);
}

void BoundaryNeumann_4thOrder::applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) {
  BoutReal delta = -(bx*metric->dx(x, y) + by*metric->dy(x, y));
  f(x, y, z) = 12.*delta/11.*val + 17./22.*f(x - bx, y - by, z) + 9./22.*f(x - 2*bx, y - 2*by, z) - 5./22.*f(x - 3*bx, y - 3*by, z) + 1./22.*f(x - 4*bx, y - 4*by, z);
}
void BoundaryNeumann_4thOrder::applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) {
  BoutReal delta = -(bx*metric->dx(x, y) + by*metric->dy(x, y));
  f(x, y, z) = 12.*delta/11.*val + 17./22.*f(x - bx, y - by, z) + 9./22.*f(x - 2*bx, y - 2*by, z) - 5./22.*f(x - 3*bx, y - 3*by, z) + 1./22.*f(x - 4*bx, y - 4*by, z);
}

void BoundaryNeumann_4thOrder::applyAtPointStaggered(Field2D &UNUSED(f), BoutReal UNUSED(val), int UNUSED(x), int UNUSED(bx), int UNUSED(y), int UNUSED(by), int UNUSED(z), Coordinates* UNUSED(metric)) {
  throw BoutException("BoundaryNeumann_4thOrder is not implemented for staggered grids.");
}
void BoundaryNeumann_4thOrder::applyAtPointStaggered(Field3D &UNUSED(f), BoutReal UNUSED(val), int UNUSED(x), int UNUSED(bx), int UNUSED(y), int UNUSED(by), int UNUSED(z), Coordinates* UNUSED(metric)) {
  throw BoutException("BoundaryNeumann_4thOrder is not implemented for staggered grids.");
}

void BoundaryNeumann_4thOrder::extrapFurther(Field2D &f, int x, int bx, int y, int by, int z) {
  // Changing this extrapolation to not depend on val, so just using grid point
  // values. Previously was:
  // f(x+bx,y+by,z) = -24.*delta*val + 27.*f(x,y,z) - 27.*f(x-bx,y-by,z) + f(x-2*bx,y-2*by,z); // The f(x-4*bx,y-4*by,z) term vanishes, so that this sets to zero the 4th order central difference first derivative at the point half way between the guard cell and the grid cell
  // - JTO 16/10/2018
  extrap5th(f, x, bx, y, by, z);
}
void BoundaryNeumann_4thOrder::extrapFurther(Field3D &f, int x, int bx, int y, int by, int z) {
  // Changing this extrapolation to not depend on val, so just using grid point
  // values. Previously was
  // f(x+bx,y+by,z) = -24.*delta*val + 27.*f(x,y,z) - 27.*f(x-bx,y-by,z) + f(x-2*bx,y-2*by,z); // The f(x-4*bx,y-4*by,z) term vanishes, so that this sets to zero the 4th order central difference first derivative at the point half way between the guard cell and the grid cell
  // - JTO 16/10/2018
  extrap5th(f, x, bx, y, by, z);
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryNeumannPar::clone(BoundaryRegion *region, const list<string> &args) {
  verifyNumPoints(region, 1);

  std::shared_ptr<FieldGenerator> newgen = nullptr;
  if(!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryNeumannPar\n";
  }
  return new BoundaryNeumannPar(region, newgen);
}

void BoundaryNeumannPar::applyAtPoint(Field2D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* metric) {
  // For each element, set equal to the next point in
  f(x, y, z) = f(x - bx, y - by, z)*sqrt(metric->g_22(x, y)/metric->g_22(x - bx, y - by));
}
void BoundaryNeumannPar::applyAtPoint(Field3D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* metric) {
  // For each element, set equal to the next point in
  f(x, y, z) = f(x - bx, y - by, z)*sqrt(metric->g_22(x, y)/metric->g_22(x - bx, y - by));
}

void BoundaryNeumannPar::applyAtPointStaggered(Field2D &UNUSED(f), BoutReal UNUSED(val), int UNUSED(x), int UNUSED(bx), int UNUSED(y), int UNUSED(by), int UNUSED(z), Coordinates* UNUSED(metric)) {
  throw BoutException("BoundaryNeumannPar is not implemented for staggered grids.");
}
void BoundaryNeumannPar::applyAtPointStaggered(Field3D &UNUSED(f), BoutReal UNUSED(val), int UNUSED(x), int UNUSED(bx), int UNUSED(y), int UNUSED(by), int UNUSED(z), Coordinates* UNUSED(metric)) {
  throw BoutException("BoundaryNeumannPar is not implemented for staggered grids.");
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryRobin::clone(BoundaryRegion *region, const list<string> &args) {
  verifyNumPoints(region, 1);
  BoutReal a = 0.5, b = 1.0, g = 0.;
  
  list<string>::const_iterator it = args.begin();
  
  if(it != args.end()) {
    // First argument is 'a'
    a = stringToReal(*it);
    it++;
    
    if(it != args.end()) {
      // Second is 'b'
      b = stringToReal(*it);
      it++;
      
      if(it != args.end()) {
	// Third is 'g'
	g = stringToReal(*it);
	it++;
	if(it != args.end()) {
	  output << "WARNING: BoundaryRobin takes maximum of 3 arguments. Ignoring extras\n";
	}
      }
    }
  }
  
  return new BoundaryRobin(region, a, b, g);
}

template<typename T>
void BoundaryRobin::applyTemplate(T &f, BoutReal UNUSED(t)) {
  if(fabs(bval) < 1.e-12) {
    for(bndry->first(); !bndry->isDone(); bndry->next())
      for(int z=0; z<f.getNz(); z++)
	f(bndry->x, bndry->y, z) = gval / aval;
  }else {
    Coordinates* metric = f.getCoordinates();
    for(bndry->first(); !bndry->isDone(); bndry->next()) {
      BoutReal delta = bndry->bx*metric->dx(bndry->x, bndry->y) + bndry->by*metric->dy(bndry->x, bndry->y);
      for(int z=0; z<f.getNz(); z++) {
        f(bndry->x, bndry->y, z) = f(bndry->x - bndry->bx, bndry->y - bndry->by, z) + (gval - aval*f(bndry->x - bndry->bx, bndry->y - bndry->by, z) ) * delta / bval;
      }
    }
  }
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryConstGradient::clone(BoundaryRegion *region, const list<string> &args) {
  verifyNumPoints(region, 2);

  std::shared_ptr<FieldGenerator> newgen = nullptr;
  if(!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryConstGradient\n";
  }
  return new BoundaryConstGradient(region, newgen);
}

void BoundaryConstGradient::applyAtPoint(Field2D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = 2.*f(x - bx, y - by, z) - f(x - 2*bx, y - 2*by, z);
}
void BoundaryConstGradient::applyAtPoint(Field3D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  f(x, y, z) = 2.*f(x - bx, y - by, z) - f(x - 2*bx, y - 2*by, z);
}

void BoundaryConstGradient::applyAtPointStaggered(Field2D &UNUSED(f), BoutReal UNUSED(val), int UNUSED(x), int UNUSED(bx), int UNUSED(y), int UNUSED(by), int UNUSED(z), Coordinates* UNUSED(metric)) {
  throw BoutException("BoundaryConstGradient is not implemented for staggered grids.");
}
void BoundaryConstGradient::applyAtPointStaggered(Field3D &UNUSED(f), BoutReal UNUSED(val), int UNUSED(x), int UNUSED(bx), int UNUSED(y), int UNUSED(by), int UNUSED(z), Coordinates* UNUSED(metric)) {
  throw BoutException("BoundaryConstGradient is not implemented for staggered grids.");
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryZeroLaplace::clone(BoundaryRegion *region, const list<string> &args) {
  verifyNumPoints(region, 2);
  if(!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryZeroLaplace\n";
  }
  return new BoundaryZeroLaplace(region);
}

void BoundaryZeroLaplace::apply(Field2D &f, BoutReal UNUSED(t)) {
  Coordinates *metric = f.getCoordinates();
  if((bndry->location != BNDRY_XIN) && (bndry->location != BNDRY_XOUT)) {
    // Can't apply this boundary condition to non-X boundaries
    throw BoutException("ERROR: Can't apply Zero Laplace condition to non-X boundaries\n");
  }
  // Constant X derivative
  int bx = bndry->bx;
  // Loop over the Y dimension
  for(bndry->first(); !bndry->isDone(); bndry->nextY()) {
    int x = bndry->x;
    int y = bndry->y;
    BoutReal g = (f(x-bx,y) - f(x-2*bx,y)) / metric->dx(x-bx,y);
    // Loop in X towards edge of domain
    do {
      f(x,y) = f(x-bx,y) + g*metric->dx(x,y);
      bndry->nextX();
      x = bndry->x; y = bndry->y;
    }while(!bndry->isDone());
  }
}

void BoundaryZeroLaplace::apply(Field3D &f, BoutReal UNUSED(t)) {
  Mesh* localmesh = f.getMesh();

  int ncz = localmesh->LocalNz;

  Coordinates *metric = f.getCoordinates();

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
    rfft(f(x - bx, y), localmesh->LocalNz, c0.begin());
    rfft(f(x - 2 * bx, y), localmesh->LocalNz, c1.begin());
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
        BoutReal kwave = jz * 2.0 * PI / metric->zlength(); // wavenumber in [rad^-1]
        c0[jz] *= exp(coef * kwave);                        // The decaying solution only
      }
      // Reverse FFT
      irfft(c0.begin(), localmesh->LocalNz, f(x, y));

      bndry->nextX();
      x = bndry->x;
      y = bndry->y;
    } while (!bndry->isDone());
  }
}

///////////////////////////////////////////////////////////////

BoundaryOp *BoundaryZeroLaplace2::clone(BoundaryRegion *region,
                                        const list<string> &args) {
  verifyNumPoints(region, 3);
  if (!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryZeroLaplace2\n";
  }
  return new BoundaryZeroLaplace2(region);
}

void BoundaryZeroLaplace2::apply(Field2D &f, BoutReal UNUSED(t)) {
  if ((bndry->location != BNDRY_XIN) && (bndry->location != BNDRY_XOUT)) {
    // Can't apply this boundary condition to non-X boundaries
    throw BoutException(
        "ERROR: Can't apply Zero Laplace condition to non-X boundaries\n");
  }

  Coordinates *metric = f.getCoordinates();

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
}

void BoundaryZeroLaplace2::apply(Field3D &f, BoutReal UNUSED(t)) {
  Mesh* localmesh = f.getMesh();

  int ncz = localmesh->LocalNz;

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
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryConstLaplace::clone(BoundaryRegion *region, const list<string> &args) {
  verifyNumPoints(region, 2);
  if(!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryConstLaplace\n";
  }
  return new BoundaryConstLaplace(region);
}

void BoundaryConstLaplace::apply(Field2D &f, BoutReal UNUSED(t)) {
  if((bndry->location != BNDRY_XIN) && (bndry->location != BNDRY_XOUT)) {
    // Can't apply this boundary condition to non-X boundaries
    throw BoutException("ERROR: Can't apply Zero Laplace condition to non-X boundaries\n");
  }
  
  // Constant X second derivative
  int bx = bndry->bx;
  // Loop over the Y dimension
  for(bndry->first(); !bndry->isDone(); bndry->nextY()) {
    int x = bndry->x;
    int y = bndry->y;
    // Calculate the Laplacian on the last point
    dcomplex la,lb,lc;
    laplace_tridag_coefs(x-2*bx, y, 0, la, lb, lc);
    dcomplex val = la*f(x-bx-1,y) + lb*f(x-2*bx,y) + lc*f(x-2*bx+1,y);
    // Loop in X towards edge of domain
    do {
      laplace_tridag_coefs(x-bx, y, 0, la, lb, lc);
      if(bx < 0) { // Lower X
	f(x,y) = ((val - lb*f(x-bx,y) + lc*f(x-2*bx,y)) / la).real();
      }else  // Upper X
	f(x,y) = ((val - lb*f(x-bx,y) + la*f(x-2*bx,y)) / lc).real();
      
      bndry->nextX();
      x = bndry->x; y = bndry->y;
    }while(!bndry->isDone());
  }
}

void BoundaryConstLaplace::apply(Field3D &f, BoutReal UNUSED(t)) {
  if((bndry->location != BNDRY_XIN) && (bndry->location != BNDRY_XOUT)) {
    // Can't apply this boundary condition to non-X boundaries
    throw BoutException("ERROR: Can't apply Zero Laplace condition to non-X boundaries\n");
  }
  
  Mesh* localmesh = f.getMesh();

  Coordinates *metric = f.getCoordinates();
  
  int ncz = localmesh->LocalNz;

  // Allocate memory
  Array<dcomplex> c0(ncz/2 + 1), c1(ncz/2 + 1), c2(ncz/2 + 1);
  
  int bx = bndry->bx;
  // Loop over the Y dimension
  for(bndry->first(); !bndry->isDone(); bndry->nextY()) {
    int x = bndry->x;
    int y = bndry->y;
    
    // Take FFT of last 3 points in domain
    rfft(f(x-bx,y), ncz, c0.begin());
    rfft(f(x-2*bx,y), ncz, c1.begin());
    rfft(f(x-3*bx,y), ncz, c2.begin());
    dcomplex k0lin = (c1[0] - c0[0])/metric->dx(x-bx,y); // for kz=0 solution
    
    // Calculate Delp2 on point MXG+1 (and put into c1)
    for(int jz=0;jz<=ncz/2;jz++) {
      dcomplex la,lb,lc;
      laplace_tridag_coefs(x-2*bx, y, jz, la, lb, lc);
      if(bx < 0) { // Inner X
	c1[jz] = la*c0[jz] + lb*c1[jz] + lc*c2[jz];
      }else { // Outer X
	c1[jz] = la*c2[jz] + lb*c1[jz] + lc*c0[jz];
      }
    }
    // Solve  metric->g11*d2f/dx2 - metric->g33*kz^2f = 0
    // Assume metric->g11, metric->g33 constant -> exponential growth or decay
    BoutReal xpos = 0.0;
    // Loop in X towards edge of domain
    do {
      // kz = 0 solution
      xpos -= metric->dx(x,y);
      c2[0] = c0[0] + k0lin*xpos + 0.5*c1[0]*xpos*xpos/metric->g11(x-bx,y);
      // kz != 0 solution
      BoutReal coef = -1.0*sqrt(metric->g33(x-bx,y) / metric->g11(x-bx,y))*metric->dx(x-bx,y);
      for(int jz=1;jz<=ncz/2;jz++) {
	BoutReal kwave=jz*2.0*PI/metric->zlength(); // wavenumber in [rad^-1]
	c0[jz] *= exp(coef*kwave); // The decaying solution only
	// Add the particular solution
	c2[jz] = c0[jz] - c1[jz]/(metric->g33(x-bx,y)*kwave*kwave); 
      }
      // Reverse FFT
      irfft(c2.begin(), ncz, f(x,y));
      
      bndry->nextX();
      x = bndry->x; y = bndry->y;
    }while(!bndry->isDone());
  }
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryDivCurl::clone(BoundaryRegion *region, const list<string> &args) {
  if(!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryDivCurl\n";
  }
  return new BoundaryDivCurl(region);
}

void BoundaryDivCurl::apply(Vector2D &UNUSED(f)) {
  throw BoutException("ERROR: DivCurl boundary not yet implemented for 2D vectors\n");
}

void BoundaryDivCurl::apply(Vector3D &var) {
  int jx, jy, jz, jzp, jzm;
  BoutReal tmp;
  
  Mesh* localmesh = var.x.getMesh();

  Coordinates *metric = localmesh->getCoordinates(var.getLocation());
  
  int ncz = localmesh->LocalNz;
  
  if(bndry->location != BNDRY_XOUT) {
    throw BoutException("ERROR: DivCurl boundary only works for outer X currently\n");
  }
  var.toCovariant();
  
  if(localmesh->xstart > 2) {
    throw BoutException("Error: Div = Curl = 0 boundary condition doesn't work for MXG > 2. Sorry\n");
  }

  jx = localmesh->xend+1;
  for(jy=1;jy<localmesh->LocalNy-1;jy++) {
    for(jz=0;jz<ncz;jz++) {
      jzp = (jz+1) % ncz;
      jzm = (jz - 1 + ncz) % ncz;

      // dB_y / dx = dB_x / dy
      
      // dB_x / dy
      tmp = (var.x(jx-1,jy+1,jz) - var.x(jx-1,jy-1,jz)) / (metric->dy(jx-1,jy-1) + metric->dy(jx-1,jy));
      
      var.y(jx,jy,jz) = var.y(jx-2,jy,jz) + (metric->dx(jx-2,jy) + metric->dx(jx-1,jy)) * tmp;
      if(localmesh->xstart == 2)
	// 4th order to get last point
	var.y(jx+1,jy,jz) = var.y(jx-3,jy,jz) + 4.*metric->dx(jx,jy)*tmp;
      
      // dB_z / dx = dB_x / dz
      
      tmp = (var.x(jx-1,jy,jzp) - var.x(jx-1,jy,jzm)) / (2.*metric->dz);
      
      var.z(jx,jy,jz) = var.z(jx-2,jy,jz) + (metric->dx(jx-2,jy) + metric->dx(jx-1,jy)) * tmp;
      if(localmesh->xstart == 2)
	var.z(jx+1,jy,jz) = var.z(jx-3,jy,jz) + 4.*metric->dx(jx,jy)*tmp;

      // d/dx( Jmetric->g11 B_x ) = - d/dx( Jmetric->g12 B_y + Jmetric->g13 B_z) 
      //                    - d/dy( JB^y ) - d/dz( JB^z )
	
      tmp = -( metric->J(jx,jy)*metric->g12(jx,jy)*var.y(jx,jy,jz) + metric->J(jx,jy)*metric->g13(jx,jy)*var.z(jx,jy,jz)
	       - metric->J(jx-2,jy)*metric->g12(jx-2,jy)*var.y(jx-2,jy,jz) + metric->J(jx-2,jy)*metric->g13(jx-2,jy)*var.z(jx-2,jy,jz) )
	/ (metric->dx(jx-2,jy) + metric->dx(jx-1,jy)); // First term (d/dx) using vals calculated above
      tmp -= (metric->J(jx-1,jy+1)*metric->g12(jx-1,jy+1)*var.x(jx-1,jy+1,jz) - metric->J(jx-1,jy-1)*metric->g12(jx-1,jy-1)*var.x(jx-1,jy-1,jz)
	      + metric->J(jx-1,jy+1)*metric->g22(jx-1,jy+1)*var.y(jx-1,jy+1,jz) - metric->J(jx-1,jy-1)*metric->g22(jx-1,jy-1)*var.y(jx-1,jy-1,jz)
	      + metric->J(jx-1,jy+1)*metric->g23(jx-1,jy+1)*var.z(jx-1,jy+1,jz) - metric->J(jx-1,jy-1)*metric->g23(jx-1,jy-1)*var.z(jx-1,jy-1,jz))
	/ (metric->dy(jx-1,jy-1) + metric->dy(jx-1,jy)); // second (d/dy)
      tmp -= (metric->J(jx-1,jy)*metric->g13(jx-1,jy)*(var.x(jx-1,jy,jzp) - var.x(jx-1,jy,jzm)) +
	      metric->J(jx-1,jy)*metric->g23(jx-1,jy)*(var.y(jx-1,jy,jzp) - var.y(jx-1,jy,jzm)) +
	      metric->J(jx-1,jy)*metric->g33(jx-1,jy)*(var.z(jx-1,jy,jzp) - var.z(jx-1,jy,jzm))) / (2.*metric->dz);
      
      var.x(jx,jy,jz) = ( metric->J(jx-2,jy)*metric->g11(jx-2,jy)*var.x(jx-2,jy,jz) + 
			  (metric->dx(jx-2,jy) + metric->dx(jx-1,jy)) * tmp ) / metric->J(jx,jy)*metric->g11(jx,jy);
      if(localmesh->xstart == 2)
	var.x(jx+1,jy,jz) = ( metric->J(jx-3,jy)*metric->g11(jx-3,jy)*var.x(jx-3,jy,jz) + 
			      4.*metric->dx(jx,jy)*tmp ) / metric->J(jx+1,jy)*metric->g11(jx+1,jy);
    }
  }
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryFree::clone(BoundaryRegion *region, const list<string> &UNUSED(args)) {
  return new BoundaryFree(region);
}

void BoundaryFree::apply(Field2D &UNUSED(f), BoutReal UNUSED(t)) {
  // Do nothing for free boundary
}

void BoundaryFree::apply(Field3D &UNUSED(f), BoutReal UNUSED(t)) {
  // Do nothing for free boundary
}

void BoundaryFree::apply_ddt(Field2D &UNUSED(f)) {
  // Do nothing for free boundary
}

void BoundaryFree::apply_ddt(Field3D &UNUSED(f)) {
  // Do nothing for free boundary
}

///////////////////////////////////////////////////////////////
// New free boundary implementation.  Uses last grid points to extrapolate into the guard cells.
// Written by L. Easy. 
///////////////////////////////////////////////////////////////

// 2nd order extrapolation:

BoundaryOp* BoundaryFree_O2::clone(BoundaryRegion *region, const list<string> &args){
  verifyNumPoints(region, 2);

  if(!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryFree_O2\n";
  }
  return new BoundaryFree_O2(region) ; 
}

void BoundaryFree_O2::applyAtPoint(Field2D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  extrap2nd(f, x, bx, y, by, z);
}
void BoundaryFree_O2::applyAtPoint(Field3D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  extrap2nd(f, x, bx, y, by, z);
}

void BoundaryFree_O2::applyAtPointStaggered(Field2D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  extrap2nd(f, x, bx, y, by, z);
}
void BoundaryFree_O2::applyAtPointStaggered(Field3D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  extrap2nd(f, x, bx, y, by, z);
}

//////////////////////////////////
// Third order extrapolation:
//////////////////////////////////
BoundaryOp* BoundaryFree_O3::clone(BoundaryRegion *region, const list<string> &args){
  verifyNumPoints(region, 3);

  if(!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryFree_O3\n";
  }
  return new BoundaryFree_O3(region) ; 
}

void BoundaryFree_O3::applyAtPoint(Field2D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  extrap3rd(f, x, bx, y, by, z);
}
void BoundaryFree_O3::applyAtPoint(Field3D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  extrap3rd(f, x, bx, y, by, z);
}

void BoundaryFree_O3::applyAtPointStaggered(Field2D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  extrap3rd(f, x, bx, y, by, z);
}
void BoundaryFree_O3::applyAtPointStaggered(Field3D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  extrap3rd(f, x, bx, y, by, z);
}

void BoundaryFree_O3::extrapFurther(Field2D &f, int x, int bx, int y, int by, int z) {
  extrap3rd(f, x, bx, y, by, z);
}
void BoundaryFree_O3::extrapFurther(Field3D &f, int x, int bx, int y, int by, int z) {
  extrap3rd(f, x, bx, y, by, z);
}

// Fourth order extrapolation:
BoundaryOp* BoundaryFree_O4::clone(BoundaryRegion *region, const list<string> &args){
  verifyNumPoints(region, 4);

  if(!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryFree_O4\n";
  }
  return new BoundaryFree_O4(region);
}

void BoundaryFree_O4::applyAtPoint(Field2D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  extrap4th(f, x, bx, y, by, z);
}
void BoundaryFree_O4::applyAtPoint(Field3D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  extrap4th(f, x, bx, y, by, z);
}

void BoundaryFree_O4::applyAtPointStaggered(Field2D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  extrap4th(f, x, bx, y, by, z);
}
void BoundaryFree_O4::applyAtPointStaggered(Field3D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  extrap4th(f, x, bx, y, by, z);
}

void BoundaryFree_O4::extrapFurther(Field2D &f, int x, int bx, int y, int by, int z) {
  extrap4th(f, x, bx, y, by, z);
}
void BoundaryFree_O4::extrapFurther(Field3D &f, int x, int bx, int y, int by, int z) {
  extrap4th(f, x, bx, y, by, z);
}

// Fifth order extrapolation:
BoundaryOp* BoundaryFree_O5::clone(BoundaryRegion *region, const list<string> &args){
  verifyNumPoints(region, 5);

  if(!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryFree_O5\n";
  }
  return new BoundaryFree_O5(region);
}

void BoundaryFree_O5::applyAtPoint(Field2D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  extrap5th(f, x, bx, y, by, z);
}
void BoundaryFree_O5::applyAtPoint(Field3D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  extrap5th(f, x, bx, y, by, z);
}

void BoundaryFree_O5::applyAtPointStaggered(Field2D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  extrap5th(f, x, bx, y, by, z);
}
void BoundaryFree_O5::applyAtPointStaggered(Field3D &f, BoutReal UNUSED(val), int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) {
  extrap5th(f, x, bx, y, by, z);
}

void BoundaryFree_O5::extrapFurther(Field2D &f, int x, int bx, int y, int by, int z) {
  extrap5th(f, x, bx, y, by, z);
}
void BoundaryFree_O5::extrapFurther(Field3D &f, int x, int bx, int y, int by, int z) {
  extrap5th(f, x, bx, y, by, z);
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryRelax::cloneMod(BoundaryOp *operation, const list<string> &args) {
  BoundaryRelax* result = new BoundaryRelax(operation, r);
  
  if(!args.empty()) {
    // First argument should be the rate
    BoutReal val = stringToReal(args.front());
    val = fabs(val); // Should always be positive
    result->r = val;
  }

  return result;
}

void BoundaryRelax::apply(Field2D &f, BoutReal t) {
  // Just apply the original boundary condition to f
  op->apply(f, t);
}

void BoundaryRelax::apply(Field3D &f, BoutReal t) {
  // Just apply the original boundary condition to f
  op->apply(f, t);
}

void BoundaryRelax::apply_ddt(Field2D &f) {
  TRACE("BoundaryRelax::apply_ddt(Field2D)");

  // Make a copy of f
  Field2D g = f;
  // Apply the boundary to g
  op->apply(g);
  
  bndry->first();
  
  // Set time-derivatives
  for(bndry->first(); !bndry->isDone(); bndry->next()) {
    ddt(f)(bndry->x, bndry->y) = r * (g(bndry->x, bndry->y) - f(bndry->x, bndry->y));
  }
}

void BoundaryRelax::apply_ddt(Field3D &f) {
  TRACE("BoundaryRelax::apply_ddt(Field3D)");
  
  Mesh* localmesh = f.getMesh();

  // Make a copy of f
  Field3D g = f; // NOTE: This is not very efficient... copying entire field
  // Apply the boundary to g
  op->apply(g);
  // Set time-derivatives
  for(bndry->first(); !bndry->isDone(); bndry->next())
    for(int z=0;z<localmesh->LocalNz;z++) {
      ddt(f)(bndry->x, bndry->y, z) = r * (g(bndry->x, bndry->y, z) - f(bndry->x, bndry->y, z));
    }
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryWidth::cloneMod(BoundaryOp *operation, const list<string> &args) {
  BoundaryWidth* result = new BoundaryWidth(operation, width);
  
  if(args.empty()) {
    output << "WARNING: BoundaryWidth expected 1 argument\n";
  }else {
    // First argument should be the rate
    int val = stringToInt(args.front());
    result->width = val;
  }
  
  return result;
}

void BoundaryWidth::apply(Field2D &f, BoutReal t) {
  // Pointer to boundary region shared between all BoundaryOp, BoundaryModifiers
  int oldwid = bndry->width;
  bndry->width = width;
  op->apply(f, t);
  bndry->width = oldwid;
}

void BoundaryWidth::apply(Field3D &f, BoutReal t) {
  int oldwid = bndry->width;
  bndry->width = width;
  op->apply(f, t);
  bndry->width = oldwid;
}
  
void BoundaryWidth::apply_ddt(Field2D &f) {
  int oldwid = bndry->width;
  bndry->width = width;
  op->apply_ddt(f);
  bndry->width = oldwid;
}
  
void BoundaryWidth::apply_ddt(Field3D &f) {
  int oldwid = bndry->width;
  bndry->width = width;
  op->apply_ddt(f);
  bndry->width = oldwid;
}

///////////////////////////////////////////////////////////////
BoundaryOp* BoundaryToFieldAligned::cloneMod(BoundaryOp *operation, const list<string> &args) {
  BoundaryToFieldAligned* result = new BoundaryToFieldAligned(operation);
  
  if(!args.empty()) {
    output << "WARNING: BoundaryToFieldAligned expected no argument\n";
    //Shouldn't we throw ?
  }
  
  return result;
}

void BoundaryToFieldAligned::apply(Field2D &f, BoutReal t) {
  op->apply(f, t);
}

void BoundaryToFieldAligned::apply(Field3D &f, BoutReal t) {
  Mesh* localmesh = f.getMesh();

  //NOTE: This is not very efficient... updating entire field
  f = localmesh->fromFieldAligned(f);

  // Apply the boundary to shifted field
  op->apply(f, t);

  //Shift back
  f = localmesh->toFieldAligned(f);

  //This is inefficient -- could instead use the shiftZ just in the bndry
  //but this is not portable to other parallel transforms -- we could instead
  //have a flag to define the region in which we want to apply to/fromFieldAligned
}
  
void BoundaryToFieldAligned::apply_ddt(Field2D &f) {
  op->apply_ddt(f);
}

void BoundaryToFieldAligned::apply_ddt(Field3D &f) {
  Mesh* localmesh = f.getMesh();

  f = localmesh->fromFieldAligned(f);
  ddt(f) = localmesh->fromFieldAligned(ddt(f));
  op->apply_ddt(f);
  ddt(f) = localmesh->toFieldAligned(ddt(f));
}


///////////////////////////////////////////////////////////////
BoundaryOp* BoundaryFromFieldAligned::cloneMod(BoundaryOp *operation, const list<string> &args) {
  BoundaryFromFieldAligned* result = new BoundaryFromFieldAligned(operation);
  
  if(!args.empty()) {
    output << "WARNING: BoundaryFromFieldAligned expected no argument\n";
    //Shouldn't we throw ?
  }
  
  return result;
}

void BoundaryFromFieldAligned::apply(Field2D &f, BoutReal t) {
  op->apply(f, t);
}

void BoundaryFromFieldAligned::apply(Field3D &f, BoutReal t) {
  Mesh* localmesh = f.getMesh();

  //NOTE: This is not very efficient... shifting entire field
  f = localmesh->toFieldAligned(f);

  // Apply the boundary to shifted field
  op->apply(f, t);

  //Shift back
  f = localmesh->fromFieldAligned(f);

  //This is inefficient -- could instead use the shiftZ just in the bndry
  //but this is not portable to other parallel transforms -- we could instead
  //have a flag to define the region in which we want to apply to/fromFieldAligned
}
  
void BoundaryFromFieldAligned::apply_ddt(Field2D &f) {
  op->apply_ddt(f);
}

void BoundaryFromFieldAligned::apply_ddt(Field3D &f) {
  Mesh* localmesh = f.getMesh();

  f = localmesh->toFieldAligned(f);
  ddt(f) = localmesh->toFieldAligned(ddt(f));
  op->apply_ddt(f);
  ddt(f) = localmesh->fromFieldAligned(ddt(f));
}
