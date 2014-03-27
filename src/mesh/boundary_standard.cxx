
#include <globals.hxx>
#include <boundary_standard.hxx>
#include <invert_laplace.hxx>
#include <fft.hxx>
#include <boutexception.hxx>
#include <utils.hxx>
#include <output.hxx>
#include <msg_stack.hxx>
#include <bout/constants.hxx>

// #define BOUNDARY_CONDITIONS_UPGRADE_EXTRAPOLATE_FOR_2ND_ORDER

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryDirichlet::clone(BoundaryRegion *region, const list<string> &args) {
  if(!args.empty()) {
    // First argument should be a value
    val = stringToReal(args.front());
    return new BoundaryDirichlet(region, val);
  }
  return new BoundaryDirichlet(region);
}

void BoundaryDirichlet::apply(Field2D &f) {
  // Just loop over all elements and set to the value
  for(bndry->first(); !bndry->isDone(); bndry->next())
    f[bndry->x][bndry->y] = val;
}

void BoundaryDirichlet::apply(Field3D &f) {
  // Just loop over all elements and set to the value
  for(bndry->first(); !bndry->isDone(); bndry->next())
    for(int z=0;z<mesh->ngz;z++)
      f[bndry->x][bndry->y][z] = val;
}

void BoundaryDirichlet::apply_ddt(Field2D &f) {
  Field2D *dt = f.timeDeriv();
  for(bndry->first(); !bndry->isDone(); bndry->next())
    (*dt)[bndry->x][bndry->y] = 0.; // Set time derivative to zero
}

void BoundaryDirichlet::apply_ddt(Field3D &f) {
  Field3D *dt = f.timeDeriv();
  for(bndry->first(); !bndry->isDone(); bndry->next())
    for(int z=0;z<mesh->ngz;z++)
      (*dt)[bndry->x][bndry->y][z] = 0.; // Set time derivative to zero
}

////////////JMAD Constructors
BndDirichlet_O2::BndDirichlet_O2(){
  bndfunc = NULL;
}

BndDirichlet_O2::BndDirichlet_O2(BoundaryRegion *region):BoundaryOp(region){
  BndDirichlet_O2();
}

BoundaryOp* BndDirichlet_O2::clone(BoundaryRegion *region, const list<string> &args){
  return new BndDirichlet_O2(region);
}

void BndDirichlet_O2::apply(Field2D &f){
  output << "Time t from physics_run must be passed to boundary operator\n BndDirichlet_O2 using Field2D::applyBoundary(BoutReal t); \n ";
  output << "applying boundary condition for t = 0.!!!!!\n";
  BndDirichlet_O2::apply(f,0.);

}

void BndDirichlet_O2::apply(Field2D &f,BoutReal t) {
  // Set (at 2nd order) the value at the mid-point between the guard cell and the grid cell to be val
  // N.B. Only first guard cells (closest to the grid) should ever be used


  /*
  BoutReal x,xb,y,yb;
  bndry->first();
  if (bndry->by == 0){//x-boundary
    if(bndry->bx == -1){ // inner boundary
      xb = 0.;
    }
    else{//outer boundary
      xb = (BoutReal)mesh->getMX()*mesh->dx[bndry->x][bndry->y]; // =Lx
      //printf("bndryx %i bndryy %i xb %f yb %f bndry->x-bndry->bx %i\n ",bndry->x,bndry->y,xb, yb,bndry->x-bndry->bx );
    }

    for(bndry->first(); !bndry->isDone(); bndry->next1d()){
      //y position
      y = (BoutReal)(bndry->y - bndry->width + 0.5)*mesh->dy[bndry->x][bndry->y];
      //RHS bndry_funcs is a map in f. Key is bdnry location, Value of the map is a function pointer which is here dereferencing the function pointer, and hence evaluate the function
      f(bndry->x,bndry->y) = 2* (*((f.bndry_funcs)[bndry->location]))(t,xb,y,0.)
        - f(bndry->x-bndry->bx,bndry->y);
      // printf("bndryx %i bndryy %i xb %f y %f bndry->x-bndry->bx %i bndry->y-bndry->by %i fghost %f zk %i\n ",bndry->x,bndry->y,xb, y,bndry->x-bndry->bx,bndry->y-bndry->by, f(bndry->x-bndry->bx,bndry->y-bndry->by,zk),zk );

    }
  }
  else{// y-boundary
    if(bndry->by == -1){ // inner boundary
      yb = 0.;
    }
    else{//outer boundary
      yb = (BoutReal)mesh->getMY()*mesh->dy[bndry->x][bndry->y]; // =Ly
    }

    for(bndry->first(); !bndry->isDone(); bndry->next1d()){
      x = (BoutReal)(bndry->x - bndry->width + 0.5)*mesh->dx[bndry->x][bndry->y];

      //RHS bndry_funcs is a map in f. Key is bdnry location, Value of the map is a function pointer which is here dereferencing the function pointer, and hence evaluate the function
      f(bndry->x,bndry->y) = 2* (*((f.bndry_funcs)[bndry->location]))(t,x,yb,0.)
        - f(bndry->x,bndry->y-bndry->by);

    }
  }
  //printf("odne apply bnd \n\n");
  */
  
  // Replacement code (BD):
  for(bndry->first(); !bndry->isDone(); bndry->next1d()) {
    // Calculate the X and Y normalised values half-way between the guard cell and grid cell 
    BoutReal xnorm = 0.5*(   mesh->GlobalX(bndry->x)  // In the guard cell
                           + mesh->GlobalX(bndry->x - bndry->bx) ); // the grid cell

    BoutReal ynorm = 0.5*(   mesh->GlobalY(bndry->y)  // In the guard cell
                           + mesh->GlobalY(bndry->y - bndry->by) ); // the grid cell
    
    //RHS bndry_funcs is a map in f. Key is bdnry location, Value of the map is a function pointer which is here dereferencing the function pointer, and hence evaluate the function
    f(bndry->x,bndry->y) = 2* (*((f.bndry_funcs)[bndry->location]))(t,xnorm,ynorm,0.)
      - f(bndry->x-bndry->bx, bndry->y-bndry->by);
  }
}


void BndDirichlet_O2::apply(Field3D &f) {
  //BndDirichlet_O2::apply(f,0.);
  output << "Time t from physics_run must be passed to boundary operator\n BndDirichlet_O2 using Field3D::applyBoundary(BoutReal t); \n ";
  output << "applying boundary condition for t = 0.!!!!!\n";
  BndDirichlet_O2::apply(f,0.);
}


void BndDirichlet_O2::apply(Field3D &f,BoutReal t) {
  // Set (at 2nd order) the value at the mid-point between the guard cell and the grid cell to be val
  // N.B. Only first guard cells (closest to the grid) should ever be used

  /*
  BoutReal x,xb,y,yb;
  bndry->first();
  //printf("MX %i My %i bx %i by %i x %i y %i \n",mesh->getMX(),mesh->getMY(),bndry->bx,bndry->by,bndry->x,bndry->y);
  if (bndry->by == 0){//x-boundary
    if(bndry->bx == -1){ // inner boundary
      xb = 0.;
    }
    else{//outer boundary
      xb = (BoutReal)mesh->getMX()*mesh->dx[bndry->x][bndry->y]; // =Lx
      //printf("bndryx %i bndryy %i xb %f yb %f bndry->x-bndry->bx %i\n ",bndry->x,bndry->y,xb, yb,bndry->x-bndry->bx );
    }

    for(bndry->first(); !bndry->isDone(); bndry->next1d()){
      //y position
      y = (BoutReal)(bndry->y - bndry->width + 0.5)*mesh->dy[bndry->x][bndry->y];
      for(int zk=0;zk<mesh->ngz;zk++) {
        //RHS bndry_funcs is a map in f. Key is bdnry location, Value of the map is a function pointer which is here dereferencing the function pointer, and hence evaluate the function
        f(bndry->x,bndry->y,zk) = 2* (*((f.bndry_funcs)[bndry->location]))(t,xb,y,zk*mesh->dz)
          - f(bndry->x-bndry->bx,bndry->y,zk);
        // printf("bndryx %i bndryy %i xb %f y %f bndry->x-bndry->bx %i bndry->y-bndry->by %i fghost %f zk %i\n ",bndry->x,bndry->y,xb, y,bndry->x-bndry->bx,bndry->y-bndry->by, f(bndry->x-bndry->bx,bndry->y-bndry->by,zk),zk );
      }
    }
  }
  else{// y-boundary
    if(bndry->by == -1){ // inner boundary
      yb = 0.;
    }
    else{//outer boundary
      yb = (BoutReal)mesh->getMY()*mesh->dy[bndry->x][bndry->y]; // =Ly
    }

    for(bndry->first(); !bndry->isDone(); bndry->next1d()){
      x = (BoutReal)(bndry->x - bndry->width + 0.5)*mesh->dx[bndry->x][bndry->y];
      for(int zk=0;zk<mesh->ngz;zk++) {
        //RHS bndry_funcs is a map in f. Key is bdnry location, Value of the map is a function pointer which is here dereferencing the function pointer, and hence evaluate the function
        f(bndry->x,bndry->y,zk) = 2* (*((f.bndry_funcs)[bndry->location]))(t,x,yb,(BoutReal)zk*mesh->dz)
          - f(bndry->x,bndry->y-bndry->by,zk);
      }
    }
  }
  //printf("odne apply bnd \n\n");
  */
  
  // Replacement code (BD):
  for(bndry->first(); !bndry->isDone(); bndry->next1d()) {
    // Calculate the X and Y normalised values half-way between the guard cell and grid cell 
    BoutReal xnorm = 0.5*(   mesh->GlobalX(bndry->x)  // In the guard cell
                           + mesh->GlobalX(bndry->x - bndry->bx) ); // the grid cell

    BoutReal ynorm = 0.5*(   mesh->GlobalY(bndry->y)  // In the guard cell
                           + mesh->GlobalY(bndry->y - bndry->by) ); // the grid cell
    
    //RHS bndry_funcs is a map in f. Key is bdnry location, Value of the map is a function pointer which is here dereferencing the function pointer, and hence evaluate the function
    for(int zk=0;zk<mesh->ngz-1;zk++) {
      f(bndry->x,bndry->y,zk) = 2* (*((f.bndry_funcs)[bndry->location]))(t,xnorm,ynorm,(BoutReal)zk*mesh->dz)
        - f(bndry->x-bndry->bx, bndry->y-bndry->by, zk);
    }
  }
}


void BndDirichlet_O2::apply_ddt(Field2D &f) {
  Field2D *dt = f.timeDeriv();
  for(bndry->first(); !bndry->isDone(); bndry->next())
    (*dt)(bndry->x,bndry->y) = 0.; // Set time derivative to zero
}

void BndDirichlet_O2::apply_ddt(Field3D &f) {
  Field3D *dt = f.timeDeriv();
  for(bndry->first(); !bndry->isDone(); bndry->next())
    for(int z=0;z<mesh->ngz;z++)
      (*dt)(bndry->x,bndry->y,z) = 0.; // Set time derivative to zero
}


/*default function*/
BoutReal default_func(BoutReal t, BoutReal  x, BoutReal  y, BoutReal  z){
  output << "defaul boundary function assigned \n";
  return 99999999.;
}


/////end JMAD




///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryDirichlet_2ndOrder::clone(BoundaryRegion *region, const list<string> &args) {
  if(!args.empty()) {
    // First argument should be a value
    val = stringToReal(args.front());
    return new BoundaryDirichlet_2ndOrder(region, val);
  }
  return new BoundaryDirichlet_2ndOrder(region);
}

void BoundaryDirichlet_2ndOrder::apply(Field2D &f) {
  // Set (at 2nd order) the value at the mid-point between the guard cell and the grid cell to be val
  // N.B. Only first guard cells (closest to the grid) should ever be used
  for(bndry->first(); !bndry->isDone(); bndry->next1d()) {
    // /* this would only give 1st order (first) derivatives? */    f(bndry->x,bndry->y) = 2.*val - f(bndry->x-bndry->bx,bndry->y-bndry->by);
    f(bndry->x,bndry->y) = 8./3.*val - 2.*f(bndry->x-bndry->bx,bndry->y-bndry->by) + 1./3.*f(bndry->x-2*bndry->bx,bndry->y-2*bndry->by);
#ifdef BOUNDARY_CONDITIONS_UPGRADE_EXTRAPOLATE_FOR_2ND_ORDER
    f(bndry->x+bndry->bx,bndry->y+bndry->by) = 3.*f(bndry->x,bndry->y) - 3.*f(bndry->x-bndry->bx,bndry->y-bndry->by) + f(bndry->x-2*bndry->bx,bndry->y-2*bndry->by);
#elif defined(CHECK)
    f(bndry->x+bndry->bx,bndry->y+bndry->by) = 1.e60;
#endif
  }
}

void BoundaryDirichlet_2ndOrder::apply(Field3D &f) {
  // Set (at 2nd order) the value at the mid-point between the guard cell and the grid cell to be val
  // N.B. Only first guard cells (closest to the grid) should ever be used
  for(bndry->first(); !bndry->isDone(); bndry->next1d())
    for(int z=0;z<mesh->ngz;z++) {
      // /* this would only give 1st order (first) derivatives? */      f(bndry->x,bndry->y,z) = 2.*val - f(bndry->x-bndry->bx,bndry->y-bndry->by,z);
      f(bndry->x,bndry->y,z) = 8./3.*val - 2.*f(bndry->x-bndry->bx,bndry->y-bndry->by,z) + 1./3.*f(bndry->x-2*bndry->bx,bndry->y-2*bndry->by,z);
#ifdef BOUNDARY_CONDITIONS_UPGRADE_EXTRAPOLATE_FOR_2ND_ORDER
      f(bndry->x+bndry->bx,bndry->y+bndry->by,z) = 3.*f(bndry->x,bndry->y,z) - 3.*f(bndry->x-bndry->bx,bndry->y-bndry->by,z) + f(bndry->x-2*bndry->bx,bndry->y-2*bndry->by,z);
#elif defined(CHECK)
      f(bndry->x+bndry->bx,bndry->y+bndry->by,z) = 1.e60;
#endif
    }
}

void BoundaryDirichlet_2ndOrder::apply_ddt(Field2D &f) {
  Field2D *dt = f.timeDeriv();
  for(bndry->first(); !bndry->isDone(); bndry->next())
    (*dt)(bndry->x,bndry->y) = 0.; // Set time derivative to zero
}

void BoundaryDirichlet_2ndOrder::apply_ddt(Field3D &f) {
  Field3D *dt = f.timeDeriv();
  for(bndry->first(); !bndry->isDone(); bndry->next())
    for(int z=0;z<mesh->ngz;z++)
      (*dt)(bndry->x,bndry->y,z) = 0.; // Set time derivative to zero
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryDirichlet_4thOrder::clone(BoundaryRegion *region, const list<string> &args) {
  if(!args.empty()) {
    // First argument should be a value
    val = stringToReal(args.front());
    return new BoundaryDirichlet_4thOrder(region, val);
  }
  return new BoundaryDirichlet_4thOrder(region);
}

void BoundaryDirichlet_4thOrder::apply(Field2D &f) {
  // Set (at 4th order) the value at the mid-point between the guard cell and the grid cell to be val
  for(bndry->first(); !bndry->isDone(); bndry->next1d()) {
    f(bndry->x,bndry->y) = 128./35.*val - 4.*f(bndry->x-bndry->bx,bndry->y-bndry->by) + 2.*f(bndry->x-2*bndry->bx,bndry->y-2*bndry->by) - 4./3.*f(bndry->x-3*bndry->bx,bndry->y-3*bndry->by) + 1./7.*f(bndry->x-4*bndry->bx,bndry->y-4*bndry->by);
    f(bndry->x+bndry->bx,bndry->y+bndry->by) = -128./5.*val + 9.*f(bndry->x,bndry->y) + 18.*f(bndry->x-bndry->bx,bndry->y-bndry->by) -4.*f(bndry->x-2*bndry->bx,bndry->y-2*bndry->by) + 3./5.*f(bndry->x-3*bndry->bx,bndry->y-3*bndry->by);
  }
}

void BoundaryDirichlet_4thOrder::apply(Field3D &f) {
  // Set (at 4th order) the value at the mid-point between the guard cell and the grid cell to be val
  for(bndry->first(); !bndry->isDone(); bndry->next1d())
    for(int z=0;z<mesh->ngz;z++) {
      f(bndry->x,bndry->y,z) = 128./35.*val - 4.*f(bndry->x-bndry->bx,bndry->y-bndry->by,z) + 2.*f(bndry->x-2*bndry->bx,bndry->y-2*bndry->by,z) - 4./3.*f(bndry->x-3*bndry->bx,bndry->y-3*bndry->by,z) + 1./7.*f(bndry->x-4*bndry->bx,bndry->y-4*bndry->by,z);
      f(bndry->x+bndry->bx,bndry->y+bndry->by,z) = -128./5.*val + 9.*f(bndry->x,bndry->y,z) + 18.*f(bndry->x-bndry->bx,bndry->y-bndry->by,z) -4.*f(bndry->x-2*bndry->bx,bndry->y-2*bndry->by,z) + 3./5.*f(bndry->x-3*bndry->bx,bndry->y-3*bndry->by,z);
    }
}

void BoundaryDirichlet_4thOrder::apply_ddt(Field2D &f) {
  Field2D *dt = f.timeDeriv();
  for(bndry->first(); !bndry->isDone(); bndry->next())
    (*dt)(bndry->x,bndry->y) = 0.; // Set time derivative to zero
}

void BoundaryDirichlet_4thOrder::apply_ddt(Field3D &f) {
  Field3D *dt = f.timeDeriv();
  for(bndry->first(); !bndry->isDone(); bndry->next())
    for(int z=0;z<mesh->ngz;z++)
      (*dt)(bndry->x,bndry->y,z) = 0.; // Set time derivative to zero
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryNeumann::clone(BoundaryRegion *region, const list<string> &args) {
  if(!args.empty()) {
    //output << "WARNING: Ignoring arguments to BoundaryNeumann\n";
    output << "WARNING: arguments is set to BoundaryNeumann None Zero Gradient\n";
    // First argument should be a value
    val = stringToReal(args.front());
    return new BoundaryNeumann(region, val);
  }
  return new BoundaryNeumann(region);
}

void BoundaryNeumann::apply(Field2D &f) {
  // Loop over all elements and set equal to the next point in
  for(bndry->first(); !bndry->isDone(); bndry->next())
    f[bndry->x][bndry->y] = f[bndry->x - bndry->bx][bndry->y - bndry->by] + val*(bndry->bx*mesh->dx[bndry->x][bndry->y]+bndry->by*mesh->dy[bndry->x][bndry->y]);
}

void BoundaryNeumann::apply(Field3D &f) {
  for(bndry->first(); !bndry->isDone(); bndry->next())
    for(int z=0;z<mesh->ngz;z++) {
      BoutReal delta = bndry->bx*mesh->dx[bndry->x][bndry->y]+bndry->by*mesh->dy[bndry->x][bndry->y];
      f[bndry->x][bndry->y][z] = f[bndry->x - bndry->bx][bndry->y - bndry->by][z] + val*delta;
    }
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryNeumann2::clone(BoundaryRegion *region, const list<string> &args) {
  if(!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryNeumann2\n";
  }
  return new BoundaryNeumann2(region);
}

void BoundaryNeumann2::apply(Field2D &f) {
  // Loop over all elements and use one-sided differences
  for(bndry->first(); !bndry->isDone(); bndry->next())
    f(bndry->x, bndry->y) = (4.*f(bndry->x - bndry->bx, bndry->y - bndry->by) - f(bndry->x - 2*bndry->bx, bndry->y - 2*bndry->by))/3.;
}

void BoundaryNeumann2::apply(Field3D &f) {
  for(bndry->first(); !bndry->isDone(); bndry->next())
    for(int z=0;z<mesh->ngz;z++)
      f(bndry->x, bndry->y, z) = (4.*f(bndry->x - bndry->bx, bndry->y - bndry->by, z) - f(bndry->x - 2*bndry->bx, bndry->y - 2*bndry->by, z))/3.;
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryNeumann_2ndOrder::clone(BoundaryRegion *region, const list<string> &args) {
  if(!args.empty()) {
    // First argument should be a value
    val = stringToReal(args.front());
    return new BoundaryNeumann_2ndOrder(region, val);
  }
  return new BoundaryNeumann_2ndOrder(region);
}

void BoundaryNeumann_2ndOrder::apply(Field2D &f) {
  // Set (at 2nd order) the gradient at the mid-point between the guard cell and the grid cell to be val
  // This sets the value of the co-ordinate derivative, i.e. DDX/DDY not Grad_par/Grad_perp.x
  // N.B. Only first guard cells (closest to the grid) should ever be used
  for(bndry->first(); !bndry->isDone(); bndry->next1d()) {
    f(bndry->x,bndry->y) = f(bndry->x-bndry->bx,bndry->y-bndry->by) + val*(bndry->bx*mesh->dx(bndry->x,bndry->y)+bndry->by*mesh->dy(bndry->x,bndry->y));
#ifdef BOUNDARY_CONDITIONS_UPGRADE_EXTRAPOLATE_FOR_2ND_ORDER
    f(bndry->x+bndry->bx,bndry->y+bndry->by) = 3.*f(bndry->x,bndry->y) - 3.*f(bndry->x-bndry->bx,bndry->y-bndry->by) + f(bndry->x-2*bndry->bx,bndry->y-2*bndry->by);
#elif defined(CHECK)
    f(bndry->x+bndry->bx,bndry->y+bndry->by) = 1.e60;
#endif
  }
}

void BoundaryNeumann_2ndOrder::apply(Field3D &f) {
  // Set (at 2nd order) the gradient at the mid-point between the guard cell and the grid cell to be val
  // This sets the value of the co-ordinate derivative, i.e. DDX/DDY not Grad_par/Grad_perp.x
  // N.B. Only first guard cells (closest to the grid) should ever be used
  for(bndry->first(); !bndry->isDone(); bndry->next1d())
    for(int z=0;z<mesh->ngz;z++) {
      BoutReal delta = bndry->bx*mesh->dx(bndry->x,bndry->y)+bndry->by*mesh->dy(bndry->x,bndry->y);
      f(bndry->x,bndry->y,z) = f(bndry->x-bndry->bx,bndry->y-bndry->by,z) + val*delta;
#ifdef BOUNDARY_CONDITIONS_UPGRADE_EXTRAPOLATE_FOR_2ND_ORDER
      f(bndry->x+bndry->bx,bndry->y+bndry->by,z) = 3.*f(bndry->x,bndry->y,z) - 3.*f(bndry->x-bndry->bx,bndry->y-bndry->by,z) + f(bndry->x-2*bndry->bx,bndry->y-2*bndry->by,z);
#elif defined(CHECK)
      f(bndry->x+bndry->bx,bndry->y+bndry->by,z) = 1.e60;
#endif
    }
}

void BoundaryNeumann_2ndOrder::apply_ddt(Field2D &f) {
  Field2D *dt = f.timeDeriv();
  for(bndry->first(); !bndry->isDone(); bndry->next())
    (*dt)(bndry->x,bndry->y) = 0.; // Set time derivative to zero
}

void BoundaryNeumann_2ndOrder::apply_ddt(Field3D &f) {
  Field3D *dt = f.timeDeriv();
  for(bndry->first(); !bndry->isDone(); bndry->next())
    for(int z=0;z<mesh->ngz;z++)
      (*dt)(bndry->x,bndry->y,z) = 0.; // Set time derivative to zero
}

///////////////////////////////////////////////////////////////

/////JMAD///
BndNeumann_O2::BndNeumann_O2(){
  bndfunc = NULL;
}

BoundaryOp* BndNeumann_O2::clone(BoundaryRegion *region, const list<string> &args){
  //return new BndNeumann_O2(region,(FuncPtr)default_func);
  return new BndNeumann_O2(region);
}

void BndNeumann_O2::apply(Field2D &f) {
  //BndDirichlet_O2::apply(f,0.);
  output << "Time t from physics_run must be passed to boundary operator\n BndNeumann_O2 using Field3D::applyBoundary(BoutReal t); \n ";
  output << "applying boundary condition for t = 0.!!!!!\n";
  BndNeumann_O2::apply(f,0.);
}


void BndNeumann_O2::apply(Field2D &f,BoutReal t) {
  // Set (at 2nd order) the value at the mid-point between the guard cell and the grid cell to be val
  // N.B. Only first guard cells (closest to the grid) should ever be used
  
  /*
  BoutReal x,xb,y,yb,delta;
  bndry->first();
  //printf("MX %i My %i bx %i by %i x %i y %i \n",mesh->getMX(),mesh->getMY(),bndry->bx,bndry->by,bndry->x,bndry->y);
  if (bndry->by == 0){//x-boundary
    if(bndry->bx == -1){ // inner boundary
      xb = 0.;
    }
    else{//outer boundary
      xb = (BoutReal)mesh->getMX()*mesh->dx[bndry->x][bndry->y]; // =Lx
      //printf("bndryx %i bndryy %i xb %f yb %f bndry->x-bndry->bx %i\n ",bndry->x,bndry->y,xb, yb,bndry->x-bndry->bx );
    }

    for(bndry->first(); !bndry->isDone(); bndry->next1d()){
      //y position
      y = (BoutReal)(bndry->y - bndry->width + 0.5)*mesh->dy[bndry->x][bndry->y];
      //RHS bndry_funcs is a map in f. Key is bdnry location, Value of the map is a function pointer which is here dereferencing the function pointer, and hence evaluate the function
      delta = bndry->bx*mesh->dx(bndry->x,bndry->y)+bndry->by*mesh->dy(bndry->x,bndry->y);
      f(bndry->x,bndry->y) = f(bndry->x-bndry->bx,bndry->y) + (*((f.bndry_funcs)[bndry->location]))(t,xb,y,0.)*delta;
      // printf("bndryx %i bndryy %i xb %f y %f bndry->x-bndry->bx %i bndry->y-bndry->by %i fghost %f zk %i\n ",bndry->x,bndry->y,xb, y,bndry->x-bndry->bx,bndry->y-bndry->by, f(bndry->x-bndry->bx,bndry->y-bndry->by,zk),zk );

    }
  }
  else{// y-boundary
    if(bndry->by == -1){ // inner boundary
      yb = 0.;
    }
    else{//outer boundary
      yb = (BoutReal)mesh->getMY()*mesh->dy[bndry->x][bndry->y]; // =Ly
    }
    for(bndry->first(); !bndry->isDone(); bndry->next1d()){
      x = (BoutReal)(bndry->x - bndry->width + 0.5)*mesh->dx[bndry->x][bndry->y];
      delta = bndry->bx*mesh->dx(bndry->x,bndry->y)+bndry->by*mesh->dy(bndry->x,bndry->y);
      f(bndry->x,bndry->y) = f(bndry->x,bndry->y-bndry->by) + (*((f.bndry_funcs)[bndry->location]))(t,x,yb,0.)*delta;

    }
  }
  */
  
  // Replacement code (BD):
  for(bndry->first(); !bndry->isDone(); bndry->next1d()) {
    // Calculate the X and Y normalised values half-way between the guard cell and grid cell 
    BoutReal xnorm = 0.5*(   mesh->GlobalX(bndry->x)  // In the guard cell
                           + mesh->GlobalX(bndry->x - bndry->bx) ); // the grid cell

    BoutReal ynorm = 0.5*(   mesh->GlobalY(bndry->y)  // In the guard cell
                           + mesh->GlobalY(bndry->y - bndry->by) ); // the grid cell
    
    BoutReal delta = bndry->bx*mesh->dx(bndry->x,bndry->y)+bndry->by*mesh->dy(bndry->x,bndry->y);

    //RHS bndry_funcs is a map in f. Key is bdnry location, Value of the map is a function pointer which is here dereferencing the function pointer, and hence evaluate the function
    f(bndry->x,bndry->y) = f(bndry->x-bndry->bx, bndry->y-bndry->by) + delta*(*((f.bndry_funcs)[bndry->location]))(t,xnorm,ynorm,0.);
  }
}


void BndNeumann_O2::apply(Field3D &f) {
  //BndDirichlet_O2::apply(f,0.);
  output << "Time t from physics_run must be passed to boundary operator\n BndNeumann_O2 using Field3D::applyBoundary(BoutReal t); \n ";
  output << "applying boundary condition for t = 0.!!!!!\n";
  BndNeumann_O2::apply(f,0.);
}


void BndNeumann_O2::apply(Field3D &f,BoutReal t) {
  /*
  BoutReal x,xb,y,yb,delta,val;
  bndry->first();
  //printf("MX %i My %i bx %i by %i x %i y %i \n",mesh->getMX(),mesh->getMY(),bndry->bx,bndry->by,bndry->x,bndry->y);
  if (bndry->by == 0){//x-boundary
    if(bndry->bx == -1){ // inner boundary
      xb = 0.;
    }
    else{//outer boundary
      xb = (BoutReal)mesh->getMX()*mesh->dx[bndry->x][bndry->y]; // =Lx
      //printf("bndryx %i bndryy %i xb %f yb %f bndry->x-bndry->bx %i\n ",bndry->x,bndry->y,xb, yb,bndry->x-bndry->bx );
    }

    for(bndry->first(); !bndry->isDone(); bndry->next1d()){
      //y position
      y = (BoutReal)(bndry->y - bndry->width + 0.5)*mesh->dy[bndry->x][bndry->y];
      for(int zk=0;zk<mesh->ngz;zk++) {
        //RHS bndry_funcs is a map in f. Key is bdnry location, Value of the map is a function pointer which is here dereferencing the function pointer, and hence evaluate the function
        //dx
        delta = bndry->bx*mesh->dx(bndry->x,bndry->y)+bndry->by*mesh->dy(bndry->x,bndry->y);
        //bnd value
        val = (*((f.bndry_funcs)[bndry->location]))(t,xb,y,(BoutReal)zk*mesh->dz);
        f(bndry->x,bndry->y,zk) = f(bndry->x-bndry->bx,bndry->y-bndry->by,zk) + val*delta;
        //printf("val %f delta %f bndryx %i bndryy %i xb %f y %f bndry->x-bndry->bx %i bndry->y-bndry->by %i fghost %f zk %i\n ",val,delta,bndry->x,bndry->y,xb, y,bndry->x-bndry->bx,bndry->y-bndry->by, f(bndry->x-bndry->bx,bndry->y-bndry->by,zk),zk );
      }
    }
  }
  else{// y-boundary
    if(bndry->by == -1){ // inner boundary
      yb = 0.;
    }
    else{//outer boundary
      yb = (BoutReal)mesh->getMY()*mesh->dy[bndry->x][bndry->y]; // =Ly
    }
    for(bndry->first(); !bndry->isDone(); bndry->next1d()){
      x = (BoutReal)(bndry->x - bndry->width + 0.5)*mesh->dx[bndry->x][bndry->y];
      for(int zk=0;zk<mesh->ngz;zk++) {
        delta = bndry->bx*mesh->dx(bndry->x,bndry->y)+bndry->by*mesh->dy(bndry->x,bndry->y);
        val = (*((f.bndry_funcs)[bndry->location]))(t,x,yb,(BoutReal)zk*mesh->dz);
        f(bndry->x,bndry->y,zk) = f(bndry->x-bndry->bx,bndry->y-bndry->by,zk) + val*delta;
      }
    }
  }
  */

  // Replacement code (BD):
  for(bndry->first(); !bndry->isDone(); bndry->next1d()) {
    // Calculate the X and Y normalised values half-way between the guard cell and grid cell 
    BoutReal xnorm = 0.5*(   mesh->GlobalX(bndry->x)  // In the guard cell
                           + mesh->GlobalX(bndry->x - bndry->bx) ); // the grid cell

    BoutReal ynorm = 0.5*(   mesh->GlobalY(bndry->y)  // In the guard cell
                           + mesh->GlobalY(bndry->y - bndry->by) ); // the grid cell
    
    BoutReal delta = bndry->bx*mesh->dx(bndry->x,bndry->y)+bndry->by*mesh->dy(bndry->x,bndry->y);

    for(int zk=0;zk<mesh->ngz-1;zk++) {
      //RHS bndry_funcs is a map in f. Key is bdnry location, Value of the map is a function pointer which is here dereferencing the function pointer, and hence evaluate the function
      f(bndry->x,bndry->y, zk) = f(bndry->x-bndry->bx, bndry->y-bndry->by, zk) + delta*(*((f.bndry_funcs)[bndry->location]))(t,xnorm,ynorm,(BoutReal)zk*mesh->dz);
    }
  }
}

void BndNeumann_O2::apply_ddt(Field2D &f) {
  Field2D *dt = f.timeDeriv();
  for(bndry->first(); !bndry->isDone(); bndry->next())
    (*dt)(bndry->x,bndry->y) = 0.; // Set time derivative to zero
}

void BndNeumann_O2::apply_ddt(Field3D &f) {
  Field3D *dt = f.timeDeriv();
  for(bndry->first(); !bndry->isDone(); bndry->next())
    for(int z=0;z<mesh->ngz;z++)
      (*dt)(bndry->x,bndry->y,z) = 0.; // Set time derivative to zero
}
////END JMAD////


BoundaryOp* BoundaryNeumann_4thOrder::clone(BoundaryRegion *region, const list<string> &args) {
  if(!args.empty()) {
    // First argument should be a value
    val = stringToReal(args.front());
    return new BoundaryNeumann_4thOrder(region, val);
  }
  return new BoundaryNeumann_4thOrder(region);
}

void BoundaryNeumann_4thOrder::apply(Field2D &f) {
  // Set (at 4th order) the gradient at the mid-point between the guard cell and the grid cell to be val
  // This sets the value of the co-ordinate derivative, i.e. DDX/DDY not Grad_par/Grad_perp.x
  for(bndry->first(); !bndry->isDone(); bndry->next1d()) {
    BoutReal delta = -(bndry->bx*mesh->dx(bndry->x,bndry->y)+bndry->by*mesh->dy(bndry->x,bndry->y));
    f(bndry->x,bndry->y) = 12.*delta/11.*val + 17./22.*f(bndry->x-bndry->bx,bndry->y-bndry->by) + 9./22.*f(bndry->x-2*bndry->bx,bndry->y-2*bndry->by) - 5./22.*f(bndry->x-3*bndry->bx,bndry->y-3*bndry->by) + 1./22.*f(bndry->x-4*bndry->bx,bndry->y-4*bndry->by);
    f(bndry->x+bndry->bx,bndry->y+bndry->by) = -24.*delta*val + 27.*f(bndry->x,bndry->y) - 27.*f(bndry->x-bndry->bx,bndry->y-bndry->by) + f(bndry->x-2*bndry->bx,bndry->y-2*bndry->by); // The f(bndry->x-4*bndry->bx,bndry->y-4*bndry->by) term vanishes, so that this sets to zero the 4th order central difference first derivative at the point half way between the guard cell and the grid cell
  }
}

void BoundaryNeumann_4thOrder::apply(Field3D &f) {
  // Set (at 4th order) the gradient at the mid-point between the guard cell and the grid cell to be val
  // This sets the value of the co-ordinate derivative, i.e. DDX/DDY not Grad_par/Grad_perp.x
  for(bndry->first(); !bndry->isDone(); bndry->next1d())
    for(int z=0;z<mesh->ngz;z++) {
      BoutReal delta = -(bndry->bx*mesh->dx(bndry->x,bndry->y)+bndry->by*mesh->dy(bndry->x,bndry->y));
      f(bndry->x,bndry->y,z) = 12.*delta/11.*val + 17./22.*f(bndry->x-bndry->bx,bndry->y-bndry->by,z) + 9./22.*f(bndry->x-2*bndry->bx,bndry->y-2*bndry->by,z) - 5./22.*f(bndry->x-3*bndry->bx,bndry->y-3*bndry->by,z) + 1./22.*f(bndry->x-4*bndry->bx,bndry->y-4*bndry->by,z);
      f(bndry->x+bndry->bx,bndry->y+bndry->by,z) = -24.*delta*val + 27.*f(bndry->x,bndry->y,z) - 27.*f(bndry->x-bndry->bx,bndry->y-bndry->by,z) + f(bndry->x-2*bndry->bx,bndry->y-2*bndry->by,z); // The f(bndry->x-4*bndry->bx,bndry->y-4*bndry->by,z) term vanishes, so that this sets to zero the 4th order central difference first derivative at the point half way between the guard cell and the grid cell
    }
}

void BoundaryNeumann_4thOrder::apply_ddt(Field2D &f) {
  Field2D *dt = f.timeDeriv();
  for(bndry->first(); !bndry->isDone(); bndry->next())
    (*dt)(bndry->x,bndry->y) = 0.; // Set time derivative to zero
}

void BoundaryNeumann_4thOrder::apply_ddt(Field3D &f) {
  Field3D *dt = f.timeDeriv();
  for(bndry->first(); !bndry->isDone(); bndry->next())
    for(int z=0;z<mesh->ngz;z++)
      (*dt)(bndry->x,bndry->y,z) = 0.; // Set time derivative to zero
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryNeumannPar::clone(BoundaryRegion *region, const list<string> &args) {
  if(!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryNeumann2\n";
  }
  return new BoundaryNeumannPar(region);
}


void BoundaryNeumannPar::apply(Field2D &f) {
  // Loop over all elements and set equal to the next point in
  for(bndry->first(); !bndry->isDone(); bndry->next())
    f(bndry->x, bndry->y) = f(bndry->x - bndry->bx, bndry->y - bndry->by)*sqrt(mesh->g_22(bndry->x, bndry->y)/mesh->g_22(bndry->x - bndry->bx, bndry->y - bndry->by));
}

void BoundaryNeumannPar::apply(Field3D &f) {
  for(bndry->first(); !bndry->isDone(); bndry->next())
    for(int z=0;z<mesh->ngz;z++)
      f(bndry->x,bndry->y,z) = f(bndry->x - bndry->bx,bndry->y - bndry->by,z)*sqrt(mesh->g_22(bndry->x, bndry->y)/mesh->g_22(bndry->x - bndry->bx, bndry->y - bndry->by));
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryRobin::clone(BoundaryRegion *region, const list<string> &args) {
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

void BoundaryRobin::apply(Field2D &f) {
  if(fabs(bval) < 1.e-12) {
    // No derivative term so just constant value
    for(bndry->first(); !bndry->isDone(); bndry->next())
      f[bndry->x][bndry->y] = gval / aval;
  }else {
    BoutReal sign = 1.;
    if( (bndry->bx < 0) || (bndry->by < 0))
      sign = -1.;
    for(bndry->first(); !bndry->isDone(); bndry->next())
      f[bndry->x][bndry->y] = f[bndry->x - bndry->bx][bndry->y - bndry->by] + sign*(gval - aval*f[bndry->x - bndry->bx][bndry->y - bndry->by] ) / bval;
  }
}

void BoundaryRobin::apply(Field3D &f) {
  if(fabs(bval) < 1.e-12) {
    for(bndry->first(); !bndry->isDone(); bndry->next())
      for(int z=0;z<mesh->ngz;z++)
	f[bndry->x][bndry->y][z] = gval / aval;
  }else {
    BoutReal sign = 1.;
    if( (bndry->bx < 0) || (bndry->by < 0))
      sign = -1.;
    for(bndry->first(); !bndry->isDone(); bndry->next())
      for(int z=0;z<mesh->ngz;z++)
	f[bndry->x][bndry->y][z] = f[bndry->x - bndry->bx][bndry->y - bndry->by][z] + sign*(gval - aval*f[bndry->x - bndry->bx][bndry->y - bndry->by][z] ) / bval;
  }
}

///////////////////////////////////////////////////////////////

void BoundaryConstGradient::apply(Field2D &f){
  // Loop over all elements and set equal to the next point in
  for(bndry->first(); !bndry->isDone(); bndry->next())
    f[bndry->x][bndry->y] = 2.*f[bndry->x - bndry->bx][bndry->y - bndry->by] - f[bndry->x - 2*bndry->bx][bndry->y - 2*bndry->by];
}

void BoundaryConstGradient::apply(Field3D &f) {
  for(bndry->first(); !bndry->isDone(); bndry->next())
    for(int z=0;z<mesh->ngz;z++)
      f[bndry->x][bndry->y][z] = 2.*f[bndry->x - bndry->bx][bndry->y - bndry->by][z] - f[bndry->x - 2*bndry->bx][bndry->y - 2*bndry->by][z];
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryConstGradient::clone(BoundaryRegion *region, const list<string> &args) {
  if(!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryConstGradient\n";
  }
  return new BoundaryConstGradient(region);
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryZeroLaplace::clone(BoundaryRegion *region, const list<string> &args) {
  if(!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryZeroLaplace\n";
  }
  return new BoundaryZeroLaplace(region);
}

void BoundaryZeroLaplace::apply(Field2D &f) {
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
    BoutReal g = (f[x-bx][y] - f[x-2*bx][y]) / mesh->dx[x-bx][y];
    // Loop in X towards edge of domain
    do {
      f[x][y] = f[x-bx][y] + g*mesh->dx[x][y];
      bndry->nextX();
      x = bndry->x; y = bndry->y;
    }while(!bndry->isDone());
  }
}

void BoundaryZeroLaplace::apply(Field3D &f) {
  static dcomplex *c0 = (dcomplex*) NULL, *c1;
  int ncz = mesh->ngz-1;
  
  if(c0 == (dcomplex*) NULL) {
    // allocate memory
    c0 = new dcomplex[ncz/2 + 1];
    c1 = new dcomplex[ncz/2 + 1];
  }
  
  if((bndry->location != BNDRY_XIN) && (bndry->location != BNDRY_XOUT)) {
    // Can't apply this boundary condition to non-X boundaries
    throw BoutException("ERROR: Can't apply Zero Laplace condition to non-X boundaries\n");
  }
  
  int bx = bndry->bx;
  // Loop over the Y dimension
  for(bndry->first(); !bndry->isDone(); bndry->nextY()) {
    // bndry->(x,y) is the first point in the boundary
    // bndry->(x-bx,y) is the last "real" point in the domain
    
    int x = bndry->x;
    int y = bndry->y;
    
    // Take FFT of last 2 points in domain
    ZFFT(f[x-bx][y], mesh->zShift[x-bx][y], c0);
    ZFFT(f[x-2*bx][y], mesh->zShift[x-2*bx][y], c1);
    c1[0] = c0[0] - c1[0]; // Only need gradient
    
    // Solve  mesh->g11*d2f/dx2 - mesh->g33*kz^2f = 0
    // Assume mesh->g11, mesh->g33 constant -> exponential growth or decay
    
    // Loop in X towards edge of domain
    do {
      // kz = 0 solution
      c0[0] += c1[0];  // Straight line
      
      // kz != 0 solution
      BoutReal coef = -1.0*sqrt(mesh->g33[x][y] / mesh->g11[x][y])*mesh->dx[x][y];
      for(int jz=1;jz<=ncz/2;jz++) {
	BoutReal kwave=jz*2.0*PI/mesh->zlength; // wavenumber in [rad^-1]
	c0[jz] *= exp(coef*kwave); // The decaying solution only
      }
      // Reverse FFT
      ZFFT_rev(c0, mesh->zShift[x][y], f[x][y]);
      
      bndry->nextX();
      x = bndry->x; y = bndry->y;
    }while(!bndry->isDone());
  }
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryZeroLaplace2::clone(BoundaryRegion *region, const list<string> &args) {
  if(!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryZeroLaplace2\n";
  }
  return new BoundaryZeroLaplace2(region);
}

void BoundaryZeroLaplace2::apply(Field2D &f) {
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
    BoutReal g = (f[x-bx][y] - f[x-2*bx][y]) / mesh->dx[x-bx][y];
    // Loop in X towards edge of domain
    do {
      f[x][y] = f[x-bx][y] + g*mesh->dx[x][y];
      bndry->nextX();
      x = bndry->x; y = bndry->y;
    }while(!bndry->isDone());
  }
}

void BoundaryZeroLaplace2::apply(Field3D &f) {
  static dcomplex *c0 = (dcomplex*) NULL, *c1, *c2;
  int ncz = mesh->ngz-1;

  if(c0 == (dcomplex*) NULL) {
    // allocate memory
    c0 = new dcomplex[ncz/2 + 1];
    c1 = new dcomplex[ncz/2 + 1];
    c2 = new dcomplex[ncz/2 + 1];
  }
  
  if((bndry->location != BNDRY_XIN) && (bndry->location != BNDRY_XOUT)) {
    // Can't apply this boundary condition to non-X boundaries
    throw BoutException("ERROR: Can't apply Zero Laplace condition to non-X boundaries\n");
  }
  
  int bx = bndry->bx;
  // Loop over the Y dimension
  for(bndry->first(); !bndry->isDone(); bndry->nextY()) {
    // bndry->(x,y) is the first point in the boundary
    // bndry->(x-bx,y) is the last "real" point in the domain
    
    int x = bndry->x;
    int y = bndry->y;
    
    // Take FFT of last 2 points in domain
    ZFFT(f[x-bx][y], mesh->zShift[x-bx][y], c1);
    ZFFT(f[x-2*bx][y], mesh->zShift[x-2*bx][y], c2);
    
    // Loop in X towards edge of domain
    do {
      for(int jz=0;jz<=ncz/2;jz++) {
        dcomplex la,lb,lc;
        laplace_tridag_coefs(x-bx, y, jz, la, lb, lc);
        if(bx > 0) {
          // Outer boundary
          swap(la, lc);
        }
        c0[jz] = -(lb*c1[jz] + lc*c2[jz]) / la;
        /*
          if((y == 2) && (x == 1)) {
          output.write("Bndry %d: (%d,%d)\n", bx, x, jz);
          output << "\t[" << la << ", " << lb << ", " << lc << "]\n";
          output << "\t[" << c0[jz] << ", " << c1[jz] << ", " << c2[jz] << "]\n";
          }
        */
      }
      // Reverse FFT
      ZFFT_rev(c0, mesh->zShift[x][y], f[x][y]);
      // cycle c0 -> c1 -> c2 -> c0
      dcomplex *tmp = c2; c2 = c1; c1 = c0; c0 = tmp;
      
      bndry->nextX();
      x = bndry->x; y = bndry->y;
    }while(!bndry->isDone());
  }
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryConstLaplace::clone(BoundaryRegion *region, const list<string> &args) {
  if(!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryConstLaplace\n";
  }
  return new BoundaryConstLaplace(region);
}

void BoundaryConstLaplace::apply(Field2D &f) {
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
    dcomplex val = la*f[x-bx-1][y] + lb*f[x-2*bx][y] + lc*f[x-2*bx+1][y];
    // Loop in X towards edge of domain
    do {
      laplace_tridag_coefs(x-bx, y, 0, la, lb, lc);
      if(bx < 0) { // Lower X
	f[x][y] = ((val - lb*f[x-bx][y] + lc*f[x-2*bx][y]) / la).Real();
      }else  // Upper X
	f[x][y] = ((val - lb*f[x-bx][y] + la*f[x-2*bx][y]) / lc).Real();
      
      bndry->nextX();
      x = bndry->x; y = bndry->y;
    }while(!bndry->isDone());
  }
}

void BoundaryConstLaplace::apply(Field3D &f) {
  if((bndry->location != BNDRY_XIN) && (bndry->location != BNDRY_XOUT)) {
    // Can't apply this boundary condition to non-X boundaries
    throw BoutException("ERROR: Can't apply Zero Laplace condition to non-X boundaries\n");
  }
  
  static dcomplex *c0 = (dcomplex*) NULL, *c1, *c2;
  int ncz = mesh->ngz-1;
  if(c0 == (dcomplex*) NULL) {
    //allocate memory
    c0 = new dcomplex[ncz/2 + 1];
    c1 = new dcomplex[ncz/2 + 1];
    c2 = new dcomplex[ncz/2 + 1];
  }
  
  int bx = bndry->bx;
  // Loop over the Y dimension
  for(bndry->first(); !bndry->isDone(); bndry->nextY()) {
    int x = bndry->x;
    int y = bndry->y;
    
    // Take FFT of last 3 points in domain
    ZFFT(f[x-bx][y], mesh->zShift[x-bx][y], c0);
    ZFFT(f[x-2*bx][y], mesh->zShift[x-2*bx][y], c1);
    ZFFT(f[x-3*bx][y], mesh->zShift[x-3*bx][y], c2);
    dcomplex k0lin = (c1[0] - c0[0])/mesh->dx[x-bx][y]; // for kz=0 solution
    
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
    // Solve  mesh->g11*d2f/dx2 - mesh->g33*kz^2f = 0
    // Assume mesh->g11, mesh->g33 constant -> exponential growth or decay
    BoutReal xpos = 0.0;
    // Loop in X towards edge of domain
    do {
      // kz = 0 solution
      xpos -= mesh->dx[x][y];
      c2[0] = c0[0] + k0lin*xpos + 0.5*c1[0]*xpos*xpos/mesh->g11[x-bx][y];
      // kz != 0 solution
      BoutReal coef = -1.0*sqrt(mesh->g33[x-bx][y] / mesh->g11[x-bx][y])*mesh->dx[x-bx][y];
      for(int jz=1;jz<=ncz/2;jz++) {
	BoutReal kwave=jz*2.0*PI/mesh->zlength; // wavenumber in [rad^-1]
	c0[jz] *= exp(coef*kwave); // The decaying solution only
	// Add the particular solution
	c2[jz] = c0[jz] - c1[jz]/(mesh->g33[x-bx][y]*kwave*kwave); 
      }
      // Reverse FFT
      ZFFT_rev(c2, mesh->zShift[x][y], f[x][y]);
      
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

void BoundaryDivCurl::apply(Vector2D &f) {
  throw BoutException("ERROR: DivCurl boundary not yet implemented for 2D vectors\n");
}

void BoundaryDivCurl::apply(Vector3D &var) {
  int jx, jy, jz, jzp, jzm;
  BoutReal tmp;
  
  int ncz = mesh->ngz-1;
  
  if(bndry->location != BNDRY_XOUT) {
    throw BoutException("ERROR: DivCurl boundary only works for outer X currently\n");
  }
  var.toCovariant();
  
  if(mesh->xstart > 2) {
    output.write("Error: Div = Curl = 0 boundary condition doesn't work for MXG > 2. Sorry\n");
    exit(1);
  }

  jx = mesh->xend+1;
  for(jy=1;jy<mesh->ngy-1;jy++) {
    for(jz=0;jz<ncz;jz++) {
      jzp = (jz+1) % ncz;
      jzm = (jz - 1 + ncz) % ncz;

      // dB_y / dx = dB_x / dy
      
      // dB_x / dy
      tmp = (var.x[jx-1][jy+1][jz] - var.x[jx-1][jy-1][jz]) / (mesh->dy[jx-1][jy-1] + mesh->dy[jx-1][jy]);
      
      var.y[jx][jy][jz] = var.y[jx-2][jy][jz] + (mesh->dx[jx-2][jy] + mesh->dx[jx-1][jy]) * tmp;
      if(mesh->xstart == 2)
	// 4th order to get last point
	var.y[jx+1][jy][jz] = var.y[jx-3][jy][jz] + 4.*mesh->dx[jx][jy]*tmp;

      // dB_z / dx = dB_x / dz

      tmp = (var.x[jx-1][jy][jzp] - var.x[jx-1][jy][jzm]) / (2.*mesh->dz);
      
      var.z[jx][jy][jz] = var.z[jx-2][jy][jz] + (mesh->dx[jx-2][jy] + mesh->dx[jx-1][jy]) * tmp;
      if(mesh->xstart == 2)
	var.z[jx+1][jy][jz] = var.z[jx-3][jy][jz] + 4.*mesh->dx[jx][jy]*tmp;

      // d/dx( Jmesh->g11 B_x ) = - d/dx( Jmesh->g12 B_y + Jmesh->g13 B_z) 
      //                    - d/dy( JB^y ) - d/dz( JB^z )
	
      tmp = -( mesh->J[jx][jy]*mesh->g12[jx][jy]*var.y[jx][jy][jz] + mesh->J[jx][jy]*mesh->g13[jx][jy]*var.z[jx][jy][jz]
	       - mesh->J[jx-2][jy]*mesh->g12[jx-2][jy]*var.y[jx-2][jy][jz] + mesh->J[jx-2][jy]*mesh->g13[jx-2][jy]*var.z[jx-2][jy][jz] )
	/ (mesh->dx[jx-2][jy] + mesh->dx[jx-1][jy]); // First term (d/dx) using vals calculated above
      tmp -= (mesh->J[jx-1][jy+1]*mesh->g12[jx-1][jy+1]*var.x[jx-1][jy+1][jz] - mesh->J[jx-1][jy-1]*mesh->g12[jx-1][jy-1]*var.x[jx-1][jy-1][jz]
	      + mesh->J[jx-1][jy+1]*mesh->g22[jx-1][jy+1]*var.y[jx-1][jy+1][jz] - mesh->J[jx-1][jy-1]*mesh->g22[jx-1][jy-1]*var.y[jx-1][jy-1][jz]
	      + mesh->J[jx-1][jy+1]*mesh->g23[jx-1][jy+1]*var.z[jx-1][jy+1][jz] - mesh->J[jx-1][jy-1]*mesh->g23[jx-1][jy-1]*var.z[jx-1][jy-1][jz])
	/ (mesh->dy[jx-1][jy-1] + mesh->dy[jx-1][jy]); // second (d/dy)
      tmp -= (mesh->J[jx-1][jy]*mesh->g13[jx-1][jy]*(var.x[jx-1][jy][jzp] - var.x[jx-1][jy][jzm]) +
	      mesh->J[jx-1][jy]*mesh->g23[jx-1][jy]*(var.y[jx-1][jy][jzp] - var.y[jx-1][jy][jzm]) +
	      mesh->J[jx-1][jy]*mesh->g33[jx-1][jy]*(var.z[jx-1][jy][jzp] - var.z[jx-1][jy][jzm])) / (2.*mesh->dz);
      
      var.x[jx][jy][jz] = ( mesh->J[jx-2][jy]*mesh->g11[jx-2][jy]*var.x[jx-2][jy][jz] + 
			    (mesh->dx[jx-2][jy] + mesh->dx[jx-1][jy]) * tmp ) / mesh->J[jx][jy]*mesh->g11[jx][jy];
      if(mesh->xstart == 2)
	var.x[jx+1][jy][jz] = ( mesh->J[jx-3][jy]*mesh->g11[jx-3][jy]*var.x[jx-3][jy][jz] + 
				4.*mesh->dx[jx][jy]*tmp ) / mesh->J[jx+1][jy]*mesh->g11[jx+1][jy];
    }
  }
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryFree::clone(BoundaryRegion *region, const list<string> &args) {
  if (region->location & BNDRY_XIN) mesh->freeboundary_xin = true;
  if (region->location & BNDRY_XOUT) mesh->freeboundary_xout = true;
  if (region->location & BNDRY_YDOWN) mesh->freeboundary_ydown = true;
  if (region->location & BNDRY_YUP) mesh->freeboundary_yup = true;
  
  if(!args.empty()) {
    // First argument should be a value
    val = stringToReal(args.front());
    return new BoundaryFree(region, val);
  }
  return new BoundaryFree(region);
}

void BoundaryFree::apply(Field2D &f) {
  // Do nothing for free boundary
}

void BoundaryFree::apply(Field3D &f) {
  // Do nothing for free boundary
}

void BoundaryFree::apply_ddt(Field2D &f) {
  // Do nothing for free boundary
}

void BoundaryFree::apply_ddt(Field3D &f) {
  // Do nothing for free boundary
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
  
void BoundaryRelax::apply(Field2D &f) {
  // Just apply the original boundary condition to f
  op->apply(f);
}

void BoundaryRelax::apply(Field3D &f) {
  // Just apply the original boundary condition to f
  op->apply(f);
}

void BoundaryRelax::apply_ddt(Field2D &f) {
#ifdef CHECK
  msg_stack.push("BoundaryRelax::apply_ddt(Field2D)");
#endif

  // Make a copy of f
  Field2D g = f;
  // Apply the boundary to g
  op->apply(g);
  
  bndry->first();
  
  // Set time-derivatives
  for(bndry->first(); !bndry->isDone(); bndry->next()) {
    /*
      BoutReal lim = r * (g[bndry->x][bndry->y] - f[bndry->x][bndry->y]);
      BoutReal val = ddt(f)[bndry->x - bndry->bx][bndry->y - bndry->by] + lim;
      if((val*lim > 0.) && (fabs(val) > fabs(lim)))
      val = lim;
    
      ddt(f)[bndry->x][bndry->y] = val;
    */
    ddt(f)[bndry->x][bndry->y] = r * (g[bndry->x][bndry->y] - f[bndry->x][bndry->y]);
  }

#ifdef CHECK
  msg_stack.pop();
#endif
}

void BoundaryRelax::apply_ddt(Field3D &f) {
#ifdef CHECK
  msg_stack.push("BoundaryRelax::apply_ddt(Field2D)");
#endif
  
  // Make a copy of f
  Field3D g = f; // NOTE: This is not very efficient... copying entire field
  // Apply the boundary to g
  op->apply(g);
  // Set time-derivatives
  for(bndry->first(); !bndry->isDone(); bndry->next())
    for(int z=0;z<mesh->ngz;z++) {
      /*
        BoutReal lim = r * (g[bndry->x][bndry->y][z] - f[bndry->x][bndry->y][z]);
        BoutReal val = ddt(f)[bndry->x - bndry->bx][bndry->y - bndry->by][z] + lim;
        if((val*lim > 0.) && (fabs(val) > fabs(lim)))
        val = lim;
         
        ddt(f)[bndry->x][bndry->y][z] = val;
      */
      ddt(f)[bndry->x][bndry->y][z] = r * (g[bndry->x][bndry->y][z] - f[bndry->x][bndry->y][z]);
    }

#ifdef CHECK
  msg_stack.pop();
#endif
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryShifted::cloneMod(BoundaryOp *operation, const list<string> &args) {
  BoundaryShifted* result = new BoundaryShifted(operation);
  
  if(!args.empty()) {
    output << "WARNING: Ignoring arguments to BoundaryShifted\n";
  }

  return result;
}
  
void BoundaryShifted::apply(Field2D &f) {
  op->apply(f); // Doesn't affect 2D boundary conditions
}

void BoundaryShifted::apply(Field3D &f) {
  if(mesh->ShiftXderivs && (mesh->ShiftOrder == 0)) {
    Field3D g = f.shiftZ(true); // Shift into orthogonal coordinates
    op->apply(g);               // Apply the boundary condition
    f = g.shiftZ(false);        // Shift back to field-aligned
  }else
    op->apply(f);
}

void BoundaryShifted::apply_ddt(Field2D &f) {
  op->apply_ddt(f);
}

void BoundaryShifted::apply_ddt(Field3D &f) {
  if(mesh->ShiftXderivs && (mesh->ShiftOrder == 0)) {
    // Shift both the field and its time-derivative into orthogonal coordinates
    // This is in case the boundary condition depends on both (e.g. relaxing)
    Field3D g = f.shiftZ(true);
    ddt(g) = ddt(f).shiftZ(true);

    op->apply_ddt(g);             // Apply the boundary condition
    
    // Shift the time-derivative back 
    ddt(f) = ddt(g).shiftZ(false);       // Shift back to field-aligned
  }else
    op->apply_ddt(f);
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

void BoundaryWidth::apply(Field2D &f) {
  // Pointer to boundary region shared between all BoundaryOp, BoundaryModifiers
  int oldwid = bndry->width;
  bndry->width = width;
  op->apply(f);
  bndry->width = oldwid;
}

void BoundaryWidth::apply(Field3D &f) {
  int oldwid = bndry->width;
  bndry->width = width;
  op->apply(f);
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

