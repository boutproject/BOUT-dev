
#include "boundary_standard.h"
#include "globals.h"

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryDirichlet::clone(BoundaryRegion *region)
{
  return new BoundaryDirichlet(region);
}

void BoundaryDirichlet::apply(Field2D &f)
{
  // Just loop over all elements and set to zero
  for(bndry->first(); !bndry->isDone(); bndry->next())
    f[bndry->x][bndry->y] = 0.;
}

void BoundaryDirichlet::apply(Field3D &f)
{
  // Just loop over all elements and set to zero
  for(bndry->first(); !bndry->isDone(); bndry->next())
    for(int z=0;z<mesh->ngz;z++)
      f[bndry->x][bndry->y][z] = 0.;
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryNeumann::clone(BoundaryRegion *region)
{
  return new BoundaryNeumann(region);
}

void BoundaryNeumann::apply(Field2D &f)
{
  // Loop over all elements and set equal to the next point in
  for(bndry->first(); !bndry->isDone(); bndry->next())
    f[bndry->x][bndry->y] = f[bndry->x - bndry->bx][bndry->y - bndry->by];
}

void BoundaryNeumann::apply(Field3D &f)
{
  for(bndry->first(); !bndry->isDone(); bndry->next())
    for(int z=0;z<mesh->ngz;z++)
      f[bndry->x][bndry->y][z] = f[bndry->x - bndry->bx][bndry->y - bndry->by][z];
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryRelax::clone(BoundaryOp *operation)
{
  BoundaryRelax* result = new BoundaryRelax(r);
  result->op = operation;
  result->bndry = operation->bndry;
  
  return result;
}
  
void BoundaryRelax::apply(Field2D &f)
{
  // Just apply the original boundary condition to f
  op->apply(f);
}

void BoundaryRelax::apply(Field3D &f)
{
  // Just apply the original boundary condition to f
  op->apply(f);
}

void BoundaryRelax::apply_ddt(Field2D &f)
{
#ifdef CHECK
  msg_stack.push("BoundaryRelax::apply_ddt(Field2D)");
#endif

  // Make a copy of f
  Field2D g = f;
  // Apply the boundary to g
  op->apply(g);
  
  bndry->first();
  exit(0);
  // Set time-derivatives
  for(bndry->first(); !bndry->isDone(); bndry->next())
    ddt(f)[bndry->x][bndry->y] = ddt(f)[bndry->x - bndry->bx][bndry->y - bndry->by] 
      + r * (g[bndry->x][bndry->y] - f[bndry->x][bndry->y]);

#ifdef CHECK
  msg_stack.pop();
#endif
}

void BoundaryRelax::apply_ddt(Field3D &f)
{
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
      ddt(f)[bndry->x][bndry->y][z] = ddt(f)[bndry->x - bndry->bx][bndry->y - bndry->by][z]
	+ r * (g[bndry->x][bndry->y][z] - f[bndry->x][bndry->y][z]);
    }

#ifdef CHECK
  msg_stack.pop();
#endif
}

