
#include "field_data.h"
#include "boundary_factory.h"
#include "globals.h"

FieldData::~FieldData()
{
  // Delete the boundary operations
  for(vector<BoundaryOp*>::iterator it = bndry_op.begin(); it != bndry_op.end(); it++)
    delete (*it);
}

void FieldData::setBoundary(const string &name)
{
  /// Get the boundary factory (singleton)
  BoundaryFactory *bfact = BoundaryFactory::getInstance();
  
  /// Get the mesh boundary regions
  vector<BoundaryRegion*> reg = mesh->getBoundaries();

  output << "Setting boundary for variable " << name << endl;
  /// Loop over the mesh boundary regions
  for(vector<BoundaryRegion*>::iterator it=reg.begin(); it != reg.end(); it++) {
    BoundaryOp* op = bfact->createFromOptions(name, (*it));
    if(op != NULL)
      bndry_op.push_back(op);
    output << endl;
  }
}

void FieldData::setBoundary(const string &region, BoundaryOp *op)
{
  /// First find if we're replacing an existing boundary
  
  /// Get the mesh boundary regions
  vector<BoundaryRegion*> reg = mesh->getBoundaries();
  
  
}
