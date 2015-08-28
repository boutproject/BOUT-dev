
#include <globals.hxx>
#include <field_data.hxx>
#include <boundary_factory.hxx>
#include <output.hxx>
#include <field_factory.hxx>

FieldData::FieldData() : boundaryIsCopy(false), boundaryIsSet(true) {

}

FieldData::~FieldData() {
  if(!boundaryIsCopy) {
    // Delete the boundary operations
    for(vector<BoundaryOp*>::iterator it = bndry_op.begin(); it != bndry_op.end(); it++)
      delete (*it);
  }
}

void FieldData::setBoundary(const string &name) {
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

  boundaryIsSet = true;
  boundaryIsCopy = false;
}

void FieldData::setBoundary(const string &region, BoundaryOp *op) {
  /// Get the mesh boundary regions
  vector<BoundaryRegion*> reg = mesh->getBoundaries();

  /// Find the region


  /// Find if we're replacing an existing boundary
  for(vector<BoundaryOp*>::iterator it = bndry_op.begin(); it != bndry_op.end(); it++) {
    if( (*it)->bndry == op->bndry ) {
      // Replacing this boundary
      output << "Replacing ";
    }
  }
}

void FieldData::copyBoundary(const FieldData &f) {
  bndry_op = f.bndry_op;
  boundaryIsCopy = true;
  boundaryIsSet = true;
}

//JMAD
void FieldData::addBndryFunction(FuncPtr userfunc, BndryLoc location){
  /// NOTE: This will allocate memory, which may never be free'd
  addBndryGenerator( new FieldFunction(userfunc), location );
}


void FieldData::addBndryGenerator(FieldGenerator* gen, BndryLoc location){
  if(location == BNDRY_ALL){
    vector<BoundaryRegion*> reg = mesh->getBoundaries();
    for(vector<BoundaryRegion*>::iterator it=reg.begin(); it != reg.end(); it++) {
      bndry_generator[(*it)->location] = gen;
    }
  }
  else{
    bndry_generator[location] = gen;
  }
}

FieldGenerator* FieldData::getBndryGenerator(BndryLoc location) {
  std::map<BndryLoc,FieldGenerator*>::iterator it = bndry_generator.find(location);
  if(it == bndry_generator.end())
    return 0;

  return it->second;
}
