
#include <globals.hxx>
#include <field_data.hxx>
#include <boundary_factory.hxx>
#include <output.hxx>
#include <field_factory.hxx>
#include "unused.hxx"

FieldData::~FieldData() {
  if(!boundaryIsCopy) {
    // Delete the boundary operations
    for(const auto& bndry : bndry_op)
      delete bndry;
  }
}

void FieldData::setBoundary(const string &name) {
  /// Get the boundary factory (singleton)
  BoundaryFactory *bfact = BoundaryFactory::getInstance();
  
  output_info << "Setting boundary for variable " << name << endl;

  /// Get rid of existing boundary ops
  for (auto &op : bndry_op) {
    delete op;
  }
  bndry_op.clear();

  /// Loop over the mesh boundary regions
  for(const auto& reg : getDataMesh()->getBoundaries()) {
    BoundaryOp* op = bfact->createFromOptions(name, reg.get());
    if (op != nullptr)
      bndry_op.push_back(op);
    output_info << endl;
  }

  /// Get rid of existing parallel boundary ops
  for (auto &op : bndry_op_par) {
    delete op;
  }
  bndry_op_par.clear();

  /// Loop over the mesh parallel boundary regions
  for(const auto& reg : getDataMesh()->getBoundariesPar()) {
    BoundaryOpPar* op = bfact->createFromOptions(name, reg.get());
    if (op != nullptr)
      bndry_op_par.push_back(op);
    output_info << endl;
  }

  boundaryIsSet = true;
  boundaryIsCopy = false;
}

void FieldData::setBoundary(const string &region, BoundaryOp *op) {
  
  output_info << "Setting " << region << " boundary for some variable" << endl;

  /// Find the region
  BoundaryRegion* region_ptr = nullptr;
  for (const auto& bndry : getDataMesh()->getBoundaries()) {
    if (bndry->label == region) {
      region_ptr = bndry.get();
    }
  }

  /// Find if we're replacing an existing boundary
  for(auto it = bndry_op.begin(); it != bndry_op.end(); it++) {
    if( (*it)->bndry == region_ptr ) {
      // Replacing this boundary

      output << "Replacing " << region_ptr->label<<endl;

      // free the memory pointed to by this element of the vector
      delete (*it);

      // remove this element from the vector
      bndry_op.erase(it);
      break;
    }
  }

  /// Create the new BoundaryOp
  BoundaryOp* new_op = op->clone(region_ptr, {}, {});

  bndry_op.push_back(new_op);
}

void FieldData::copyBoundary(const FieldData &f) {
  bndry_op = f.bndry_op;
  bndry_op_par = f.bndry_op_par;
  boundaryIsCopy = true;
  boundaryIsSet = true;
}

//JMAD
void FieldData::addBndryFunction(FuncPtr userfunc, BndryLoc location){
  /// NOTE: This will allocate memory, which may never be free'd
  addBndryGenerator(std::make_shared<FieldFunction>(userfunc), location);
}

void FieldData::addBndryGenerator(FieldGeneratorPtr gen, BndryLoc location) {
  if(location == BNDRY_ALL){
    for(const auto& reg : getDataMesh()->getBoundaries()) {
      bndry_generator[reg->location] = gen;
    }
  } else {
    bndry_generator[location] = gen;
  }
}

FieldGeneratorPtr FieldData::getBndryGenerator(BndryLoc location) {
  std::map<BndryLoc, FieldGeneratorPtr>::iterator it = bndry_generator.find(location);
  if(it == bndry_generator.end())
    return nullptr;

  return it->second;
}
