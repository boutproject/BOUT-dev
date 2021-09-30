
#include <bout/mesh.hxx>
#include <globals.hxx>
#include <field_data.hxx>
#include <boundary_factory.hxx>
#include <output.hxx>
#include <field_factory.hxx>
#include "unused.hxx"

namespace bout {
/// Make sure \p location is a sensible value for \p mesh
///
/// Throws if checks are enabled and trying to use a staggered
/// location on a non-staggered mesh
CELL_LOC normaliseLocation(CELL_LOC location, Mesh* mesh) {
  AUTO_TRACE();

  // CELL_DEFAULT always means CELL_CENTRE
  if (location == CELL_DEFAULT) {
    return CELL_CENTRE;
  }

  // No mesh means we can't check if we're using staggered grids, so
  // we'll have to trust the user in this case. This can happen if
  // we're making a field before the global mesh has been initialised
  // -- probably not good, but possible.
  if (mesh == nullptr) {
    return location;
  }

  if (mesh->StaggerGrids) {
    if (location == CELL_VSHIFT) {
      throw BoutException(
          "Field: CELL_VSHIFT cell location only makes sense for vectors");
    }
    return location;
  }
#if CHECK > 0
  if (location != CELL_CENTRE) {
    throw BoutException("Field: Trying to set off-centre location on "
                        "non-staggered grid\n"
                        "         Did you mean to enable staggered grids?");
  }
#endif
  return CELL_CENTRE;
}
} // namespace bout

FieldData::FieldData(Mesh* localmesh, CELL_LOC location_in)
    : fieldmesh(localmesh == nullptr ? bout::globals::mesh : localmesh),
      location(bout::normaliseLocation(
          location_in, fieldmesh)) { // Need to check for nullptr again, because the
                                     // fieldmesh might still be
  // nullptr if the global mesh hasn't been initialized yet
  if (fieldmesh != nullptr) {
    // sets fieldCoordinates by getting Coordinates for our location from
    // fieldmesh
    getCoordinates();
  }
}

FieldData::FieldData(const FieldData& other) {
  copyBoundary(other);
  fieldmesh = other.fieldmesh;
  location = other.location;
  fieldCoordinates = other.fieldCoordinates;
}

FieldData::~FieldData() {
  if(!boundaryIsCopy) {
    // Delete the boundary operations
    for(const auto& bndry : bndry_op)
      delete bndry;
  }
}

FieldData& FieldData::operator=(const FieldData& other) {
  if (this == &other) {
    return *this;
  }

  // Note we don't copy the boundaries here!
  fieldmesh = other.fieldmesh;
  location = other.location;
  fieldCoordinates = other.fieldCoordinates;
  return *this;
}

FieldData& FieldData::operator=(FieldData&& other) noexcept {
  // Note we don't copy the boundaries here!
  fieldmesh = other.fieldmesh;
  location = other.location;
  fieldCoordinates = std::move(other.fieldCoordinates);
  return *this;
}

namespace {
// Mark all the boundary condition input parameters as being
// "conditionally used".
//
// We now (soft) error if there are any unused input parameters, but
// because the set of boundary regions can be different on different
// processors, we might not "use" (that is, read the value of) the
// corresponding input parameters on the other processors!
//
// This function is a bit of a workaround for the above: essentially
// marks all possible boundary condition inputs as having been used on
// each processor. This way, we should still catch invalid inputs and
// misspellings (for example, `bdry_oll`) while still allowing valid
// boundaries on other processors.
void markBoundariesAsConditionallyUsed(const Mesh* mesh, Options& field_section) {
  // These are the "dynamic" boundaries that depend on the particular grid
  for (const auto& possible_boundary : mesh->getPossibleBoundaries()) {
    field_section["bndry_" + possible_boundary].setConditionallyUsed();
  }

  // These aren't currently handled by
  // `Mesh::getPossibleboundaries`. Potential fix: make `BndryLoc` a
  // `BOUT_ENUM_CLASS` for easy conversion to `string`, and then also
  // add `region->location` to possible boundaries
  using namespace std::string_literals;
  for (const auto& standard_boundary :
       {"xin"s, "xout"s, "ydown"s, "yup"s, "par_yup"s, "par_ydown"s, "all"s}) {
    field_section["bndry_" + standard_boundary].setConditionallyUsed();
  }
}
} // namespace

void FieldData::setBoundary(const std::string& name) {
  /// Get the boundary factory (singleton)
  BoundaryFactory* bfact = BoundaryFactory::getInstance();

  // It's plausible `fieldmesh` hasn't actually been set yet. This
  // seems most likely to happen if the field was declared as a
  // global. A different fix would be to move the call to
  // `setBoundary` to after the field initialisation in `Solver::add`
  Mesh* mesh = getMesh();

  markBoundariesAsConditionallyUsed(mesh, Options::root()[name]);
  markBoundariesAsConditionallyUsed(mesh, Options::root()["all"]);

  output_info << "Setting boundary for variable " << name << endl;
  /// Loop over the mesh boundary regions
  for (const auto& reg : mesh->getBoundaries()) {
    auto* op = dynamic_cast<BoundaryOp*>(bfact->createFromOptions(name, reg));
    if (op != nullptr)
      bndry_op.push_back(op);
    output_info << endl;
  }

  /// Get the mesh boundary regions
  std::vector<BoundaryRegionPar*> par_reg = mesh->getBoundariesPar();
  /// Loop over the mesh parallel boundary regions
  for (const auto& reg : mesh->getBoundariesPar()) {
    auto* op = dynamic_cast<BoundaryOpPar*>(bfact->createFromOptions(name, reg));
    if (op != nullptr)
      bndry_op_par.push_back(op);
    output_info << endl;
  }

  boundaryIsSet = true;
  boundaryIsCopy = false;
}

void FieldData::copyBoundary(const FieldData &f) {
  bndry_op = f.bndry_op;
  bndry_op_par = f.bndry_op_par;
  boundaryIsCopy = true;
  boundaryIsSet = true;
}

//JMAD
void FieldData::addBndryFunction(FuncPtr userfunc, BndryLoc location){
  addBndryGenerator(std::make_shared<FieldFunction>(userfunc), location);
}

void FieldData::addBndryGenerator(FieldGeneratorPtr gen, BndryLoc location) {
  // It's plausible `fieldmesh` hasn't actually been set yet. This
  // seems most likely to happen if the field was declared as a
  // global. A different fix would be to move the call to
  // `setBoundary` to after the field initialisation in `Solver::add`
  Mesh* mesh = getMesh();

  if (location == BNDRY_ALL) {
    for (const auto& reg : mesh->getBoundaries()) {
      bndry_generator[reg->location] = gen;
    }
  } else {
    bndry_generator[location] = std::move(gen);
  }
}

FieldGeneratorPtr FieldData::getBndryGenerator(BndryLoc location) {
  auto it = bndry_generator.find(location);
  if(it == bndry_generator.end())
    return nullptr;

  return it->second;
}

Mesh* FieldData::getMesh() const {
  if (fieldmesh != nullptr) {
    return fieldmesh;
  }
  // Don't set fieldmesh=mesh here, so that fieldmesh==nullptr until
  // allocate() is called in one of the derived classes. fieldmesh==nullptr
  // indicates that some initialization that would be done in the
  // constructor if fieldmesh was a valid Mesh object still needs to be
  // done.
  return bout::globals::mesh;
}

FieldData& FieldData::setLocation(CELL_LOC new_location) {
  AUTO_TRACE();

  location = bout::normaliseLocation(new_location, getMesh());

  fieldCoordinates.reset();
  // Sets correct fieldCoordinates pointer and ensures Coordinates object is
  // initialized for this Field's location
  getCoordinates();

  return *this;
}

CELL_LOC FieldData::getLocation() const {
  AUTO_TRACE();
  return location;
}

Coordinates* FieldData::getCoordinates() const {
  auto fieldCoordinates_shared = fieldCoordinates.lock();
  if (fieldCoordinates_shared) {
    return fieldCoordinates_shared.get();
  }
  fieldCoordinates = getMesh()->getCoordinatesSmart(FieldData::getLocation());
  return fieldCoordinates.lock().get();
}

Coordinates* FieldData::getCoordinates(CELL_LOC loc) const {
  if (loc == CELL_DEFAULT) {
    return getCoordinates();
  }
  return getMesh()->getCoordinates(loc);
}

void swap(FieldData& first, FieldData& second) noexcept {
  using std::swap;
  swap(first.fieldmesh, second.fieldmesh);
  swap(first.fieldCoordinates, second.fieldCoordinates);
  swap(first.location, second.location);
  swap(first.bndry_op, second.bndry_op);
  swap(first.boundaryIsCopy, second.boundaryIsCopy);
  swap(first.boundaryIsSet, second.boundaryIsSet);
  swap(first.bndry_op_par, second.bndry_op_par);
  swap(first.bndry_generator, second.bndry_generator);
}
