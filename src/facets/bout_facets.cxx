/**
 * @file	bout_facets.cxx
 *
 * @brief	Implementation of calling BOUT through FACETS
 *
 */

// std includes
#include <string>

// txbase includes
// #include <TxDebugExcept.h>

// interface includes
#include <bout_facets.hxx>

// bout includes
#include <bout.hxx>

BoutIfc::BoutIfc(int petsc, int mpi) {
  std::cerr << "BoutIfc ctor called with petsc = " << petsc <<
    ", and mpi = " << mpi << "." << std::endl;
}

BoutIfc::~BoutIfc() {
}

int
BoutIfc::setLogFile(const std::string& fname) {
/** UNIMPLEMENTED ***/
  return 1;
}

int
BoutIfc::initialize() {
  char* argv[5] = {(char*)"",(char*)"-d",(char*)".",(char*)"-f",(char*)paramFile.c_str()};
  return bout_init(5, argv);
}

int
BoutIfc::setMpiComm(long fint) {
  BoutComm::getInstance()->setComm(fint);
  return 0;
}

// Only quasi-implemented. Still need to break-out the reading of params from bout_init
int
BoutIfc::readParams(const std::string& name) {
  std::cerr << "paramFile: " << name << std::endl;
  paramFile = name;
  return 1;
}

int
BoutIfc::buildData() {
/** UNIMPLEMENTED ***/
  return 1;
}

int
BoutIfc::buildUpdaters() {
/** UNIMPLEMENTED ***/
  return 1;
}

int
BoutIfc::getRankOfInterface(const std::string& ifc, size_t& rank) const {
/** UNIMPLEMENTED ***/
  return 1;
}

int
BoutIfc::set0dDouble(const std::string& name, double val) {
/** UNIMPLEMENTED ***/
  return 1;
}

int
BoutIfc::get0dDouble(const std::string& name, double& val) const {
/** UNIMPLEMENTED ***/
  return 1;
}

int
BoutIfc::set0dInt(const std::string& name, int value) {
/*** UNIMPLEMENTED ***/
  return 1;
}

int
BoutIfc::get0dInt(const std::string& name, int& value) const {
/*** UNIMPLEMENTED ***/
  return 1;
}

int
BoutIfc::set0dString(const std::string& name, const std::string& val) {
/** UNIMPLEMENTED ***/
  return 1;
}

int
BoutIfc::get0dString(const std::string& name, std::string& val) const {
/** UNIMPLEMENTED ***/
  return 1;
}

int
BoutIfc::dumpToFile(const std::string& name) const {
/** UNIMPLEMENTED ***/
  return 1;
}

int
BoutIfc::restoreFromFile(const std::string& name) {
/** UNIMPLEMENTED ***/
  return 1;
}

int
BoutIfc::update(double t) {
  bout_run();
/** Need to try and catch ***/
  return 0;
}

int
BoutIfc::revert() {
/** UNIMPLEMENTED ***/
  return 1;
}

int
BoutIfc::complete() {
  bout_finish();
  return 0;
}

// This is physics from advect1d, just an example
// If different physics is wanted, place it here

// Evolving variables 
Field3D N, P; // Density, Pressure
Vector3D V;   // velocity

// parameters
BoutReal gamma_ratio;   // Ratio of specific heats
BoutReal nu;      // Viscosity
bool include_viscosity;

Vector2D g; // Acceleration

int physics_init(bool restarting)
{
  // 2D initial profiles
  Field2D N0, P0;
  Vector2D V0;
  BoutReal v0_multiply;

  // Read initial conditions

  mesh->get(N0, "density");
  mesh->get(P0, "pressure");
  V0.covariant = false; // Read contravariant components
  V.covariant = false; // Evolve contravariant components
  mesh->get(V0, "v");
  g.covariant = false;
  mesh->get(g, "g");

  // read options

  Options *options = Options::getRoot();
  options->getSection("gas");
  options->get("gamma", gamma_ratio, 5./3.);
  options->get("viscosity", nu, 0.1);
  options->get("include_viscosity", include_viscosity, false);
  options->get("v0_multiply", v0_multiply, 1.0);

  V0 *= v0_multiply;

  // Set evolving variables

  bout_solve(N, "density");
  bout_solve(P, "pressure");
  bout_solve(V, "v");

  if(!restarting) {
    // Set variables to these values (+ the initial perturbation)
    // NOTE: This must be after the calls to bout_solve
    N += N0;
    P += P0;
    V += V0;
  }

  return 0;
}

int physics_run(BoutReal t)
{
  // Run communications
  mesh->communicate(N,P,V);

  // Density

  ddt(N) = -V_dot_Grad(V, N) - N*Div(V);

  // Velocity 

  ddt(V) = -V_dot_Grad(V, V) - Grad(P)/N + g;

  if(include_viscosity) {
    // Add viscosity

    ddt(V).y += nu*Laplacian(V.y);
    ddt(V).z += nu*Laplacian(V.z);
  }

  // Pressure

  ddt(P) = -V_dot_Grad(V, P) - gamma_ratio*P*Div(V);

  return 0;
}

