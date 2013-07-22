
#include <globals.hxx>
#include <boundary_factory.hxx>
#include <boundary_standard.hxx>
#include <utils.hxx>

#include <list>
#include <string>
using std::list;
using std::string;

#include <output.hxx>

BoundaryFactory* BoundaryFactory::instance = NULL;

BoundaryFactory::BoundaryFactory() {
  add(new BoundaryDirichlet(), "dirichlet");
  add(new BoundaryNeumann(), "neumann");
  add(new BoundaryNeumann2(), "neumann2");
  add(new BoundaryNeumannPar(), "neumannpar");
  add(new BoundaryRobin(), "robin");
  add(new BoundaryConstGradient(), "constgradient");
  add(new BoundaryZeroLaplace(), "zerolaplace");
  add(new BoundaryZeroLaplace2(), "zerolaplace2");
  add(new BoundaryConstLaplace(), "constlaplace");
  addMod(new BoundaryRelax(), "relax");
  addMod(new BoundaryShifted(), "shifted");
  addMod(new BoundaryWidth(), "width");
}

BoundaryFactory::~BoundaryFactory() {
  // Free any boundaries
  for(map<string, BoundaryOp*>::iterator it = opmap.begin(); it != opmap.end(); it++) {
    delete it->second;
  }
  for(map<string, BoundaryModifier*>::iterator it = modmap.begin(); it != modmap.end(); it++) {
    delete it->second;
  }
}

BoundaryFactory* BoundaryFactory::getInstance() {
  if(instance == NULL) {
    // Create the singleton object
    instance = new BoundaryFactory();
  }
  return instance;
}

void BoundaryFactory::cleanup() {
  if(instance == NULL)
    return;
  
  // Just delete the instance
  delete instance;
  instance = NULL;
}

BoundaryOp* BoundaryFactory::create(const string &name, BoundaryRegion *region) {
  //output <<  name;
  
  // Search for a string of the form: modifier(operation)  
  int pos = name.find('(');
  if(pos == string::npos) {
    // No more (opening) brackets. Should be a boundary operation
    // Need to strip whitespace
    
    if( (name == "null") || (name == "none") )
      return NULL;

    BoundaryOp *op = findBoundaryOp(trim(name));
    if(op == NULL)
      throw BoutException("Could not find boundary condition '%s'",  name.c_str());
    
    // Clone the boundary operation, passing the region to operate over and an empty args list
    list<string> args;
    return op->clone(region, args); 
  }
  // Contains a bracket. Find the last bracket and remove
  int pos2 = name.rfind(')');
  if(pos2 == string::npos) {
    output << "\tWARNING: Unmatched brackets in boundary condition: " << name << endl; 
  }
  
  // Find the function name before the bracket
  string func = trim(name.substr(0,pos));
  // And the argument inside the bracket
  string arg = trim(name.substr(pos+1, pos2-pos-1));
  // Split the argument on commas
  list<string> arglist = strsplit(arg, ',');
  for(list<string>::iterator it=arglist.begin(); it != arglist.end(); it++) {
    // Trim each argument
    (*it) = trim(*it);
  }
  
  // Test if func is a modifier
  BoundaryModifier *mod = findBoundaryMod(func);
  if(mod != NULL) {
    // The first argument should be an operation
    BoundaryOp *op = create(arglist.front(), region);
    if(op == NULL)
      return NULL;
    
    // Remove the first element (name of operation)
    arglist.pop_front();
    
    // Clone the modifier, passing in the operator and remaining strings as argument
    return mod->cloneMod(op, arglist);
  }
  
  BoundaryOp *op = findBoundaryOp(trim(func));
  if(op != NULL) {
    // An operation with arguments
    return op->clone(region, arglist); 
  }
  // Otherwise nothing matches
  output << "  Boundary setting is neither an operation nor modifier: " << func << endl;
  
  return NULL;
}

BoundaryOp* BoundaryFactory::create(const char* name, BoundaryRegion *region) {
  return create(string(name), region);
}

BoundaryOp* BoundaryFactory::createFromOptions(const string &varname, BoundaryRegion *region) {
  if(region == NULL)
    return NULL;
  
  output << "\t" << region->label << " region: ";
  
  string prefix("bndry_");
  
  string side("all");
  switch(region->location) {
  case BNDRY_XIN: {
    side = "xin";
    break;
  }
  case BNDRY_XOUT: {
    side = "xout";
    break;
  }
  case BNDRY_YDOWN: {
    side = "ydown";
    break;
  }
  case BNDRY_YUP: {
    side = "yup";
    break;
  }
  }
  
  // Get options
  Options *options = Options::getRoot();
  
  // Get variable options
  Options *varOpts = options->getSection(varname);
  string set;
  
  /// First try looking for (var, region)
  if(varOpts->isSet(prefix+region->label)) {
    varOpts->get(prefix+region->label, set, "");
    return create(set, region);
  }
  
  /// Then (var, side)
  if(varOpts->isSet(prefix+side)) {
    varOpts->get(prefix+side, set, "");
    return create(set, region);
  }
  
  /// Then (var, all)
  if(varOpts->isSet(prefix+"all")) {
    varOpts->get(prefix+"all", set, "");
    return create(set, region);
  }
  
  // Get the "all" options
  varOpts = options->getSection("All");

  /// Then (all, region)
  if(varOpts->isSet(prefix+region->label)) {
    varOpts->get(prefix+region->label, set, "");
    return create(set, region);
  }
  
  /// Then (all, side)
  if(varOpts->isSet(prefix+side)) {
    varOpts->get(prefix+side, set, "");
    return create(set, region);
  }
  
  /// Then (all, all)
  varOpts->get(prefix+"all", set, "dirichlet");
  return create(set, region);
  // Defaults to Dirichlet conditions, to prevent undefined boundary
  // values. If a user want to override, specify "none" or "null"
}

BoundaryOp* BoundaryFactory::createFromOptions(const char* varname, BoundaryRegion *region) {
  return createFromOptions(string(varname), region);
}

void BoundaryFactory::add(BoundaryOp* bop, const string &name) {
  if( (findBoundaryMod(name) != NULL) || (findBoundaryOp(name) != NULL) ) {
    // error - already exists
    output << "ERROR: Trying to add an already existing boundary: " << name << endl;
    return;
  }
  opmap[lowercase(name)] = bop;
}

void BoundaryFactory::add(BoundaryOp* bop, const char *name) {
  add(bop, string(name));
}

void BoundaryFactory::addMod(BoundaryModifier* bmod, const string &name) {
  if( (findBoundaryMod(name) != NULL) || (findBoundaryOp(name) != NULL) ) {
    // error - already exists
    output << "ERROR: Trying to add an already existing boundary modifier: " << name << endl;
    return;
  }
  modmap[lowercase(name)] = bmod;
}

void BoundaryFactory::addMod(BoundaryModifier* bmod, const char *name) {
  addMod(bmod, string(name));
}

BoundaryOp* BoundaryFactory::findBoundaryOp(const string &s) {
  map<string,BoundaryOp*>::iterator it;
  it = opmap.find(lowercase(s));
  if(it == opmap.end())
    return NULL;
  return it->second;
}

BoundaryModifier* BoundaryFactory::findBoundaryMod(const string &s) {
  map<string,BoundaryModifier*>::iterator it;
  it = modmap.find(lowercase(s));
  if(it == modmap.end())
    return NULL;
  return it->second;
}
