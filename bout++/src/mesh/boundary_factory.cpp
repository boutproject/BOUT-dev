
#include "boundary_factory.h"

BoundaryFactory* BoundaryFactory::instance = NULL;

BoundaryFactory* BoundaryFactory::getInstance()
{
  if(instance == NULL) {
    // Create the singleton object
    instance = new BoundaryFactory();
  }
  return instance;
}
  
BoundaryOp* BoundaryFactory::create(const string &name, BoundaryRegion *region)
{
  output << "BOUNDARY FACTORY: Creating '" << name << "'" << endl;
  
  // Search for a string of the form: modifier(operation)  
  int pos = name.find('(');
  if(pos == string::npos) {
    // No more (opening) brackets. Should be a boundary operation
    // Need to strip whitespace
    
    BoundaryOp *op = findBoundaryOp(name);
    if(op == NULL) {
      output << "ERROR: Could not find boundary condition '" << name << "'" << endl;
      return NULL;
    }
    
    // Clone the boundary operation, passing the region to operate over
    return op->clone(region); 
  }
  // Contains a bracket. Find the last bracket and remove
  int pos2 = name.rfind(')');
  if(pos2 == string::npos) {
    output << "WARNING: Unmatched brackets in boundary condition: " << name << endl; 
  }
  
  // Find the modifier
  string modname = name.substr(0,pos);
  output << "BOUNDARY FACTORY: Looking for modifier " << modname  << endl;
  BoundaryModifier *mod = findBoundaryMod(modname);
  if(mod == NULL) {
    output << "ERROR: Could not find boundary modifier '" << modname << "'" << endl;
    return NULL;
  }

  // Create the operation inside the brackets
  BoundaryOp *op = create(name.substr(pos+1, pos2-pos-1), region);
  if(op == NULL) {
    // some error - abort
    return NULL;
  }
  
  // Clone the modifier, passing in the operator as argument
  return mod->clone(op);
}

BoundaryOp* BoundaryFactory::create(const char* name, BoundaryRegion *region)
{
  return create(string(name), region);
}


void BoundaryFactory::add(BoundaryOp* bop, const string &name)
{
  if(findBoundaryOp(name) != NULL) {
    // error - already exists
    output << "ERROR: Trying to add an already existing boundary: " << name << endl;
    return;
  }
  opmap[name] = bop;
}

void BoundaryFactory::add(BoundaryOp* bop, const char *name)
{
  add(bop, string(name));
}

void BoundaryFactory::addMod(BoundaryModifier* bmod, const string &name)
{
  if(findBoundaryMod(name) != NULL) {
    // error - already exists
    output << "ERROR: Trying to add an already existing boundary modifier: " << name << endl;
    return;
  }
  modmap[name] = bmod;
}

void BoundaryFactory::addMod(BoundaryModifier* bmod, const char *name)
{
  addMod(bmod, string(name));
}

BoundaryOp* BoundaryFactory::findBoundaryOp(string s)
{
  map<string,BoundaryOp*>::iterator it;
  it = opmap.find(s);
  if(it == opmap.end())
    return NULL;
  return it->second;
}

BoundaryModifier* BoundaryFactory::findBoundaryMod(string s)
{
  map<string,BoundaryModifier*>::iterator it;
  it = modmap.find(s);
  if(it == modmap.end())
    return NULL;
  return it->second;
}

