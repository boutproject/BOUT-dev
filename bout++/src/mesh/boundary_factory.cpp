
#include "boundary_factory.h"
#include "utils.h"

#include <list>
#include <string>
using std::list;
using std::string;

BoundaryFactory* BoundaryFactory::instance = NULL;

BoundaryFactory::~BoundaryFactory()
{
  // Free any boundaries
  for(map<string, BoundaryOp*>::iterator it = opmap.begin(); it != opmap.end(); it++) {
    delete it->second;
  }
  for(map<string, BoundaryModifier*>::iterator it = modmap.begin(); it != modmap.end(); it++) {
    delete it->second;
  }
}

BoundaryFactory* BoundaryFactory::getInstance()
{
  if(instance == NULL) {
    // Create the singleton object
    instance = new BoundaryFactory();
  }
  return instance;
}

void BoundaryFactory::cleanup()
{
  if(instance == NULL)
    return;
  
  // Just delete the instance
  delete instance;
  instance = NULL;
}

BoundaryOp* BoundaryFactory::create(const string &name, BoundaryRegion *region)
{
  output <<  name;
  
  // Search for a string of the form: modifier(operation)  
  int pos = name.find('(');
  if(pos == string::npos) {
    // No more (opening) brackets. Should be a boundary operation
    // Need to strip whitespace
    
    if( (name == "null") || (name == "none") )
      return NULL;

    BoundaryOp *op = findBoundaryOp(trim(name));
    if(op == NULL) {
      output << "\tERROR: Could not find boundary condition '" << name << "'" << endl;
      return NULL;
    }
    
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
    return mod->clone(op, arglist);
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

BoundaryOp* BoundaryFactory::create(const char* name, BoundaryRegion *region)
{
  return create(string(name), region);
}

BoundaryOp* BoundaryFactory::createFromOptions(const string &varname, BoundaryRegion *region)
{
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
  
  /// First try looking for (var, region)
  string set;
  options.get<string>(varname, prefix+region->label, set, "");
  
  if(!set.empty())
    return create(set, region);
  
  /// Then (var, side)
  options.get<string>(varname, prefix+side, set, "");
  if(!set.empty())
    return create(set, region);
  
  /// Then (var, all)
  options.get<string>(varname, prefix+"all", set, "");
  if(!set.empty())
    return create(set, region);
  
  /// Then (all, region)
  options.get<string>("all", prefix+region->label, set, "");
  if(!set.empty())
    return create(set, region);
  
  /// Then (all, side)
  options.get<string>("all", prefix+side, set, "");
  if(!set.empty())
    return create(set, region);
  
  /// Then (all, all)
  options.get<string>("all", prefix+"all", set, "");
  if(!set.empty())
    return create(set, region);
  
  output << "NONE" << endl;
  return NULL;
}

BoundaryOp* BoundaryFactory::createFromOptions(const char* varname, BoundaryRegion *region)
{
  return createFromOptions(string(varname), region);
}

void BoundaryFactory::add(BoundaryOp* bop, const string &name)
{
  if( (findBoundaryMod(name) != NULL) || (findBoundaryOp(name) != NULL) ) {
    // error - already exists
    output << "ERROR: Trying to add an already existing boundary: " << name << endl;
    return;
  }
  opmap[lowercase(name)] = bop;
}

void BoundaryFactory::add(BoundaryOp* bop, const char *name)
{
  add(bop, string(name));
}

void BoundaryFactory::addMod(BoundaryModifier* bmod, const string &name)
{
  if( (findBoundaryMod(name) != NULL) || (findBoundaryOp(name) != NULL) ) {
    // error - already exists
    output << "ERROR: Trying to add an already existing boundary modifier: " << name << endl;
    return;
  }
  modmap[lowercase(name)] = bmod;
}

void BoundaryFactory::addMod(BoundaryModifier* bmod, const char *name)
{
  addMod(bmod, string(name));
}

BoundaryOp* BoundaryFactory::findBoundaryOp(const string &s)
{
  map<string,BoundaryOp*>::iterator it;
  it = opmap.find(lowercase(s));
  if(it == opmap.end())
    return NULL;
  return it->second;
}

BoundaryModifier* BoundaryFactory::findBoundaryMod(const string &s)
{
  map<string,BoundaryModifier*>::iterator it;
  it = modmap.find(lowercase(s));
  if(it == modmap.end())
    return NULL;
  return it->second;
}

// Strips leading and trailing spaces from a string
const string BoundaryFactory::trim(const string &str, const string &c)
{
  string s(str);
  // Find the first character position after excluding leading blank spaces
  size_t startpos = s.find_first_not_of(c);
  // Find the first character position from reverse af
  size_t endpos = s.find_last_not_of(c);

  // if all spaces or empty, then return an empty string
  if(( startpos == string::npos ) || ( endpos == string::npos ))
  {
    s = "";
  }
  else
  {
    s = s.substr(startpos, endpos-startpos+1);
  }
  return s;
}
