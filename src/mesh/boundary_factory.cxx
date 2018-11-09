#include <globals.hxx>
#include <boundary_factory.hxx>
#include <boundary_standard.hxx>
#include <utils.hxx>

#include <list>
#include <string>
#include <map>
using std::list;
using std::string;

#include <output.hxx>

BoundaryFactory *BoundaryFactory::instance = nullptr;

BoundaryFactory::BoundaryFactory() {
  add(new BoundaryDirichlet(), "dirichlet");
  add(new BoundaryDirichlet(), "dirichlet_o2"); // Synonym for "dirichlet"
  add(new BoundaryDirichlet_2ndOrder(), "dirichlet_2ndorder"); // Deprecated
  add(new BoundaryDirichlet_O3(), "dirichlet_o3");
  add(new BoundaryDirichlet_O4(), "dirichlet_o4");
  add(new BoundaryDirichlet_O5(), "dirichlet_o5");
  add(new BoundaryDirichlet_O5(), "dirichlet_4thorder"); // Synonym for "dirichlet_o5"
  add(new BoundaryDirichlet_smooth(), "dirichlet_smooth");
  add(new BoundaryNeumann(), "neumann");
  add(new BoundaryNeumann(), "neumann_O2"); // Synonym for "neumann"
  add(new BoundaryNeumann2(), "neumann2"); // Deprecated
  add(new BoundaryNeumann_2ndOrder(), "neumann_2ndorder"); // Deprecated
  add(new BoundaryNeumann_4thOrder(), "neumann_4thorder"); // Deprecated: Less good version of neumann_O4
  add(new BoundaryNeumann_O4(), "neumann_O4");
  add(new BoundaryNeumannPar(), "neumannpar");
  add(new BoundaryNeumann_NonOrthogonal(), "neumann_nonorthogonal");
  add(new BoundaryRobin(), "robin");
  add(new BoundaryConstGradient(), "constgradient");
  add(new BoundaryZeroLaplace(), "zerolaplace");
  add(new BoundaryZeroLaplace2(), "zerolaplace2");
  add(new BoundaryConstLaplace(), "constlaplace");
  add(new BoundaryFree(), "free");
  add(new BoundaryFree_O2(), "free_o2");
  add(new BoundaryFree_O3(), "free_o3");
  add(new BoundaryFree_O4(), "free_o4");
  add(new BoundaryFree_O5(), "free_o5");
  
  addMod(new BoundaryRelax(), "relax");
  addMod(new BoundaryWidth(), "width");
  addMod(new BoundaryToFieldAligned(), "toFieldAligned");
  addMod(new BoundaryFromFieldAligned(), "fromFieldAligned");

  // Parallel boundaries
  add(new BoundaryOpPar_dirichlet(), "parallel_dirichlet");
  add(new BoundaryOpPar_dirichlet_O3(), "parallel_dirichlet_O3");
  add(new BoundaryOpPar_dirichlet_interp(), "parallel_dirichlet_interp");
  add(new BoundaryOpPar_neumann(), "parallel_neumann");
}

BoundaryFactory::~BoundaryFactory() {
  // Free any boundaries
  for (const auto &it : opmap) {
    delete it.second;
  }
  for (const auto &it : modmap) {
    delete it.second;
  }
  for (const auto &it : par_opmap) {
    delete it.second;
  }
}

BoundaryFactory* BoundaryFactory::getInstance() {
  if (instance == nullptr) {
    // Create the singleton object
    instance = new BoundaryFactory();
  }
  return instance;
}

void BoundaryFactory::cleanup() {
  if (instance == nullptr)
    return;

  // Just delete the instance
  delete instance;
  instance = nullptr;
}

template<typename T>
BoundaryRegionOp<T>* BoundaryFactory::create(const string &name, T* region) {

  // Search for a string of the form: modifier(operation)
  auto pos = name.find('(');
  if (pos == string::npos) {
    // No more (opening) brackets. Should be a boundary operation
    // Need to strip whitespace

    if( (name == "null") || (name == "none") )
      return nullptr;

    BoundaryRegionOp<T> *op = findBoundaryOp< BoundaryRegionOp<T> >(trim(name));
    if (op == nullptr) {
      throw BoutException("Could not find boundary condition '%s'",  name.c_str());
    }

    // Clone the boundary operation, passing the region to operate over,
    // an empty args list and empty keyword map
    list<string> args;
    return op->clone(region, args, {});
  }

  // Contains a bracket. Find the last bracket and remove
  auto pos2 = name.rfind(')');
  if(pos2 == string::npos) {
    output_warn << "\tWARNING: Unmatched brackets in boundary condition: " << name << endl;
  }

  // Find the function name before the bracket
  string func = trim(name.substr(0,pos));
  // And the argument inside the bracket
  string arg = trim(name.substr(pos+1, pos2-pos-1));
  // Split the argument on commas
  // NOTE: Commas could be part of sub-expressions, so
  //       need to take account of brackets
  list<string> arglist;
  std::map<std::string, std::string> keywords;
  int level = 0;
  int start = 0;
  for (string::size_type i = 0; i < arg.length(); i++) {
    switch (arg[i]) {
    case '(':
    case '[':
    case '<':
      level++;
      break;
    case ')':
    case ']':
    case '>':
      level--;
      break;
    case ',': {
      if (level == 0) {
        string s = arg.substr(start, i-start);

        // Check if s contains '=', and if so treat as a keyword
        auto poseq = s.find('=');
        if (poseq != string::npos) {
          keywords[trim(s.substr(0, poseq))] = trim(s.substr(poseq + 1));
        } else {
          // No '=', so a positional argument
          arglist.push_back(trim(s));
        }
        start = i + 1;
      }
      break;
    }
    };
  }
  string s = arg.substr(start);
  auto poseq = s.find('=');
  if (poseq != string::npos) {
    keywords[trim(s.substr(0,poseq))] = trim(s.substr(poseq+1));
  } else {
    // No '=', so a positional argument
    arglist.push_back(trim(s));
  }

  if (std::is_same<BoundaryRegionOp<T>, BoundaryOp>::value) {
    // Test if func is a modifier
    BoundaryModifier *mod = findBoundaryMod(func);
    if (mod != nullptr) {
      // The first argument should be an operation
      BoundaryRegionOp<T> *op = create(arglist.front(), region);
      if (op == nullptr) {
        return nullptr;
      }

      // Remove the first element (name of operation)
      arglist.pop_front();

      // Clone the modifier, passing in the operator and remaining strings as argument
      return mod->cloneMod(op, arglist);
    }

    BoundaryRegionOp<T> *op = findBoundaryOp< BoundaryRegionOp<T> >(trim(func));
    if (op != nullptr) {
      // An operation with arguments
      return op->clone(region, arglist, keywords);
    }
  } else {
    // Parallel boundary
    BoundaryRegionOp<T> *pop = findBoundaryOp< BoundaryRegionOp<T> >(trim(func));
    if (pop != nullptr) {
      // An operation with arguments
      return pop->clone(region, arglist, keywords);
    }
  }

  // Otherwise nothing matches
  throw BoutException("  Boundary setting is neither an operation nor modifier: %s\n",func.c_str());

  return nullptr;
}

template<typename T>
BoundaryRegionOp<T>* BoundaryFactory::create(const char* name, T* region) {
  return create(string(name), region);
}
// const char* version calls the string version, so this should instantiate
// both:
template
BoundaryOp* BoundaryFactory::create(const char* name, BoundaryRegion* region);
template
BoundaryOpPar* BoundaryFactory::create(const char* name, BoundaryRegionPar* region);

template<typename T>
BoundaryRegionOp<T>* BoundaryFactory::createFromOptions(const string &varname, T* region) {
  if (region == nullptr)
    return nullptr;

  output_info << "\t" << region->label << " region: ";

  string prefix("bndry_");

  string side;
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
  case BNDRY_PAR_FWD: {
    side = "par_yup";
    break;
  }
  case BNDRY_PAR_BKWD: {
    side = "par_ydown";
    break;
  }
  default: {
    side = "all";
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
  if(std::is_same<BoundaryRegionOp<T>, BoundaryOpPar>::value) {
    if(varOpts->isSet(prefix+"par_all")) {
      varOpts->get(prefix+"par_all", set, "");
      return create(set, region);
    }
  } else {
    if(varOpts->isSet(prefix+"all")) {
      varOpts->get(prefix+"all", set, "");
      return create(set, region);
    }
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
  if(std::is_same<BoundaryRegionOp<T>, BoundaryOpPar>::value) {
    // Different default for parallel boundary regions
    varOpts->get(prefix+"par_all", set, "parallel_dirichlet");
  } else {
    varOpts->get(prefix+"all", set, "dirichlet");
  }
  return create(set, region);
  // Defaults to Dirichlet conditions, to prevent undefined boundary
  // values. If a user want to override, specify "none" or "null"
}

template<typename T>
BoundaryRegionOp<T>* BoundaryFactory::createFromOptions(const char* varname, T* region) {
  return createFromOptions(string(varname), region);
}
// const char* version calls the string version, so this should instantiate
// both:
template
BoundaryOp* BoundaryFactory::createFromOptions(const char* name, BoundaryRegion* region);
template
BoundaryOpPar* BoundaryFactory::createFromOptions(const char* name, BoundaryRegionPar* region);

void BoundaryFactory::add(BoundaryOp* bop, const string &name) {
  if ((findBoundaryMod(name) != nullptr) || (findBoundaryOp<BoundaryOp>(name) != nullptr)) {
    // error - already exists
    output_error << "ERROR: Trying to add an already existing boundary: " << name << endl;
    return;
  }
  opmap[lowercase(name)] = bop;
}

void BoundaryFactory::add(BoundaryOp* bop, const char *name) {
  add(bop, string(name));
}

void BoundaryFactory::add(BoundaryOpPar* bop, const string &name) {
  if (findBoundaryOp<BoundaryOpPar>(name) != nullptr) {
    // error - already exists
    output_error << "ERROR: Trying to add an already existing boundary: " << name << endl;
    return;
  }
  par_opmap[lowercase(name)] = bop;
}

void BoundaryFactory::add(BoundaryOpPar* bop, const char *name) {
  add(bop, string(name));
}

void BoundaryFactory::addMod(BoundaryModifier* bmod, const string &name) {
  if ((findBoundaryMod(name) != nullptr) || (findBoundaryOp<BoundaryOp>(name) != nullptr)) {
    // error - already exists
    output_error << "ERROR: Trying to add an already existing boundary modifier: " << name << endl;
    return;
  }
  modmap[lowercase(name)] = bmod;
}

void BoundaryFactory::addMod(BoundaryModifier* bmod, const char *name) {
  addMod(bmod, string(name));
}

template<>
BoundaryOp* BoundaryFactory::findBoundaryOp<BoundaryOp>(const string &s) {
  map<string,BoundaryOp*>::iterator it;
  it = opmap.find(lowercase(s));
  if(it == opmap.end())
    return nullptr;
  return it->second;
}

template<>
BoundaryOpPar* BoundaryFactory::findBoundaryOp<BoundaryOpPar>(const string &s) {
  map<string,BoundaryOpPar*>::iterator it;
  it = par_opmap.find(lowercase(s));
  if(it == par_opmap.end())
    return nullptr;
  return it->second;
}

BoundaryModifier* BoundaryFactory::findBoundaryMod(const string &s) {
  map<string,BoundaryModifier*>::iterator it;
  it = modmap.find(lowercase(s));
  if(it == modmap.end())
    return nullptr;
  return it->second;
}
