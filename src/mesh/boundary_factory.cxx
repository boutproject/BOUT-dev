#include "bout/parallel_boundary_op.hxx"
#include "bout/parallel_boundary_region.hxx"
#include <bout/boundary_factory.hxx>
#include <bout/boundary_standard.hxx>
#include <bout/globals.hxx>
#include <bout/options.hxx>
#include <bout/utils.hxx>

#include <list>
#include <map>
#include <string>
using std::list;
using std::string;

#include <bout/output.hxx>

BoundaryFactory* BoundaryFactory::instance = nullptr;

BoundaryFactory::BoundaryFactory() {
  add(new BoundaryDirichlet(), "dirichlet");
  add(new BoundaryDirichlet(), "dirichlet_o2"); // Synonym for "dirichlet"
  add(new BoundaryDirichlet_O3(), "dirichlet_o3");
  add(new BoundaryDirichlet_O4(), "dirichlet_o4");
  add(new BoundaryDirichlet_4thOrder(), "dirichlet_4thorder");
  add(new BoundaryNeumann(), "neumann");
  add(new BoundaryNeumann(), "neumann_O2"); // Synonym for "neumann"
  add(new BoundaryNeumann_4thOrder(), "neumann_4thorder");
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

  addMod(new BoundaryRelax(), "relax");
  addMod(new BoundaryWidth(), "width");
  addMod(new BoundaryToFieldAligned(), "toFieldAligned");
  addMod(new BoundaryFromFieldAligned(), "fromFieldAligned");

  // Parallel boundaries
  add(new BoundaryOpPar_dirichlet_o1(), "parallel_dirichlet_o1");
  add(new BoundaryOpPar_dirichlet_o2(), "parallel_dirichlet_o2");
  add(new BoundaryOpPar_dirichlet_o3(), "parallel_dirichlet_o3");
  add(new BoundaryOpPar_neumann_o1(), "parallel_neumann_o1");
  add(new BoundaryOpPar_neumann_o2(), "parallel_neumann_o2");
  add(new BoundaryOpPar_neumann_o3(), "parallel_neumann_o3");
}

BoundaryFactory::~BoundaryFactory() {
  // Free any boundaries
  for (const auto& it : opmap) {
    delete it.second;
  }
  for (const auto& it : modmap) {
    delete it.second;
  }
  for (const auto& it : par_opmap) {
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
  if (instance == nullptr) {
    return;
  }

  // Just delete the instance
  delete instance;
  instance = nullptr;
}

BoundaryOpBase* BoundaryFactory::create(const string& name, BoundaryRegionBase* region) {

  // Search for a string of the form: modifier(operation)
  auto pos = name.find('(');
  if (pos == string::npos) {
    // No more (opening) brackets. Should be a boundary operation
    // Need to strip whitespace

    if ((name == "null") || (name == "none")) {
      return nullptr;
    }

    if (region->isParallel) {
      // Parallel boundary
      BoundaryOpPar* pop = findBoundaryOpPar(trim(name));
      if (pop == nullptr) {
        throw BoutException("Could not find parallel boundary condition '{:s}'", name);
      }

      // Clone the boundary operation, passing the region to operate over,
      // an empty args list and empty keyword map
      list<string> args;
      return pop->clone(dynamic_cast<BoundaryRegionPar*>(region), args, {});
    } else {
      // Perpendicular boundary
      BoundaryOp* op = findBoundaryOp(trim(name));
      if (op == nullptr) {
        throw BoutException("Could not find boundary condition '{:s}'", name);
      }

      // Clone the boundary operation, passing the region to operate over,
      // an empty args list and empty keyword map
      list<string> args;
      return op->clone(dynamic_cast<BoundaryRegion*>(region), args, {});
    }
  }
  // Contains a bracket. Find the last bracket and remove
  auto pos2 = name.rfind(')');
  if (pos2 == string::npos) {
    output_warn << "\tWARNING: Unmatched brackets in boundary condition: " << name
                << endl;
  }

  // Find the function name before the bracket
  string func = trim(name.substr(0, pos));
  // And the argument inside the bracket
  string arg = trim(name.substr(pos + 1, pos2 - pos - 1));
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
        string s = arg.substr(start, i - start);

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
  std::string s = arg.substr(start);
  auto poseq = s.find('=');
  if (poseq != string::npos) {
    keywords[trim(s.substr(0, poseq))] = trim(s.substr(poseq + 1));
  } else {
    // No '=', so a positional argument
    arglist.push_back(trim(s));
  }

  // Test if func is a modifier
  BoundaryModifier* mod = findBoundaryMod(func);
  if (mod != nullptr) {
    // The first argument should be an operation
    auto* op = dynamic_cast<BoundaryOp*>(create(arglist.front(), region));
    if (op == nullptr) {
      return nullptr;
    }

    // Remove the first element (name of operation)
    arglist.pop_front();

    // Clone the modifier, passing in the operator and remaining strings as argument
    return mod->cloneMod(op, arglist);
  }

  if (region->isParallel) {
    // Parallel boundary
    BoundaryOpPar* pop = findBoundaryOpPar(trim(func));
    if (pop != nullptr) {
      // An operation with arguments
      return pop->clone(dynamic_cast<BoundaryRegionPar*>(region), arglist, keywords);
    }
  } else {
    // Perpendicular boundary
    BoundaryOp* op = findBoundaryOp(trim(func));
    if (op != nullptr) {
      // An operation with arguments
      return op->clone(dynamic_cast<BoundaryRegion*>(region), arglist, keywords);
    }
  }

  // Otherwise nothing matches
  throw BoutException("  Boundary setting is neither an operation nor modifier: {:s}\n",
                      func);

  return nullptr;
}

BoundaryOpBase* BoundaryFactory::create(const char* name, BoundaryRegionBase* region) {
  return create(string(name), region);
}

BoundaryOpBase* BoundaryFactory::createFromOptions(const string& varname,
                                                   BoundaryRegionBase* region) {
  if (region == nullptr) {
    return nullptr;
  }

  output_info << "\t" << region->label << " region: ";

  string prefix("bndry_");

  string side;
  switch (region->location) {
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
  case BNDRY_PAR_FWD_XIN: {
    side = "par_yup_xin";
    break;
  }
  case BNDRY_PAR_FWD_XOUT: {
    side = "par_yup_xout";
    break;
  }
  case BNDRY_PAR_BKWD_XIN: {
    side = "par_ydown_xin";
    break;
  }
  case BNDRY_PAR_BKWD_XOUT: {
    side = "par_ydown_xout";
    break;
  }
  default: {
    side = "all";
    break;
  }
  }

  // Get options
  Options* options = Options::getRoot();

  // Get variable options
  Options* varOpts = options->getSection(varname);
  string set;

  /// First try looking for (var, region)
  if (varOpts->isSet(prefix + region->label)) {
    varOpts->get(prefix + region->label, set, "");
    return create(set, region);
  }

  /// Then (var, side)
  if (varOpts->isSet(prefix + side)) {
    varOpts->get(prefix + side, set, "");
    return create(set, region);
  }

  /// Then (var, all)
  if (region->isParallel) {
    if (varOpts->isSet(prefix + "par_all")) {
      varOpts->get(prefix + "par_all", set, "");
      return create(set, region);
    }
  } else {
    if (varOpts->isSet(prefix + "all")) {
      varOpts->get(prefix + "all", set, "");
      return create(set, region);
    }
  }

  // Get the "all" options
  varOpts = options->getSection("all");

  /// Then (all, region)
  if (varOpts->isSet(prefix + region->label)) {
    varOpts->get(prefix + region->label, set, "");
    return create(set, region);
  }

  /// Then (all, side)
  if (varOpts->isSet(prefix + side)) {
    varOpts->get(prefix + side, set, "");
    return create(set, region);
  }

  /// Then (all, all)
  if (region->isParallel) {
    // Different default for parallel boundary regions
    varOpts->get(prefix + "par_all", set, "parallel_dirichlet_o2");
  } else {
    varOpts->get(prefix + "all", set, "dirichlet");
  }
  return create(set, region);
  // Defaults to Dirichlet conditions, to prevent undefined boundary
  // values. If a user want to override, specify "none" or "null"
}

BoundaryOpBase* BoundaryFactory::createFromOptions(const char* varname,
                                                   BoundaryRegionBase* region) {
  return createFromOptions(string(varname), region);
}

void BoundaryFactory::add(BoundaryOp* bop, const string& name) {
  if ((findBoundaryMod(name) != nullptr) || (findBoundaryOp(name) != nullptr)) {
    // error - already exists
    output_error << "ERROR: Trying to add an already existing boundary: " << name << endl;
    return;
  }
  opmap[lowercase(name)] = bop;
}

void BoundaryFactory::add(BoundaryOp* bop, const char* name) { add(bop, string(name)); }

void BoundaryFactory::add(BoundaryOpPar* bop, const string& name) {
  if (findBoundaryOpPar(name) != nullptr) {
    // error - already exists
    output_error << "ERROR: Trying to add an already existing boundary: " << name << endl;
    return;
  }
  par_opmap[lowercase(name)] = bop;
}

void BoundaryFactory::add(BoundaryOpPar* bop, const char* name) {
  add(bop, string(name));
}

void BoundaryFactory::addMod(BoundaryModifier* bmod, const string& name) {
  if ((findBoundaryMod(name) != nullptr) || (findBoundaryOp(name) != nullptr)) {
    // error - already exists
    output_error << "ERROR: Trying to add an already existing boundary modifier: " << name
                 << endl;
    return;
  }
  modmap[lowercase(name)] = bmod;
}

void BoundaryFactory::addMod(BoundaryModifier* bmod, const char* name) {
  addMod(bmod, string(name));
}

BoundaryOp* BoundaryFactory::findBoundaryOp(const string& s) {
  auto it = opmap.find(lowercase(s));
  if (it == opmap.end()) {
    return nullptr;
  }
  return it->second;
}

BoundaryModifier* BoundaryFactory::findBoundaryMod(const string& s) {
  auto it = modmap.find(lowercase(s));
  if (it == modmap.end()) {
    return nullptr;
  }
  return it->second;
}

BoundaryOpPar* BoundaryFactory::findBoundaryOpPar(const string& s) {
  auto it = par_opmap.find(lowercase(s));
  if (it == par_opmap.end()) {
    return nullptr;
  }
  return it->second;
}
