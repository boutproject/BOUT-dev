
#include <bout/griddata.hxx>

#include <field_factory.hxx>

#include <bout/constants.hxx>

#include <boutexception.hxx>

#include <output.hxx>

bool GridFromOptions::hasVar(const string &name) {
  return options->isSet(name);
}

/*!
 * Reads integers from options. Assumes that integers are not
 * expressions.
 * 
 * Inputs
 * ======
 *
 * m      [Mesh pointer] Not used
 * name   [string] containing name of variable
 * 
 * Outputs
 * =======
 * 
 * ival    [integer] Always given a value, defaults to 0
 *
 * Returns
 * =======
 *
 * True if option is set, false if ival is default (0)
 */
bool GridFromOptions::get(Mesh *m, int &ival,      const string &name) {
  if(!hasVar(name)) {
    ival = 0;
    return false;
  }
  
  options->get(name, ival, 0);
  return true;
}

bool GridFromOptions::get(Mesh *m, BoutReal &rval, const string &name) {
  if(!hasVar(name)) {
    rval = 0.0;
    return false;
  }
  
  // Fetch expression as a string 
  std::string expr;
  options->get(name, expr, "0");

  // Parse, and evaluate with x,y,z,t = 0
  FieldGenerator* gen = FieldFactory::get()->parse(expr, options);
  rval = gen->generate(0.0, 0.0, 0.0, 0.0);

  return true;
}

bool GridFromOptions::get(Mesh *m, Field2D &var,   const string &name, BoutReal def) {
  if(!hasVar(name)) {
    output.write("Variable '%s' not in mesh options. Setting to %e\n", name.c_str(), def);
    var = def;
    return false;
  }

  var = FieldFactory::get()->create2D(name, options, m);
  return true;
}

bool GridFromOptions::get(Mesh *m, Field3D &var,   const string &name, BoutReal def) {
  if(!hasVar(name)) {
    var = def;
    return false;
  }

  var = FieldFactory::get()->create3D(name, options, m);
  return true;
}

bool GridFromOptions::get(Mesh *m, vector<int> &var,      const string &name, int len, int offset, GridDataSource::Direction dir) {
  // Integers not expressions yet

  int ival;
  if(!get(m, ival, name))
    return false;

  var.resize(len);

  for(int i=0; i<var.size(); i++)
    var[i] = ival;

  return true;
}

bool GridFromOptions::get(Mesh *m, vector<BoutReal> &var, const string &name, int len, int offset, GridDataSource::Direction dir) {
  
  if(!hasVar(name)) {
    return false;
  }
  
  // Fetch expression as a string 
  std::string expr;
  options->get(name, expr, "0");
  
  // Parse, and evaluate with x,y,z,t = 0
  FieldGenerator* gen = FieldFactory::get()->parse(expr, options);

  var.resize(len);
  
  switch(dir) {
  case GridDataSource::X: {
    for(int x=0;x<var.size();x++){
      var[x] = gen->generate(m->GlobalX(x - m->OffsetX + offset), 0.0, 0.0, 0.0);
    }
    break;
  }
  case GridDataSource::Y : {
    for(int y=0;y<var.size();y++){
      var[y] = gen->generate(0.0, TWOPI*m->GlobalY(y - m->OffsetY + offset), 0.0, 0.0);
    }
  }
  case GridDataSource::Z : {
    for(int z=0;z<var.size();z++){
      var[z] = gen->generate(0.0, 0.0, TWOPI*((BoutReal) z + offset) / ((BoutReal) (m->ngz-1)), 0.0);
    }
  }
  default: {
    throw BoutException("Invalid direction argument");
  }
  }
  return true;
}

