
#include <bout/griddata.hxx>

#include <field_factory.hxx>

#include <bout/constants.hxx>

#include <boutexception.hxx>

#include <output.hxx>

#include <unused.hxx>

bool GridFromOptions::hasVar(const string &name) {
  return options->isSet(name);
}

bool GridFromOptions::get(Mesh *UNUSED(m), int &ival, const string &name) {
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
  std::shared_ptr<FieldGenerator> gen = FieldFactory::get()->parse(expr, options);
  DataIterator i(0,0,0,0,0,0);
  rval = gen->generate(0.0, 0.0, 0.0, 0.0, i ,m);

  return true;
}

bool GridFromOptions::get(Mesh *m, Field2D &var, const string &name, BoutReal def) {
  if (!hasVar(name)) {
    output_warn.write("Variable '%s' not in mesh options. Setting to %e\n", name.c_str(), def);
    var = def;
    return false;
  }

  var = FieldFactory::get()->create2D(name, options, m);
  return true;
}

bool GridFromOptions::get(Mesh *m, Field3D &var, const string &name, BoutReal def) {
  if(!hasVar(name)) {
    var = def;
    return false;
  }

  var = FieldFactory::get()->create3D(name, options, m);
  return true;
}

bool GridFromOptions::get(Mesh *m, vector<int> &var, const string &name, int len,
                          int UNUSED(offset), GridDataSource::Direction UNUSED(dir)) {
  // Integers not expressions yet

  int ival;
  if(!get(m, ival, name))
    return false;

  var.resize(len, ival);

  return true;
}

bool GridFromOptions::get(Mesh *m, vector<BoutReal> &var, const string &name, int len,
                          int offset, GridDataSource::Direction dir) {

  if(!hasVar(name)) {
    return false;
  }

  // Fetch expression as a string
  std::string expr;
  options->get(name, expr, "0");

  // Parse, and evaluate with x,y,z,t = 0
  std::shared_ptr<FieldGenerator> gen = FieldFactory::get()->parse(expr, options);

  var.resize(len);

  switch(dir) {
  case GridDataSource::X: {
    DataIterator i(0,0,0,len,0,0);
    for(int x=0;x<len;x++){
      var[x] = gen->generate(m->GlobalX(x - m->OffsetX + offset), 0.0, 0.0, 0.0,i,m);
    }
    break;
  }
  case GridDataSource::Y : {
    DataIterator i(0,0,0,0,len,0);
    for(int y=0;y<len;y++){
      var[y] = gen->generate(0.0, TWOPI*m->GlobalY(y - m->OffsetY + offset), 0.0, 0.0,i,m);
    }
    break;
  }
  case GridDataSource::Z : {
    DataIterator i(0,0,0,0,0,len);
    for(int z=0;z<len;z++){
      var[z] = gen->generate(0.0, 0.0, TWOPI*(static_cast<BoutReal>(z) + offset) / static_cast<BoutReal>(m->LocalNz), 0.0,i,m);
    }
    break;
  }
  default: {
    throw BoutException("Invalid direction argument");
  }
  }
  return true;
}

