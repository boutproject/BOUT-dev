#include <bout/constants.hxx>
#include <bout/griddata.hxx>
#include <boutexception.hxx>
#include <field_factory.hxx>
#include <output.hxx>
#include <unused.hxx>

bool GridFromOptions::hasVar(const std::string& name) { return options->isSet(name); }

bool GridFromOptions::get(Mesh* UNUSED(m), std::string& sval, const std::string& name) {
  if (!hasVar(name)) {
    const std::string def{};
    output_warn.write("Variable '%s' not in mesh options. Setting to \"%s\"\n",
                      name.c_str(), def.c_str());
    sval = def;
    return false;
  }

  options->get(name, sval, "");
  return true;
}

bool GridFromOptions::get(Mesh* UNUSED(m), int& ival, const std::string& name, int def) {
  if (!hasVar(name)) {
    output_warn.write("Variable '%s' not in mesh options. Setting to %d\n", name.c_str(),
                      def);
    ival = def;
    return false;
  }

  options->get(name, ival, 0);
  return true;
}

bool GridFromOptions::get(Mesh* UNUSED(m), BoutReal& rval, const std::string& name) {
  if (!hasVar(name)) {
    constexpr BoutReal def{0.0};
    output_warn.write("Variable '%s' not in mesh options. Setting to %e\n", name.c_str(),
                      def);
    rval = def;
    return false;
  }

  rval = (*options)[name].withDefault(0.0);

  return true;
}

bool GridFromOptions::get(Mesh* m, Field2D& var, const std::string& name, BoutReal def) {
  if (!hasVar(name)) {
    output_warn.write("Variable '%s' not in mesh options. Setting to %e\n", name.c_str(),
                      def);
    var = def;
    return false;
  }

  var = FieldFactory::get()->create2D(name, options, m);
  return true;
}

bool GridFromOptions::get(Mesh* m, Field3D& var, const std::string& name, BoutReal def) {
  if (!hasVar(name)) {
    output_warn.write("Variable '%s' not in mesh options. Setting to %e\n", name.c_str(),
                      def);
    var = def;
    return false;
  }

  var = FieldFactory::get()->create3D(name, options, m);
  return true;
}

bool GridFromOptions::get(Mesh* m, FieldPerp& var, const std::string& name, BoutReal def) {
  // Cannot set attributes from options at the moment, so don't know what 'yindex' this
  // FieldPerp should have: just set to 0 for now, and create FieldPerp on all processors
  // (note: this is different to behaviour of GridFromFile which will only create the
  // FieldPerp at a single global y-index).

  if (!hasVar(name)) {
    output_warn.write("Variable '%s' not in mesh options. Setting to %e\n", name.c_str(),
                      def);
    var = def;
    var.setIndex(0);
    return false;
  }

  var = FieldFactory::get()->createPerp(name, options, m);
  var.setIndex(0);

  return true;
}

bool GridFromOptions::get(Mesh* m, std::vector<int>& var, const std::string& name,
                          int len, int UNUSED(offset),
                          GridDataSource::Direction UNUSED(dir)) {
  if (!hasVar(name)) {
    std::vector<int> def{};
    output_warn.write("Variable '%s' not in mesh options. Setting to empty vector\n",
                      name.c_str());
    var = def;
    return false;
  }

  // FIXME: actually implement this!
  throw BoutException("not implemented");
  int ival;
  get(m, ival, name);
  var.resize(len, ival);

  return true;
}

bool GridFromOptions::get(Mesh* m, std::vector<BoutReal>& var, const std::string& name,
                          int len, int offset, GridDataSource::Direction dir) {
  if (!hasVar(name)) {
    std::vector<BoutReal> def{};
    output_warn.write("Variable '%s' not in mesh options. Setting to empty vector\n",
                      name.c_str());
    var = def;
    return false;
  }

  auto expr = (*options)[name].withDefault(std::string{"0"});
  auto gen = FieldFactory::get()->parse(expr, options);

  var.resize(len);

  switch (dir) {
  case GridDataSource::X: {
    for (int x = 0; x < len; x++) {
      var[x] = gen->generate(m->GlobalX(x - m->OffsetX + offset), 0.0, 0.0, 0.0);
    }
    break;
  }
  case GridDataSource::Y: {
    for (int y = 0; y < len; y++) {
      var[y] = gen->generate(0.0, TWOPI * m->GlobalY(y - m->OffsetY + offset), 0.0, 0.0);
    }
    break;
  }
  case GridDataSource::Z: {
    for (int z = 0; z < len; z++) {
      var[z] = gen->generate(
          0.0, 0.0,
          (TWOPI * (z - m->OffsetZ + offset)) / static_cast<BoutReal>(m->LocalNz), 0.0);
    }
    break;
  }
  default: { throw BoutException("Invalid direction argument"); }
  }
  return true;
}
