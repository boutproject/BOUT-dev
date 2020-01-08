#include <bout/constants.hxx>
#include <bout/griddata.hxx>
#include <boutexception.hxx>
#include <field_factory.hxx>
#include <output.hxx>
#include <unused.hxx>

bool GridFromOptions::hasVar(const std::string& name) { return options->isSet(name); }

namespace {
/// Return value of \p name in \p options, using \p def as a default
/// if it doesn't exist
template <class T>
auto getWithDefault(const Options& options, const std::string& name, const T& def) -> T {
  const bool has_var = options.isSet(name);
  if (!has_var) {
    output_warn.write("Variable '{:s}' not in mesh options. Setting to ", name);
    output_warn << def << "\n";
  }
  // Note! We don't use `Options::withDefault` here because that
  // records the default in the `Options` object and we don't know if
  // we'll actually end up using that value. This is because
  // `GridFromOptions::get` is probably being called via `Mesh::get`,
  // and the calling site may use the return value of that to set its
  // own default
  return has_var ? options[name].as<T>() : def;
}
} // namespace

bool GridFromOptions::get(Mesh*, std::string& sval, const std::string& name,
                          const std::string& def) {
  sval = getWithDefault(*options, name, def);
  return hasVar(name);
}

bool GridFromOptions::get(Mesh*, int& ival, const std::string& name, int def) {
  ival = getWithDefault(*options, name, def);
  return hasVar(name);
}

bool GridFromOptions::get(Mesh*, BoutReal& rval, const std::string& name, BoutReal def) {
  rval = getWithDefault(*options, name, def);
  return hasVar(name);
}

bool GridFromOptions::get(Mesh* m, Field2D& var, const std::string& name, BoutReal def) {
  if (!hasVar(name)) {
    output_warn.write("Variable '{:s}' not in mesh options. Setting to {:e}\n", name,
                      def);
    var = def;
    return false;
  }

  var = FieldFactory::get()->create2D(name, options, m);
  return true;
}

bool GridFromOptions::get(Mesh* m, Field3D& var, const std::string& name, BoutReal def) {
  if (!hasVar(name)) {
    output_warn.write("Variable '{:s}' not in mesh options. Setting to {:e}\n", name,
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
    output_warn.write("Variable '{:s}' not in mesh options. Setting to {:e}\n", name,
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
    output_warn.write("Variable '{:s}' not in mesh options. Setting to empty vector\n", name);
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
    output_warn.write("Variable '{:s}' not in mesh options. Setting to empty vector\n", name);
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
