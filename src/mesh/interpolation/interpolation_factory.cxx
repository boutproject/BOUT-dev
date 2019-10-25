#include "interpolation_factory.hxx"
#include "output.hxx"

InterpolationFactory* InterpolationFactory::instance = nullptr;

InterpolationFactory::InterpolationFactory() {
  add(HermiteSpline::CreateHermiteSpline, "hermitespline");
  add(MonotonicHermiteSpline::CreateMonotonicHermiteSpline, "monotonichermitespline");
  add(Lagrange4pt::CreateLagrange4pt, "lagrange4pt");
  add(Bilinear::CreateBilinear, "bilinear");
}

InterpolationFactory* InterpolationFactory::getInstance() {
  if (instance == nullptr) {
    // Create the singleton object
    instance = new InterpolationFactory();
  }
  return instance;
}

void InterpolationFactory::cleanup() {
  if (instance == nullptr)
    return;

  // Just delete the instance
  delete instance;
  instance = nullptr;
}

Interpolation* InterpolationFactory::create(Options *options, Mesh *mesh) {
  // Get the default interpolation type
  std::string type = getDefaultInterpType();

  // If no options section passed (e.g. for a variable), then use the
  // "interpolation" section
  if (options == nullptr)
    options = Options::getRoot()->getSection("interpolation");

  std::string interp_option = (*options)["type"].withDefault(type);

  if (!interp_option.empty()) type = interp_option.c_str();

  return create(type, options, mesh);
}

Interpolation* InterpolationFactory::create(const std::string &name, Options *options, Mesh *localmesh) {
  // If no options section passed (e.g. for a variable), then use the
  // "interpolation" section
  if (options == nullptr) {
    options = Options::getRoot()->getSection("interpolation");
  }

  // Use the global mesh if none passed
  if (localmesh == nullptr) {
    localmesh = bout::globals::mesh;
  }

  auto interp = findInterpolation(name);
  if (interp == nullptr) {
    throw BoutException("Could not find interpolation method '%s'", name.c_str());
  }

  return interp(localmesh);
}

void InterpolationFactory::add(CreateInterpCallback interp, const std::string &name) {
  if (findInterpolation(name) != nullptr) {
    output_warn << "ERROR: Trying to add an already existing interpolation: " << name << endl;
    return;
  }
  interp_map[lowercase(name)] = interp;
}

InterpolationFactory::CreateInterpCallback InterpolationFactory::findInterpolation(const std::string &name) {
  auto interp = interp_map.find(lowercase(name));
  if (interp == end(interp_map))
    return nullptr;
  return interp->second;
}
