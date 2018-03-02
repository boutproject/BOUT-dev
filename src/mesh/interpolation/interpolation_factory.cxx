#include "interpolation_factory.hxx"
#include "output.hxx"

InterpolationFactory* InterpolationFactory::instance = nullptr;

/**
 * Add the available interpolation methods to the internal map
 *
 */
InterpolationFactory::InterpolationFactory() {
  add(HermiteSpline::CreateHermiteSpline, "hermitespline");
  add(Lagrange4pt::CreateLagrange4pt, "lagrange4pt");
  add(Bilinear::CreateBilinear, "bilinear");
}

inline string InterpolationFactory::getDefaultInterpType() {
  return "hermitespline";
}

/**
 * Create or get the singleton instance of the factory
 *
 *
 * @return The singleton instance of the InterpolationFactory
 */
InterpolationFactory* InterpolationFactory::getInstance() {
  if (instance == nullptr) {
    // Create the singleton object
    instance = new InterpolationFactory();
  }
  return instance;
}

/**
 * Destroy the singleton instance
 *
 */
void InterpolationFactory::cleanup() {
  if (instance == nullptr)
    return;

  // Just delete the instance
  delete instance;
  instance = nullptr;
}

/**
 * Create an Interpolation object given its name and an Options object
 *
 * @param name The name of the interpolation method
 * @param opt An Options object (e.g. an input file)
 *
 * @return A new copy of an Interpolation object
 */
Interpolation* InterpolationFactory::create(Options *options, Mesh *mesh) {
  // Get the default interpolation type
  string type = getDefaultInterpType();

  // If no options section passed (e.g. for a variable), then use the
  // "interpolation" section
  if (options == nullptr)
    options = Options::getRoot()->getSection("interpolation");

  string interp_option;
  options->get("type", interp_option, type);

  if (!interp_option.empty()) type = interp_option.c_str();

  return create(type, options, mesh);
}

Interpolation* InterpolationFactory::create(const string &name, Options *options, Mesh *localmesh) {
  // If no options section passed (e.g. for a variable), then use the
  // "interpolation" section
  if (options == nullptr)
    options = Options::getRoot()->getSection("interpolation");

  // Use the global mesh if none passed
  if (localmesh == nullptr)
    localmesh = mesh;

  auto interp = findInterpolation(name);
  if (interp == nullptr)
    throw BoutException("Could not find interpolation method '%s'", name.c_str());

  return interp(localmesh);
}

void InterpolationFactory::add(CreateInterpCallback interp, const string &name) {
  if ((findInterpolation(name)) != nullptr) {
    // error - already exists
    output << "ERROR: Trying to add an already existing interpolation: " << name << endl;
    return;
  }
  interp_map[lowercase(name)] = interp;
}

/**
 * Find an interpolation method in the list of available methods
 *
 * @param name Name of the interpolation method
 *
 * @return A pointer to the Interpolation object in the map
 */
InterpolationFactory::CreateInterpCallback InterpolationFactory::findInterpolation(const string &name) {
  auto interp = interp_map.find(lowercase(name));
  if (interp == end(interp_map))
    return nullptr;
  return interp->second;
}
