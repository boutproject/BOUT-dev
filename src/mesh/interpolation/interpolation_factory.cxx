#include "interpolation_factory.hxx"

InterpolationFactory* InterpolationFactory::instance = nullptr;

InterpolationFactory::InterpolationFactory() {
  add(new HermiteSpline(), "hermitespline");
}

InterpolationFactory::~InterpolationFactory() {
  // Free any interpolation objects
  for(auto& interp : interp_map) {
    delete interp->second;
  }
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
  instance == nullptr;
}

Interpolation* InterpolationFactory::create(const string &name, Options *opt) {
  
}

void InterpolationFactory::add(Interpolation* interp, const string &name) {
  if ((findInterpolation(name)) != nullptr) {
    // error - already exists
    output << "ERROR: Trying to add an already existing interpolation: " << name << endl;
    return;
  }
  interp_map[lowercase(name)] = interp;
}

Interpolation* InterpolationFactory::findInterpolation(const string &name) {
  auto interp = interp_map.find(lowercase(name));
  if (interp == interp_map.end())
    return nullptr;
  return interp->second;
}
    
