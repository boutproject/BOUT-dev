
class LaplaceFactory;

#ifndef __LAPLACE_FACTORY_H__
#define __LAPLACE_FACTORY_H__

#include <invert_laplace.hxx>
#include "bout/generic_factory.hxx"

template<>
struct StandardFactoryTraits<Laplacian> {
  static constexpr auto type_name = "Laplacian";
  static constexpr auto section_name = "laplace";
  static constexpr auto option_name = "type";
  static std::string getDefaultType();
};

class LaplaceFactory
    : public StandardFactory<
          Laplacian, LaplaceFactory,
          std::function<std::unique_ptr<Laplacian>(Options*, CELL_LOC, Mesh*)>> {
public:
  ReturnType create(Options* options = nullptr, CELL_LOC loc = CELL_CENTRE,
                    Mesh* mesh = nullptr) {
    return StandardFactory::create(getType(options), options, loc, mesh);
  }
};

/// Simpler name for Factory registration helper class
///
/// Usage:
///
///     #include <bout/laplacefactory.hxx>
///     namespace {
///     RegisterLaplace<MyLaplace> registerlaplacemine("mylaplace");
///     }
template <class DerivedType>
class RegisterInStandardFactory<Laplacian, DerivedType, LaplaceFactory> {
public:
  RegisterInStandardFactory(const std::string& name) {
    LaplaceFactory::getInstance().add(
      name, [](Options* options, CELL_LOC loc, Mesh* mesh) -> std::unique_ptr<Laplacian> {
        return std::make_unique<DerivedType>(options, loc, mesh);
      });
  }
};

template <class DerivedType>
using RegisterLaplace = RegisterInStandardFactory<Laplacian, DerivedType, LaplaceFactory>;

#endif // __LAPLACE_FACTORY_H__
