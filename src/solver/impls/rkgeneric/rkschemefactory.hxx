#ifndef __RKSCHEME_FACTORY_H__
#define __RKSCHEME_FACTORY_H__

#include <bout/rkscheme.hxx>
#include "bout/generic_factory.hxx"

template<>
struct StandardFactoryTraits<RKScheme> {
  static constexpr auto type_name = "RKScheme";
  static constexpr auto section_name = "solver";
  static constexpr auto option_name = "scheme";
  static std::string getDefaultType();
};

class RKSchemeFactory : public StandardFactory<RKScheme, RKSchemeFactory>{};

/// Simpler name for Factory registration helper class
///
/// Usage:
///
///     #include <bout/rkschemefactory.hxx>
///     namespace {
///     RegisterRKScheme<MyRKScheme> registerrkschememine("myrkscheme");
///     }
template <typename DerivedType>
using RegisterRKScheme = RegisterInStandardFactory<RKScheme, DerivedType, RKSchemeFactory>;

#endif // __RKSCHEME_FACTORY_H__

