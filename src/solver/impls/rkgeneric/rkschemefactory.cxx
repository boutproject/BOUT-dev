#include "rkschemefactory.hxx"

#include "impls/rkf45/rkf45.hxx"
#include "impls/cashkarp/cashkarp.hxx"
#include "impls/rk4simple/rk4simple.hxx"
#include "impls/rkf34/rkf34.hxx"

#include <boutexception.hxx>

std::string StandardFactoryTraits<RKScheme>::getDefaultType() {
  return RKSCHEME_RKF45;
}
