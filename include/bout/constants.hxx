/***************************************************************************
 * Define some useful constants
 * 
 **************************************************************************/

#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <bout_types.hxx>

/// Mathematical constant pi
const BoutReal PI = 3.141592653589793;
const BoutReal TWOPI = 2 * PI;

namespace SI {
  // Constants in SI system of units
  const BoutReal c   = 299792458;       ///< Speed of light in vacuum
  const BoutReal mu0 = 4.e-7*PI;        ///< Permeability of free space
  const BoutReal e0  = 1/(c*c*mu0);     ///< Permittivity of free space
  const BoutReal qe  = 1.602176634e-19; ///< Electron charge
  const BoutReal Me  = 9.10938356e-31;  ///< Electron mass
  const BoutReal Mp  = 1.672621898e-27; ///< Proton mass
  const BoutReal kb  = 1.38064852e-23;  ///< Boltzmanns constant
}

#endif // __CONSTANTS_H__
