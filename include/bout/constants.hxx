/***************************************************************************
 * Define some useful constants
 * 
 **************************************************************************/

#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <bout_types.hxx>

/// Mathematical constant pi
const BoutReal PI = 3.141592653589793;
/// Mathematical constant 2 * pi
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
  const BoutReal amu = 1.660539040e-27;          ///< Unified atomic mass unit
  const BoutReal M_Hydrogen = 1.008 * amu;       ///< Mass of a Hydrogen atom
  const BoutReal M_Deuterium = 2.01410178 * amu; ///< Mass of a Deuterium atom
  const BoutReal M_Tritium = 3.0160492 * amu;    ///< Mass of a Tritium atom
}

#endif // __CONSTANTS_H__
