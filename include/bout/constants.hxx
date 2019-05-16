/***************************************************************************
 * Define some useful constants
 * 
 **************************************************************************/

#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <bout_types.hxx>

/// Mathematical constant pi
constexpr BoutReal PI = 3.141592653589793;
/// Mathematical constant 2 * pi
constexpr BoutReal TWOPI = 2 * PI;

namespace SI {
  // Constants in SI system of units
  constexpr BoutReal c   = 299792458;       ///< Speed of light in vacuum
  constexpr BoutReal mu0 = 4.e-7*PI;        ///< Permeability of free space
  constexpr BoutReal e0  = 1/(c*c*mu0);     ///< Permittivity of free space
  constexpr BoutReal qe  = 1.602176634e-19; ///< Electron charge
  constexpr BoutReal Me  = 9.10938356e-31;  ///< Electron mass
  constexpr BoutReal Mp  = 1.672621898e-27; ///< Proton mass
  constexpr BoutReal kb  = 1.38064852e-23;  ///< Boltzmanns constant
  constexpr BoutReal amu = 1.660539040e-27;          ///< Unified atomic mass unit
  constexpr BoutReal M_Hydrogen = 1.008 * amu;       ///< Mass of a Hydrogen atom
  constexpr BoutReal M_Deuterium = 2.01410178 * amu; ///< Mass of a Deuterium atom
  constexpr BoutReal M_Tritium = 3.0160492 * amu;    ///< Mass of a Tritium atom
}

#endif // __CONSTANTS_H__
