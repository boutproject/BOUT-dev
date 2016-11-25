/***************************************************************************
 * Define some useful constants
 * 
 **************************************************************************/

#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <bout_types.hxx>

/// Mathematical constant pi
const BoutReal PI = 3.141592653589793;
const BoutReal TWOPI = 6.2831853071795;

namespace SI {
  // Constants in SI system of units
  
  const BoutReal e0  = 8.854e-12;      ///< Permittivity of free space
  const BoutReal mu0 = 4.e-7*PI;       ///< Permeability of free space
  const BoutReal qe  = 1.60217646e-19;      ///< Electron charge
  const BoutReal Me  = 9.1093816e-31;      ///< Electron mass
  const BoutReal Mp  = 1.67262158e-27; ///< Proton mass
}

#endif // __CONSTANTS_H__
