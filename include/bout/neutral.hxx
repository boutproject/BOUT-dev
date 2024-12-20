/**************************************************************************
 * Neutral reaction rates based on Reiter's polynomial fitting formula 
 * (http://www.eirene.de/amjuel.pdf)
 *
 **************************************************************************
 * by B. Zhu (02/05/2020)
 *
 **************************************************************************/

#ifndef __NEUTRAL_H__
#define __NEUTRAL_H__

#include "bout/field3d.hxx"
#include "bout/field2d.hxx"

// iz: ionization; cx: charge exchange; rc: recombination; sec: cross section
const Field3D iz_rate(const Field3D &tin, BoutReal AA);
const Field3D cx_rate(const Field3D &tin, BoutReal AA);
const Field3D rc_rate(const Field3D &tin, const Field3D &nin, BoutReal AA);
const Field3D cx_sect(const Field3D &tin, BoutReal AA);

// 2D functions have not been implemented yet
const Field2D iz_rate(const Field2D &tin, BoutReal AA);
const Field2D cx_rate(const Field2D &tin, BoutReal AA);
const Field2D rc_rate(const Field2D &tin, Field2D &nin, BoutReal AA);
const Field2D cx_sect(const Field2D &tin, BoutReal AA);


#endif // __NEUTRAL_H__

