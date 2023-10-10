/**************************************************************
 * Radial mask operators
 **************************************************************/

#ifndef __MASKX_H__
#define __MASKX_H__

#include "bout/field3d.hxx"

// create a radial buffer zone to set jpar zero near radial boundary
const Field3D mask_x(const Field3D& f, bool BoutRealspace = true);
const Field2D source_tanhx(const Field2D& f, BoutReal swidth, BoutReal slength);
const Field2D source_expx2(const Field2D& f, BoutReal swidth, BoutReal slength);
const Field3D sink_tanhx(const Field2D& f0, const Field3D& f, BoutReal swidth,
                         BoutReal slength, bool BoutRealspace = true);

const Field3D sink_tanhxl(const Field2D& f0, const Field3D& f, BoutReal swidth,
                          BoutReal slength, bool BoutRealspace = true);
const Field3D sink_tanhxr(const Field2D& f0, const Field3D& f, BoutReal swidth,
                          BoutReal slength, bool BoutRealspace = true);

const Field3D buff_x(const Field3D& f, bool BoutRealspace = true);

#endif // __MASKX_H__
