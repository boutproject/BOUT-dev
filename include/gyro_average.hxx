/**************************************************************
 * Gyro-averaging operators
 *
 *
 * 2010-09-03 Ben Dudson <bd512@york.ac.uk>
 *    * Initial version, simple averaging operator
 * 
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 * 
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************/

#ifndef __GYRO_AVERAGE_H__
#define __GYRO_AVERAGE_H__

#include "field3d.hxx"
#include "invert_laplace.hxx"

const int GYRO_FLAGS = INVERT_BNDRY_ONE | INVERT_IN_RHS | INVERT_OUT_RHS;

const Field3D gyroTaylor0(const Field3D &f, const Field3D &rho);

const Field3D gyroPade0(const Field3D &f, const Field3D &rho, 
                        int flags=GYRO_FLAGS);
const Field3D gyroPade1(const Field3D &f, const Field3D &rho, 
                        int flags=GYRO_FLAGS);
const Field3D gyroPade2(const Field3D &f, const Field3D &rho, 
                        int flags=GYRO_FLAGS);

const Field3D gyroPade0(const Field3D &f, const Field2D &rho, 
                        int flags=GYRO_FLAGS);
const Field3D gyroPade1(const Field3D &f, const Field2D &rho, 
                        int flags=GYRO_FLAGS);
const Field3D gyroPade2(const Field3D &f, const Field2D &rho, 
                        int flags=GYRO_FLAGS);

const Field3D gyroPade0(const Field3D &f, BoutReal rho, 
                        int flags=GYRO_FLAGS);
const Field3D gyroPade1(const Field3D &f, BoutReal rho, 
                        int flags=GYRO_FLAGS);
const Field3D gyroPade2(const Field3D &f, BoutReal rho, 
                        int flags=GYRO_FLAGS);

const Field2D gyroPade1(const Field2D &f, const Field2D &rho,
                        int flags=GYRO_FLAGS);

#endif // __GYRO_AVERAGE_H__
