/*!
 * \file mesh_facefield_comm.hxx
 *
 * \brief Extension to Mesh class for face-centered field communication
 *
 * \author BOUT++ Team
 *
 **************************************************************************
 * Copyright 2024 BOUT++ Contributors
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
 */

#pragma once
#ifndef BOUT_MESH_FACEFIELD_COMM_H
#define BOUT_MESH_FACEFIELD_COMM_H

#include "bout/facefield3d.hxx"
#include "bout/mesh.hxx"

/*!
 * \brief Extension methods for Mesh to handle FaceField3D communication
 * 
 * This file provides convenience methods for communicating face-centered
 * fields across processor boundaries. The communication handles the
 * staggered grid locations appropriately.
 */

/*!
 * \brief Communicate all components of a FaceField3D
 * 
 * This method communicates the x, y, and z components of the face field,
 * handling the staggered locations correctly. Each component is communicated
 * with its appropriate CELL_*LOW location.
 * 
 * \param mesh The mesh to use for communication
 * \param field The FaceField3D to communicate
 */
inline void communicate(Mesh* mesh, FaceField3D& field) {
  // Communicate each component separately
  // Each component already has the correct staggered location
  mesh->communicate(field.x(), field.y(), field.z());
}

/*!
 * \brief Communicate only XZ components of a FaceField3D
 * 
 * \param mesh The mesh to use for communication
 * \param field The FaceField3D to communicate
 */
inline void communicateXZ(Mesh* mesh, FaceField3D& field) {
  // Only communicate X and Z components
  mesh->communicateXZ(field.x(), field.z());
}

/*!
 * \brief Communicate only YZ components of a FaceField3D
 * 
 * \param mesh The mesh to use for communication
 * \param field The FaceField3D to communicate
 */
inline void communicateYZ(Mesh* mesh, FaceField3D& field) {
  // Only communicate Y and Z components
  mesh->communicateYZ(field.y(), field.z());
}

// Convenience methods that use the mesh from the field
inline void communicate(FaceField3D& field) {
  communicate(field.getMesh(), field);
}

inline void communicateXZ(FaceField3D& field) {
  communicateXZ(field.getMesh(), field);
}

inline void communicateYZ(FaceField3D& field) {
  communicateYZ(field.getMesh(), field);
}

#endif // BOUT_MESH_FACEFIELD_COMM_H