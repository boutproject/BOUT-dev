/********************************************************
 * BOUT++ Library - Write fluid simulations in curviilinear geometry
 * Copyright (C) 2016, 2017, 2018 David Schwörer
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
 *****************************************************************/

#include "aiolos_init.hxx"
#include "aiolosmesh.hxx"

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexDDX_norm(const Field3D &f, CELL_LOC outloc,
                                        DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_First][0];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_x_norm(f);
    break;
  case DIFF_W2:
    return DDX_CWENO2_x_norm(f);
    break;
  case DIFF_C4:
    return DDX_C4_x_norm(f);
    break;
  case DIFF_S2:
    return DDX_S2_x_norm(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDX_norm unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_W2"
                        " * DIFF_C4"
                        " * DIFF_S2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexDDY_norm(const Field3D &f, CELL_LOC outloc,
                                        DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_First][1];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_y_norm(f);
    break;
  case DIFF_W2:
    return DDX_CWENO2_y_norm(f);
    break;
  case DIFF_C4:
    return DDX_C4_y_norm(f);
    break;
  case DIFF_S2:
    return DDX_S2_y_norm(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDY_norm unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_W2"
                        " * DIFF_C4"
                        " * DIFF_S2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexDDZ_norm(const Field3D &f, CELL_LOC outloc,
                                        DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_First][2];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_z_norm(f);
    break;
  case DIFF_W2:
    return DDX_CWENO2_z_norm(f);
    break;
  case DIFF_C4:
    return DDX_C4_z_norm(f);
    break;
  case DIFF_S2:
    return DDX_S2_z_norm(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDZ_norm unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_W2"
                        " * DIFF_C4"
                        " * DIFF_S2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexDDX_norm(const Field2D &f, CELL_LOC outloc,
                                        DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_First][0];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_x_norm(f);
    break;
  case DIFF_W2:
    return DDX_CWENO2_x_norm(f);
    break;
  case DIFF_C4:
    return DDX_C4_x_norm(f);
    break;
  case DIFF_S2:
    return DDX_S2_x_norm(f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexDDX_norm unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_W2"
                        " * DIFF_C4"
                        " * DIFF_S2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexDDY_norm(const Field2D &f, CELL_LOC outloc,
                                        DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_First][1];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_y_norm(f);
    break;
  case DIFF_W2:
    return DDX_CWENO2_y_norm(f);
    break;
  case DIFF_C4:
    return DDX_C4_y_norm(f);
    break;
  case DIFF_S2:
    return DDX_S2_y_norm(f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexDDY_norm unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_W2"
                        " * DIFF_C4"
                        " * DIFF_S2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexD2DX2_norm(const Field3D &f, CELL_LOC outloc,
                                          DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_Second][0];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_x_norm(f);
    break;
  case DIFF_C4:
    return D2DX2_C4_x_norm(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DX2_norm unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexD2DY2_norm(const Field3D &f, CELL_LOC outloc,
                                          DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_Second][1];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_y_norm(f);
    break;
  case DIFF_C4:
    return D2DX2_C4_y_norm(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DY2_norm unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexD2DZ2_norm(const Field3D &f, CELL_LOC outloc,
                                          DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_Second][2];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_z_norm(f);
    break;
  case DIFF_C4:
    return D2DX2_C4_z_norm(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DZ2_norm unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexD2DX2_norm(const Field2D &f, CELL_LOC outloc,
                                          DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_Second][0];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_x_norm(f);
    break;
  case DIFF_C4:
    return D2DX2_C4_x_norm(f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexD2DX2_norm unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexD2DY2_norm(const Field2D &f, CELL_LOC outloc,
                                          DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_Second][1];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_y_norm(f);
    break;
  case DIFF_C4:
    return D2DX2_C4_y_norm(f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexD2DY2_norm unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexVDDX_norm(const Field3D &v, const Field3D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_Upwind][0];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_x_norm(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_x_norm(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_x_norm(v, f);
    break;
  case DIFF_U3:
    return VDDX_U3_x_norm(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_x_norm(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDX_norm unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_U2"
                        " * DIFF_C2"
                        " * DIFF_U3"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexVDDY_norm(const Field3D &v, const Field3D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_Upwind][1];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_y_norm(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_y_norm(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_y_norm(v, f);
    break;
  case DIFF_U3:
    return VDDX_U3_y_norm(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_y_norm(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDY_norm unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_U2"
                        " * DIFF_C2"
                        " * DIFF_U3"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexVDDZ_norm(const Field3D &v, const Field3D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_Upwind][2];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_z_norm(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_z_norm(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_z_norm(v, f);
    break;
  case DIFF_U3:
    return VDDX_U3_z_norm(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_z_norm(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDZ_norm unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_U2"
                        " * DIFF_C2"
                        " * DIFF_U3"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexVDDX_norm(const Field2D &v, const Field2D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_Upwind][0];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_x_norm(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_x_norm(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_x_norm(v, f);
    break;
  case DIFF_U3:
    return VDDX_U3_x_norm(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_x_norm(v, f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexVDDX_norm unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_U2"
                        " * DIFF_C2"
                        " * DIFF_U3"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexVDDY_norm(const Field2D &v, const Field2D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_Upwind][1];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_y_norm(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_y_norm(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_y_norm(v, f);
    break;
  case DIFF_U3:
    return VDDX_U3_y_norm(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_y_norm(v, f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexVDDY_norm unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_U2"
                        " * DIFF_C2"
                        " * DIFF_U3"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexFDDX_norm(const Field3D &v, const Field3D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_Flux][0];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_x_norm(v, f);
    break;
  case DIFF_C2:
    return FDDX_C2_x_norm(v, f);
    break;
  case DIFF_C4:
    return FDDX_C4_x_norm(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDX_norm unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexFDDY_norm(const Field3D &v, const Field3D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_Flux][1];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_y_norm(v, f);
    break;
  case DIFF_C2:
    return FDDX_C2_y_norm(v, f);
    break;
  case DIFF_C4:
    return FDDX_C4_y_norm(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDY_norm unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexFDDZ_norm(const Field3D &v, const Field3D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_Flux][2];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_z_norm(v, f);
    break;
  case DIFF_C2:
    return FDDX_C2_z_norm(v, f);
    break;
  case DIFF_C4:
    return FDDX_C4_z_norm(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDZ_norm unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexFDDX_norm(const Field2D &v, const Field2D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_Flux][0];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_x_norm(v, f);
    break;
  case DIFF_C2:
    return FDDX_C2_x_norm(v, f);
    break;
  case DIFF_C4:
    return FDDX_C4_x_norm(v, f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexFDDX_norm unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexFDDY_norm(const Field2D &v, const Field2D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_Flux][1];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_y_norm(v, f);
    break;
  case DIFF_C2:
    return FDDX_C2_y_norm(v, f);
    break;
  case DIFF_C4:
    return FDDX_C4_y_norm(v, f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexFDDY_norm unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexDDX_CtoL(const Field3D &f, CELL_LOC outloc,
                                        DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FirstStag][0];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_stag_x_CtoL(f);
    break;
  case DIFF_C4:
    return DDX_C4_stag_x_CtoL(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDX_CtoL unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexDDX_LtoC(const Field3D &f, CELL_LOC outloc,
                                        DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FirstStag][0];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_stag_x_LtoC(f);
    break;
  case DIFF_C4:
    return DDX_C4_stag_x_LtoC(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDX_LtoC unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexDDY_CtoL(const Field3D &f, CELL_LOC outloc,
                                        DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FirstStag][1];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_stag_y_CtoL(f);
    break;
  case DIFF_C4:
    return DDX_C4_stag_y_CtoL(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDY_CtoL unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexDDY_LtoC(const Field3D &f, CELL_LOC outloc,
                                        DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FirstStag][1];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_stag_y_LtoC(f);
    break;
  case DIFF_C4:
    return DDX_C4_stag_y_LtoC(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDY_LtoC unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexDDZ_CtoL(const Field3D &f, CELL_LOC outloc,
                                        DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FirstStag][2];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_stag_z_CtoL(f);
    break;
  case DIFF_C4:
    return DDX_C4_stag_z_CtoL(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDZ_CtoL unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexDDZ_LtoC(const Field3D &f, CELL_LOC outloc,
                                        DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FirstStag][2];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_stag_z_LtoC(f);
    break;
  case DIFF_C4:
    return DDX_C4_stag_z_LtoC(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDZ_LtoC unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexDDX_CtoL(const Field2D &f, CELL_LOC outloc,
                                        DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FirstStag][0];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_stag_x_CtoL(f);
    break;
  case DIFF_C4:
    return DDX_C4_stag_x_CtoL(f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexDDX_CtoL unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexDDX_LtoC(const Field2D &f, CELL_LOC outloc,
                                        DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FirstStag][0];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_stag_x_LtoC(f);
    break;
  case DIFF_C4:
    return DDX_C4_stag_x_LtoC(f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexDDX_LtoC unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexDDY_CtoL(const Field2D &f, CELL_LOC outloc,
                                        DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FirstStag][1];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_stag_y_CtoL(f);
    break;
  case DIFF_C4:
    return DDX_C4_stag_y_CtoL(f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexDDY_CtoL unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexDDY_LtoC(const Field2D &f, CELL_LOC outloc,
                                        DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FirstStag][1];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_stag_y_LtoC(f);
    break;
  case DIFF_C4:
    return DDX_C4_stag_y_LtoC(f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexDDY_LtoC unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexD2DX2_CtoL(const Field3D &f, CELL_LOC outloc,
                                          DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_SecondStag][0];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_stag_x_CtoL(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DX2_CtoL unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexD2DX2_LtoC(const Field3D &f, CELL_LOC outloc,
                                          DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_SecondStag][0];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_stag_x_LtoC(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DX2_LtoC unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexD2DY2_CtoL(const Field3D &f, CELL_LOC outloc,
                                          DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_SecondStag][1];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_stag_y_CtoL(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DY2_CtoL unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexD2DY2_LtoC(const Field3D &f, CELL_LOC outloc,
                                          DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_SecondStag][1];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_stag_y_LtoC(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DY2_LtoC unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexD2DZ2_CtoL(const Field3D &f, CELL_LOC outloc,
                                          DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_SecondStag][2];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_stag_z_CtoL(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DZ2_CtoL unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexD2DZ2_LtoC(const Field3D &f, CELL_LOC outloc,
                                          DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_SecondStag][2];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_stag_z_LtoC(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DZ2_LtoC unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexD2DX2_CtoL(const Field2D &f, CELL_LOC outloc,
                                          DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_SecondStag][0];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_stag_x_CtoL(f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexD2DX2_CtoL unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexD2DX2_LtoC(const Field2D &f, CELL_LOC outloc,
                                          DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_SecondStag][0];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_stag_x_LtoC(f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexD2DX2_LtoC unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexD2DY2_CtoL(const Field2D &f, CELL_LOC outloc,
                                          DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_SecondStag][1];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_stag_y_CtoL(f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexD2DY2_CtoL unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexD2DY2_LtoC(const Field2D &f, CELL_LOC outloc,
                                          DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_SecondStag][1];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_stag_y_LtoC(f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexD2DY2_LtoC unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexVDDX_CtoL(const Field3D &v, const Field3D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_UpwindStag][0];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_stag_x_CtoL(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_stag_x_CtoL(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_stag_x_CtoL(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_stag_x_CtoL(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDX_CtoL unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_U2"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexVDDX_LtoC(const Field3D &v, const Field3D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_UpwindStag][0];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_stag_x_LtoC(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_stag_x_LtoC(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_stag_x_LtoC(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_stag_x_LtoC(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDX_LtoC unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_U2"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexVDDY_CtoL(const Field3D &v, const Field3D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_UpwindStag][1];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_stag_y_CtoL(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_stag_y_CtoL(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_stag_y_CtoL(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_stag_y_CtoL(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDY_CtoL unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_U2"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexVDDY_LtoC(const Field3D &v, const Field3D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_UpwindStag][1];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_stag_y_LtoC(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_stag_y_LtoC(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_stag_y_LtoC(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_stag_y_LtoC(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDY_LtoC unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_U2"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexVDDZ_CtoL(const Field3D &v, const Field3D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_UpwindStag][2];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_stag_z_CtoL(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_stag_z_CtoL(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_stag_z_CtoL(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_stag_z_CtoL(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDZ_CtoL unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_U2"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexVDDZ_LtoC(const Field3D &v, const Field3D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_UpwindStag][2];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_stag_z_LtoC(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_stag_z_LtoC(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_stag_z_LtoC(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_stag_z_LtoC(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDZ_LtoC unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_U2"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexVDDX_CtoL(const Field2D &v, const Field2D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_UpwindStag][0];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_stag_x_CtoL(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_stag_x_CtoL(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_stag_x_CtoL(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_stag_x_CtoL(v, f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexVDDX_CtoL unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_U2"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexVDDX_LtoC(const Field2D &v, const Field2D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_UpwindStag][0];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_stag_x_LtoC(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_stag_x_LtoC(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_stag_x_LtoC(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_stag_x_LtoC(v, f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexVDDX_LtoC unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_U2"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexVDDY_CtoL(const Field2D &v, const Field2D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_UpwindStag][1];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_stag_y_CtoL(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_stag_y_CtoL(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_stag_y_CtoL(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_stag_y_CtoL(v, f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexVDDY_CtoL unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_U2"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexVDDY_LtoC(const Field2D &v, const Field2D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_UpwindStag][1];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_stag_y_LtoC(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_stag_y_LtoC(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_stag_y_LtoC(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_stag_y_LtoC(v, f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexVDDY_LtoC unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_U2"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexFDDX_CtoL(const Field3D &v, const Field3D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FluxStag][0];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_stag_x_CtoL(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDX_CtoL unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexFDDX_LtoC(const Field3D &v, const Field3D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FluxStag][0];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_stag_x_LtoC(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDX_LtoC unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexFDDY_CtoL(const Field3D &v, const Field3D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FluxStag][1];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_stag_y_CtoL(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDY_CtoL unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexFDDY_LtoC(const Field3D &v, const Field3D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FluxStag][1];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_stag_y_LtoC(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDY_LtoC unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexFDDZ_CtoL(const Field3D &v, const Field3D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FluxStag][2];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_stag_z_CtoL(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDZ_CtoL unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexFDDZ_LtoC(const Field3D &v, const Field3D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FluxStag][2];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_stag_z_LtoC(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDZ_LtoC unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexFDDX_CtoL(const Field2D &v, const Field2D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FluxStag][0];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_stag_x_CtoL(v, f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexFDDX_CtoL unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexFDDX_LtoC(const Field2D &v, const Field2D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FluxStag][0];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_stag_x_LtoC(v, f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexFDDX_LtoC unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexFDDY_CtoL(const Field2D &v, const Field2D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FluxStag][1];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_stag_y_CtoL(v, f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexFDDY_CtoL unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexFDDY_LtoC(const Field2D &v, const Field2D &f,
                                         CELL_LOC outloc, DIFF_METHOD method) const {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FluxStag][1];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_stag_y_LtoC(v, f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexFDDY_LtoC unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexDDX(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                                   REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (this->LocalNx == 1) {
    Field3D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if ((outloc == CELL_XLOW) && (f.getLocation() != CELL_XLOW)) {
    // we are going onto a staggered grid
    ASSERT1(f.getLocation() == CELL_CENTRE);
    return indexDDX_CtoL(f, outloc, method);
  } else if ((outloc != CELL_XLOW) && (f.getLocation() == CELL_XLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexDDX_LtoC(f, outloc, method);
  } else {
    ASSERT1(outloc == f.getLocation());
    return indexDDX_norm(f, outloc, method);
  }
}

// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexDDY(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                                   REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (this->LocalNy == 1) {
    Field3D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if ((outloc == CELL_YLOW) && (f.getLocation() != CELL_YLOW)) {
    // we are going onto a staggered grid
    ASSERT1(f.getLocation() == CELL_CENTRE);
    return indexDDY_CtoL(f, outloc, method);
  } else if ((outloc != CELL_YLOW) && (f.getLocation() == CELL_YLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexDDY_LtoC(f, outloc, method);
  } else {
    ASSERT1(outloc == f.getLocation());
    return indexDDY_norm(f, outloc, method);
  }
}

// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexDDZ(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                                   REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (this->LocalNz == 1) {
    Field3D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if ((outloc == CELL_ZLOW) && (f.getLocation() != CELL_ZLOW)) {
    // we are going onto a staggered grid
    ASSERT1(f.getLocation() == CELL_CENTRE);
    return indexDDZ_CtoL(f, outloc, method);
  } else if ((outloc != CELL_ZLOW) && (f.getLocation() == CELL_ZLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexDDZ_LtoC(f, outloc, method);
  } else {
    ASSERT1(outloc == f.getLocation());
    return indexDDZ_norm(f, outloc, method);
  }
}

// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexDDX(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                                   REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (this->LocalNx == 1) {
    Field2D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if ((outloc == CELL_XLOW) && (f.getLocation() != CELL_XLOW)) {
    // we are going onto a staggered grid
    ASSERT1(f.getLocation() == CELL_CENTRE);
    return indexDDX_CtoL(f, outloc, method);
  } else if ((outloc != CELL_XLOW) && (f.getLocation() == CELL_XLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexDDX_LtoC(f, outloc, method);
  } else {
    ASSERT1(outloc == f.getLocation());
    return indexDDX_norm(f, outloc, method);
  }
}

// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexDDY(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                                   REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (this->LocalNy == 1) {
    Field2D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if ((outloc == CELL_YLOW) && (f.getLocation() != CELL_YLOW)) {
    // we are going onto a staggered grid
    ASSERT1(f.getLocation() == CELL_CENTRE);
    return indexDDY_CtoL(f, outloc, method);
  } else if ((outloc != CELL_YLOW) && (f.getLocation() == CELL_YLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexDDY_LtoC(f, outloc, method);
  } else {
    ASSERT1(outloc == f.getLocation());
    return indexDDY_norm(f, outloc, method);
  }
}

// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexD2DX2(const Field3D &f, CELL_LOC outloc,
                                     DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (this->LocalNx == 1) {
    Field3D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if ((outloc == CELL_XLOW) && (f.getLocation() != CELL_XLOW)) {
    // we are going onto a staggered grid
    ASSERT1(f.getLocation() == CELL_CENTRE);
    return indexD2DX2_CtoL(f, outloc, method);
  } else if ((outloc != CELL_XLOW) && (f.getLocation() == CELL_XLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexD2DX2_LtoC(f, outloc, method);
  } else {
    ASSERT1(outloc == f.getLocation());
    return indexD2DX2_norm(f, outloc, method);
  }
}

// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexD2DY2(const Field3D &f, CELL_LOC outloc,
                                     DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (this->LocalNy == 1) {
    Field3D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if ((outloc == CELL_YLOW) && (f.getLocation() != CELL_YLOW)) {
    // we are going onto a staggered grid
    ASSERT1(f.getLocation() == CELL_CENTRE);
    return indexD2DY2_CtoL(f, outloc, method);
  } else if ((outloc != CELL_YLOW) && (f.getLocation() == CELL_YLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexD2DY2_LtoC(f, outloc, method);
  } else {
    ASSERT1(outloc == f.getLocation());
    return indexD2DY2_norm(f, outloc, method);
  }
}

// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexD2DZ2(const Field3D &f, CELL_LOC outloc,
                                     DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (this->LocalNz == 1) {
    Field3D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if ((outloc == CELL_ZLOW) && (f.getLocation() != CELL_ZLOW)) {
    // we are going onto a staggered grid
    ASSERT1(f.getLocation() == CELL_CENTRE);
    return indexD2DZ2_CtoL(f, outloc, method);
  } else if ((outloc != CELL_ZLOW) && (f.getLocation() == CELL_ZLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexD2DZ2_LtoC(f, outloc, method);
  } else {
    ASSERT1(outloc == f.getLocation());
    return indexD2DZ2_norm(f, outloc, method);
  }
}

// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexD2DX2(const Field2D &f, CELL_LOC outloc,
                                     DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (this->LocalNx == 1) {
    Field2D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if ((outloc == CELL_XLOW) && (f.getLocation() != CELL_XLOW)) {
    // we are going onto a staggered grid
    ASSERT1(f.getLocation() == CELL_CENTRE);
    return indexD2DX2_CtoL(f, outloc, method);
  } else if ((outloc != CELL_XLOW) && (f.getLocation() == CELL_XLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexD2DX2_LtoC(f, outloc, method);
  } else {
    ASSERT1(outloc == f.getLocation());
    return indexD2DX2_norm(f, outloc, method);
  }
}

// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexD2DY2(const Field2D &f, CELL_LOC outloc,
                                     DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (this->LocalNy == 1) {
    Field2D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if ((outloc == CELL_YLOW) && (f.getLocation() != CELL_YLOW)) {
    // we are going onto a staggered grid
    ASSERT1(f.getLocation() == CELL_CENTRE);
    return indexD2DY2_CtoL(f, outloc, method);
  } else if ((outloc != CELL_YLOW) && (f.getLocation() == CELL_YLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexD2DY2_LtoC(f, outloc, method);
  } else {
    ASSERT1(outloc == f.getLocation());
    return indexD2DY2_norm(f, outloc, method);
  }
}

// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexVDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                    DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (this->LocalNx == 1) {
    Field3D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if (outloc != f.getLocation()) {
    throw BoutException("AiolosMesh::index?DDX: Unhandled case for "
                        "shifting.\nf.getLocation()==outloc is required!");
  }
  if ((outloc == CELL_XLOW) && (v.getLocation() != CELL_XLOW)) {
    // we are going onto a staggered grid
    ASSERT1(v.getLocation() == CELL_CENTRE);
    return indexVDDX_CtoL(v, f, outloc, method);
  } else if ((outloc != CELL_XLOW) && (v.getLocation() == CELL_XLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexVDDX_LtoC(v, f, outloc, method);
  } else {
    ASSERT1(outloc == v.getLocation());
    return indexVDDX_norm(v, f, outloc, method);
  }
}

// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexVDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                    DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (this->LocalNy == 1) {
    Field3D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if (outloc != f.getLocation()) {
    throw BoutException("AiolosMesh::index?DDX: Unhandled case for "
                        "shifting.\nf.getLocation()==outloc is required!");
  }
  if ((outloc == CELL_YLOW) && (v.getLocation() != CELL_YLOW)) {
    // we are going onto a staggered grid
    ASSERT1(v.getLocation() == CELL_CENTRE);
    return indexVDDY_CtoL(v, f, outloc, method);
  } else if ((outloc != CELL_YLOW) && (v.getLocation() == CELL_YLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexVDDY_LtoC(v, f, outloc, method);
  } else {
    ASSERT1(outloc == v.getLocation());
    return indexVDDY_norm(v, f, outloc, method);
  }
}

// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexVDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                    DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (this->LocalNz == 1) {
    Field3D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if (outloc != f.getLocation()) {
    throw BoutException("AiolosMesh::index?DDX: Unhandled case for "
                        "shifting.\nf.getLocation()==outloc is required!");
  }
  if ((outloc == CELL_ZLOW) && (v.getLocation() != CELL_ZLOW)) {
    // we are going onto a staggered grid
    ASSERT1(v.getLocation() == CELL_CENTRE);
    return indexVDDZ_CtoL(v, f, outloc, method);
  } else if ((outloc != CELL_ZLOW) && (v.getLocation() == CELL_ZLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexVDDZ_LtoC(v, f, outloc, method);
  } else {
    ASSERT1(outloc == v.getLocation());
    return indexVDDZ_norm(v, f, outloc, method);
  }
}

// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexVDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                                    DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (this->LocalNx == 1) {
    Field2D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if (outloc != f.getLocation()) {
    throw BoutException("AiolosMesh::index?DDX: Unhandled case for "
                        "shifting.\nf.getLocation()==outloc is required!");
  }
  if ((outloc == CELL_XLOW) && (v.getLocation() != CELL_XLOW)) {
    // we are going onto a staggered grid
    ASSERT1(v.getLocation() == CELL_CENTRE);
    return indexVDDX_CtoL(v, f, outloc, method);
  } else if ((outloc != CELL_XLOW) && (v.getLocation() == CELL_XLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexVDDX_LtoC(v, f, outloc, method);
  } else {
    ASSERT1(outloc == v.getLocation());
    return indexVDDX_norm(v, f, outloc, method);
  }
}

// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexVDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                                    DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (this->LocalNy == 1) {
    Field2D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if (outloc != f.getLocation()) {
    throw BoutException("AiolosMesh::index?DDX: Unhandled case for "
                        "shifting.\nf.getLocation()==outloc is required!");
  }
  if ((outloc == CELL_YLOW) && (v.getLocation() != CELL_YLOW)) {
    // we are going onto a staggered grid
    ASSERT1(v.getLocation() == CELL_CENTRE);
    return indexVDDY_CtoL(v, f, outloc, method);
  } else if ((outloc != CELL_YLOW) && (v.getLocation() == CELL_YLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexVDDY_LtoC(v, f, outloc, method);
  } else {
    ASSERT1(outloc == v.getLocation());
    return indexVDDY_norm(v, f, outloc, method);
  }
}

// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexFDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                    DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (this->LocalNx == 1) {
    Field3D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if (outloc != f.getLocation()) {
    throw BoutException("AiolosMesh::index?DDX: Unhandled case for "
                        "shifting.\nf.getLocation()==outloc is required!");
  }
  if ((outloc == CELL_XLOW) && (v.getLocation() != CELL_XLOW)) {
    // we are going onto a staggered grid
    ASSERT1(v.getLocation() == CELL_CENTRE);
    return indexFDDX_CtoL(v, f, outloc, method);
  } else if ((outloc != CELL_XLOW) && (v.getLocation() == CELL_XLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexFDDX_LtoC(v, f, outloc, method);
  } else {
    ASSERT1(outloc == v.getLocation());
    return indexFDDX_norm(v, f, outloc, method);
  }
}

// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexFDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                    DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (this->LocalNy == 1) {
    Field3D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if (outloc != f.getLocation()) {
    throw BoutException("AiolosMesh::index?DDX: Unhandled case for "
                        "shifting.\nf.getLocation()==outloc is required!");
  }
  if ((outloc == CELL_YLOW) && (v.getLocation() != CELL_YLOW)) {
    // we are going onto a staggered grid
    ASSERT1(v.getLocation() == CELL_CENTRE);
    return indexFDDY_CtoL(v, f, outloc, method);
  } else if ((outloc != CELL_YLOW) && (v.getLocation() == CELL_YLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexFDDY_LtoC(v, f, outloc, method);
  } else {
    ASSERT1(outloc == v.getLocation());
    return indexFDDY_norm(v, f, outloc, method);
  }
}

// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative

// This file is auto-generated - do not edit!
const Field3D AiolosMesh::indexFDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                    DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (this->LocalNz == 1) {
    Field3D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if (outloc != f.getLocation()) {
    throw BoutException("AiolosMesh::index?DDX: Unhandled case for "
                        "shifting.\nf.getLocation()==outloc is required!");
  }
  if ((outloc == CELL_ZLOW) && (v.getLocation() != CELL_ZLOW)) {
    // we are going onto a staggered grid
    ASSERT1(v.getLocation() == CELL_CENTRE);
    return indexFDDZ_CtoL(v, f, outloc, method);
  } else if ((outloc != CELL_ZLOW) && (v.getLocation() == CELL_ZLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexFDDZ_LtoC(v, f, outloc, method);
  } else {
    ASSERT1(outloc == v.getLocation());
    return indexFDDZ_norm(v, f, outloc, method);
  }
}

// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexFDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                                    DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (this->LocalNx == 1) {
    Field2D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if (outloc != f.getLocation()) {
    throw BoutException("AiolosMesh::index?DDX: Unhandled case for "
                        "shifting.\nf.getLocation()==outloc is required!");
  }
  if ((outloc == CELL_XLOW) && (v.getLocation() != CELL_XLOW)) {
    // we are going onto a staggered grid
    ASSERT1(v.getLocation() == CELL_CENTRE);
    return indexFDDX_CtoL(v, f, outloc, method);
  } else if ((outloc != CELL_XLOW) && (v.getLocation() == CELL_XLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexFDDX_LtoC(v, f, outloc, method);
  } else {
    ASSERT1(outloc == v.getLocation());
    return indexFDDX_norm(v, f, outloc, method);
  }
}

// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative

// This file is auto-generated - do not edit!
const Field2D AiolosMesh::indexFDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                                    DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (this->LocalNy == 1) {
    Field2D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if (outloc != f.getLocation()) {
    throw BoutException("AiolosMesh::index?DDX: Unhandled case for "
                        "shifting.\nf.getLocation()==outloc is required!");
  }
  if ((outloc == CELL_YLOW) && (v.getLocation() != CELL_YLOW)) {
    // we are going onto a staggered grid
    ASSERT1(v.getLocation() == CELL_CENTRE);
    return indexFDDY_CtoL(v, f, outloc, method);
  } else if ((outloc != CELL_YLOW) && (v.getLocation() == CELL_YLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexFDDY_LtoC(v, f, outloc, method);
  } else {
    ASSERT1(outloc == v.getLocation());
    return indexFDDY_norm(v, f, outloc, method);
  }
}
