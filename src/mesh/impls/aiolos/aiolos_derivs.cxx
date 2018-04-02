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
#include "aiolos_stencils.hxx"
#include "aiolosmesh.hxx"

// This file is auto-generated - do not edit!
const Field3D indexDDX_norm(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
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
const Field3D indexDDY_norm(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
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
const Field3D indexDDZ_norm(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
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
const Field2D indexDDX_norm(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
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
const Field2D indexDDY_norm(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
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
const Field3D indexD2DX2_norm(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
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
const Field3D indexD2DY2_norm(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
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
const Field3D indexD2DZ2_norm(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
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
const Field2D indexD2DX2_norm(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
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
const Field2D indexD2DY2_norm(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
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
const Field3D indexVDDX_norm(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                             DIFF_METHOD method) {
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
const Field3D indexVDDY_norm(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                             DIFF_METHOD method) {
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
const Field3D indexVDDZ_norm(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                             DIFF_METHOD method) {
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
const Field2D indexVDDX_norm(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                             DIFF_METHOD method) {
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
const Field2D indexVDDY_norm(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                             DIFF_METHOD method) {
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
const Field3D indexFDDX_norm(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                             DIFF_METHOD method) {
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
const Field3D indexFDDY_norm(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                             DIFF_METHOD method) {
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
const Field3D indexFDDZ_norm(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                             DIFF_METHOD method) {
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
const Field2D indexFDDX_norm(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                             DIFF_METHOD method) {
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
const Field2D indexFDDY_norm(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                             DIFF_METHOD method) {
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
const Field3D indexDDX_on(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FirstStag][0];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_stag_x_on(f);
    break;
  case DIFF_C4:
    return DDX_C4_stag_x_on(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDX_on unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexDDX_off(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FirstStag][0];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_stag_x_off(f);
    break;
  case DIFF_C4:
    return DDX_C4_stag_x_off(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDX_off unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexDDY_on(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FirstStag][1];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_stag_y_on(f);
    break;
  case DIFF_C4:
    return DDX_C4_stag_y_on(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDY_on unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexDDY_off(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FirstStag][1];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_stag_y_off(f);
    break;
  case DIFF_C4:
    return DDX_C4_stag_y_off(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDY_off unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexDDZ_on(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FirstStag][2];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_stag_z_on(f);
    break;
  case DIFF_C4:
    return DDX_C4_stag_z_on(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDZ_on unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexDDZ_off(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FirstStag][2];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_stag_z_off(f);
    break;
  case DIFF_C4:
    return DDX_C4_stag_z_off(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDZ_off unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexDDX_on(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FirstStag][0];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_stag_x_on(f);
    break;
  case DIFF_C4:
    return DDX_C4_stag_x_on(f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexDDX_on unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexDDX_off(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FirstStag][0];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_stag_x_off(f);
    break;
  case DIFF_C4:
    return DDX_C4_stag_x_off(f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexDDX_off unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexDDY_on(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FirstStag][1];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_stag_y_on(f);
    break;
  case DIFF_C4:
    return DDX_C4_stag_y_on(f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexDDY_on unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexDDY_off(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FirstStag][1];
  }
  switch (method) {
  case DIFF_C2:
    return DDX_C2_stag_y_off(f);
    break;
  case DIFF_C4:
    return DDX_C4_stag_y_off(f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexDDY_off unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexD2DX2_on(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_SecondStag][0];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_stag_x_on(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DX2_on unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexD2DX2_off(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_SecondStag][0];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_stag_x_off(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DX2_off unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexD2DY2_on(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_SecondStag][1];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_stag_y_on(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DY2_on unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexD2DY2_off(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_SecondStag][1];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_stag_y_off(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DY2_off unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexD2DZ2_on(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_SecondStag][2];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_stag_z_on(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DZ2_on unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexD2DZ2_off(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_SecondStag][2];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_stag_z_off(f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DZ2_off unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexD2DX2_on(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_SecondStag][0];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_stag_x_on(f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexD2DX2_on unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexD2DX2_off(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_SecondStag][0];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_stag_x_off(f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexD2DX2_off unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexD2DY2_on(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_SecondStag][1];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_stag_y_on(f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexD2DY2_on unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexD2DY2_off(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_SecondStag][1];
  }
  switch (method) {
  case DIFF_C2:
    return D2DX2_C2_stag_y_off(f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexD2DY2_off unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexVDDX_on(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                           DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_UpwindStag][0];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_stag_x_on(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_stag_x_on(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_stag_x_on(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_stag_x_on(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDX_on unknown method %d.\n"
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
const Field3D indexVDDX_off(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                            DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_UpwindStag][0];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_stag_x_off(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_stag_x_off(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_stag_x_off(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_stag_x_off(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDX_off unknown method %d.\n"
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
const Field3D indexVDDY_on(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                           DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_UpwindStag][1];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_stag_y_on(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_stag_y_on(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_stag_y_on(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_stag_y_on(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDY_on unknown method %d.\n"
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
const Field3D indexVDDY_off(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                            DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_UpwindStag][1];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_stag_y_off(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_stag_y_off(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_stag_y_off(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_stag_y_off(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDY_off unknown method %d.\n"
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
const Field3D indexVDDZ_on(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                           DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_UpwindStag][2];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_stag_z_on(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_stag_z_on(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_stag_z_on(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_stag_z_on(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDZ_on unknown method %d.\n"
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
const Field3D indexVDDZ_off(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                            DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_UpwindStag][2];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_stag_z_off(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_stag_z_off(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_stag_z_off(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_stag_z_off(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDZ_off unknown method %d.\n"
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
const Field2D indexVDDX_on(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                           DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_UpwindStag][0];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_stag_x_on(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_stag_x_on(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_stag_x_on(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_stag_x_on(v, f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexVDDX_on unknown method %d.\n"
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
const Field2D indexVDDX_off(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                            DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_UpwindStag][0];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_stag_x_off(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_stag_x_off(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_stag_x_off(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_stag_x_off(v, f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexVDDX_off unknown method %d.\n"
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
const Field2D indexVDDY_on(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                           DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_UpwindStag][1];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_stag_y_on(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_stag_y_on(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_stag_y_on(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_stag_y_on(v, f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexVDDY_on unknown method %d.\n"
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
const Field2D indexVDDY_off(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                            DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_UpwindStag][1];
  }
  switch (method) {
  case DIFF_U1:
    return VDDX_U1_stag_y_off(v, f);
    break;
  case DIFF_U2:
    return VDDX_U2_stag_y_off(v, f);
    break;
  case DIFF_C2:
    return VDDX_C2_stag_y_off(v, f);
    break;
  case DIFF_C4:
    return VDDX_C4_stag_y_off(v, f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexVDDY_off unknown method %d.\n"
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
const Field3D indexFDDX_on(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                           DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FluxStag][0];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_stag_x_on(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDX_on unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexFDDX_off(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                            DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FluxStag][0];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_stag_x_off(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDX_off unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexFDDY_on(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                           DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FluxStag][1];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_stag_y_on(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDY_on unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexFDDY_off(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                            DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FluxStag][1];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_stag_y_off(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDY_off unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexFDDZ_on(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                           DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FluxStag][2];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_stag_z_on(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDZ_on unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexFDDZ_off(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                            DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FluxStag][2];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_stag_z_off(v, f);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDZ_off unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexFDDX_on(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                           DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FluxStag][0];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_stag_x_on(v, f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexFDDX_on unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexFDDX_off(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                            DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FluxStag][0];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_stag_x_off(v, f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexFDDX_off unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexFDDY_on(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                           DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FluxStag][1];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_stag_y_on(v, f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexFDDY_on unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexFDDY_off(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                            DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_stencil[AIOLOS_FluxStag][1];
  }
  switch (method) {
  case DIFF_U1:
    return FDDX_U1_stag_y_off(v, f);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexFDDY_off unknown method %d.\n"
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
    return indexDDX_on(f, outloc, method);
  } else if ((outloc != CELL_XLOW) && (f.getLocation() == CELL_XLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexDDX_off(f, outloc, method);
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
    return indexDDY_on(f, outloc, method);
  } else if ((outloc != CELL_YLOW) && (f.getLocation() == CELL_YLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexDDY_off(f, outloc, method);
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
    return indexDDZ_on(f, outloc, method);
  } else if ((outloc != CELL_ZLOW) && (f.getLocation() == CELL_ZLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexDDZ_off(f, outloc, method);
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
    return indexDDX_on(f, outloc, method);
  } else if ((outloc != CELL_XLOW) && (f.getLocation() == CELL_XLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexDDX_off(f, outloc, method);
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
    return indexDDY_on(f, outloc, method);
  } else if ((outloc != CELL_YLOW) && (f.getLocation() == CELL_YLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexDDY_off(f, outloc, method);
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
    return indexD2DX2_on(f, outloc, method);
  } else if ((outloc != CELL_XLOW) && (f.getLocation() == CELL_XLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexD2DX2_off(f, outloc, method);
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
    return indexD2DY2_on(f, outloc, method);
  } else if ((outloc != CELL_YLOW) && (f.getLocation() == CELL_YLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexD2DY2_off(f, outloc, method);
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
    return indexD2DZ2_on(f, outloc, method);
  } else if ((outloc != CELL_ZLOW) && (f.getLocation() == CELL_ZLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexD2DZ2_off(f, outloc, method);
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
    return indexD2DX2_on(f, outloc, method);
  } else if ((outloc != CELL_XLOW) && (f.getLocation() == CELL_XLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexD2DX2_off(f, outloc, method);
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
    return indexD2DY2_on(f, outloc, method);
  } else if ((outloc != CELL_YLOW) && (f.getLocation() == CELL_YLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexD2DY2_off(f, outloc, method);
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
    return indexVDDX_on(v, f, outloc, method);
  } else if ((outloc != CELL_XLOW) && (v.getLocation() == CELL_XLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexVDDX_off(v, f, outloc, method);
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
    return indexVDDY_on(v, f, outloc, method);
  } else if ((outloc != CELL_YLOW) && (v.getLocation() == CELL_YLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexVDDY_off(v, f, outloc, method);
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
    return indexVDDZ_on(v, f, outloc, method);
  } else if ((outloc != CELL_ZLOW) && (v.getLocation() == CELL_ZLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexVDDZ_off(v, f, outloc, method);
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
    return indexVDDX_on(v, f, outloc, method);
  } else if ((outloc != CELL_XLOW) && (v.getLocation() == CELL_XLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexVDDX_off(v, f, outloc, method);
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
    return indexVDDY_on(v, f, outloc, method);
  } else if ((outloc != CELL_YLOW) && (v.getLocation() == CELL_YLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexVDDY_off(v, f, outloc, method);
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
    return indexFDDX_on(v, f, outloc, method);
  } else if ((outloc != CELL_XLOW) && (v.getLocation() == CELL_XLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexFDDX_off(v, f, outloc, method);
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
    return indexFDDY_on(v, f, outloc, method);
  } else if ((outloc != CELL_YLOW) && (v.getLocation() == CELL_YLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexFDDY_off(v, f, outloc, method);
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
    return indexFDDZ_on(v, f, outloc, method);
  } else if ((outloc != CELL_ZLOW) && (v.getLocation() == CELL_ZLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexFDDZ_off(v, f, outloc, method);
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
    return indexFDDX_on(v, f, outloc, method);
  } else if ((outloc != CELL_XLOW) && (v.getLocation() == CELL_XLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexFDDX_off(v, f, outloc, method);
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
    return indexFDDY_on(v, f, outloc, method);
  } else if ((outloc != CELL_YLOW) && (v.getLocation() == CELL_YLOW)) {
    // we are coming from a staggered grid
    ASSERT1(outloc == CELL_CENTRE);
    return indexFDDY_off(v, f, outloc, method);
  } else {
    ASSERT1(outloc == v.getLocation());
    return indexFDDY_norm(v, f, outloc, method);
  }
}
