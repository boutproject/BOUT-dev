
// This file is auto-generated - do not edit!
const Field3D indexDDX_non_stag(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_x_FirstDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_C2:
    return interp_to(indexDDX_norm_DIFF_C2(f), outloc);
    break;
  case DIFF_W2:
    return interp_to(indexDDX_norm_DIFF_W2(f), outloc);
    break;
  case DIFF_C4:
    return interp_to(indexDDX_norm_DIFF_C4(f), outloc);
    break;
  case DIFF_S2:
    return interp_to(indexDDX_norm_DIFF_S2(f), outloc);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDX_non_stag unknown method %d.\n"
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
const Field3D indexDDY_non_stag(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_y_FirstDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_C2:
    return interp_to(indexDDY_norm_DIFF_C2(f), outloc);
    break;
  case DIFF_W2:
    return interp_to(indexDDY_norm_DIFF_W2(f), outloc);
    break;
  case DIFF_C4:
    return interp_to(indexDDY_norm_DIFF_C4(f), outloc);
    break;
  case DIFF_S2:
    return interp_to(indexDDY_norm_DIFF_S2(f), outloc);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDY_non_stag unknown method %d.\n"
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
const Field3D indexDDZ_non_stag(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_z_FirstDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_C2:
    return interp_to(indexDDZ_norm_DIFF_C2(f), outloc);
    break;
  case DIFF_W2:
    return interp_to(indexDDZ_norm_DIFF_W2(f), outloc);
    break;
  case DIFF_C4:
    return interp_to(indexDDZ_norm_DIFF_C4(f), outloc);
    break;
  case DIFF_S2:
    return interp_to(indexDDZ_norm_DIFF_S2(f), outloc);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDZ_non_stag unknown method %d.\n"
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
const Field2D indexDDX_non_stag(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_x_FirstDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_C2:
    return interp_to(indexDDX_norm_DIFF_C2(f), outloc);
    break;
  case DIFF_W2:
    return interp_to(indexDDX_norm_DIFF_W2(f), outloc);
    break;
  case DIFF_C4:
    return interp_to(indexDDX_norm_DIFF_C4(f), outloc);
    break;
  case DIFF_S2:
    return interp_to(indexDDX_norm_DIFF_S2(f), outloc);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexDDX_non_stag unknown method %d.\n"
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
const Field2D indexDDY_non_stag(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_y_FirstDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_C2:
    return interp_to(indexDDY_norm_DIFF_C2(f), outloc);
    break;
  case DIFF_W2:
    return interp_to(indexDDY_norm_DIFF_W2(f), outloc);
    break;
  case DIFF_C4:
    return interp_to(indexDDY_norm_DIFF_C4(f), outloc);
    break;
  case DIFF_S2:
    return interp_to(indexDDY_norm_DIFF_S2(f), outloc);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexDDY_non_stag unknown method %d.\n"
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
const Field3D indexD2DX2_non_stag(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_x_SecondDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_C2:
    return interp_to(indexD2DX2_norm_DIFF_C2(f), outloc);
    break;
  case DIFF_C4:
    return interp_to(indexD2DX2_norm_DIFF_C4(f), outloc);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DX2_non_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexD2DY2_non_stag(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_y_SecondDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_C2:
    return interp_to(indexD2DY2_norm_DIFF_C2(f), outloc);
    break;
  case DIFF_C4:
    return interp_to(indexD2DY2_norm_DIFF_C4(f), outloc);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DY2_non_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexD2DZ2_non_stag(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_z_SecondDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_C2:
    return interp_to(indexD2DZ2_norm_DIFF_C2(f), outloc);
    break;
  case DIFF_C4:
    return interp_to(indexD2DZ2_norm_DIFF_C4(f), outloc);
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DZ2_non_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexD2DX2_non_stag(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_x_SecondDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_C2:
    return interp_to(indexD2DX2_norm_DIFF_C2(f), outloc);
    break;
  case DIFF_C4:
    return interp_to(indexD2DX2_norm_DIFF_C4(f), outloc);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexD2DX2_non_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexD2DY2_non_stag(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_y_SecondDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_C2:
    return interp_to(indexD2DY2_norm_DIFF_C2(f), outloc);
    break;
  case DIFF_C4:
    return interp_to(indexD2DY2_norm_DIFF_C4(f), outloc);
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexD2DY2_non_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexVDDX_non_stag(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                 DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_x_UpwindDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_U1:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDX_norm_DIFF_U1(v, f), outloc);
    } else {
      return interp_to(
          indexVDDX_norm_DIFF_U1(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_U2:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDX_norm_DIFF_U2(v, f), outloc);
    } else {
      return interp_to(
          indexVDDX_norm_DIFF_U2(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_C2:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDX_norm_DIFF_C2(v, f), outloc);
    } else {
      return interp_to(
          indexVDDX_norm_DIFF_C2(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_U3:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDX_norm_DIFF_U3(v, f), outloc);
    } else {
      return interp_to(
          indexVDDX_norm_DIFF_U3(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_C4:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDX_norm_DIFF_C4(v, f), outloc);
    } else {
      return interp_to(
          indexVDDX_norm_DIFF_C4(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDX_non_stag unknown method %d.\n"
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
const Field3D indexVDDY_non_stag(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                 DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_y_UpwindDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_U1:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDY_norm_DIFF_U1(v, f), outloc);
    } else {
      return interp_to(
          indexVDDY_norm_DIFF_U1(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_U2:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDY_norm_DIFF_U2(v, f), outloc);
    } else {
      return interp_to(
          indexVDDY_norm_DIFF_U2(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_C2:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDY_norm_DIFF_C2(v, f), outloc);
    } else {
      return interp_to(
          indexVDDY_norm_DIFF_C2(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_U3:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDY_norm_DIFF_U3(v, f), outloc);
    } else {
      return interp_to(
          indexVDDY_norm_DIFF_U3(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_C4:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDY_norm_DIFF_C4(v, f), outloc);
    } else {
      return interp_to(
          indexVDDY_norm_DIFF_C4(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDY_non_stag unknown method %d.\n"
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
const Field3D indexVDDZ_non_stag(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                 DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_z_UpwindDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_U1:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDZ_norm_DIFF_U1(v, f), outloc);
    } else {
      return interp_to(
          indexVDDZ_norm_DIFF_U1(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_U2:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDZ_norm_DIFF_U2(v, f), outloc);
    } else {
      return interp_to(
          indexVDDZ_norm_DIFF_U2(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_C2:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDZ_norm_DIFF_C2(v, f), outloc);
    } else {
      return interp_to(
          indexVDDZ_norm_DIFF_C2(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_U3:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDZ_norm_DIFF_U3(v, f), outloc);
    } else {
      return interp_to(
          indexVDDZ_norm_DIFF_U3(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_C4:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDZ_norm_DIFF_C4(v, f), outloc);
    } else {
      return interp_to(
          indexVDDZ_norm_DIFF_C4(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDZ_non_stag unknown method %d.\n"
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
const Field2D indexVDDX_non_stag(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                                 DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_x_UpwindDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_U1:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDX_norm_DIFF_U1(v, f), outloc);
    } else {
      return interp_to(
          indexVDDX_norm_DIFF_U1(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_U2:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDX_norm_DIFF_U2(v, f), outloc);
    } else {
      return interp_to(
          indexVDDX_norm_DIFF_U2(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_C2:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDX_norm_DIFF_C2(v, f), outloc);
    } else {
      return interp_to(
          indexVDDX_norm_DIFF_C2(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_U3:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDX_norm_DIFF_U3(v, f), outloc);
    } else {
      return interp_to(
          indexVDDX_norm_DIFF_U3(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_C4:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDX_norm_DIFF_C4(v, f), outloc);
    } else {
      return interp_to(
          indexVDDX_norm_DIFF_C4(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexVDDX_non_stag unknown method %d.\n"
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
const Field2D indexVDDY_non_stag(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                                 DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_y_UpwindDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_U1:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDY_norm_DIFF_U1(v, f), outloc);
    } else {
      return interp_to(
          indexVDDY_norm_DIFF_U1(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_U2:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDY_norm_DIFF_U2(v, f), outloc);
    } else {
      return interp_to(
          indexVDDY_norm_DIFF_U2(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_C2:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDY_norm_DIFF_C2(v, f), outloc);
    } else {
      return interp_to(
          indexVDDY_norm_DIFF_C2(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_U3:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDY_norm_DIFF_U3(v, f), outloc);
    } else {
      return interp_to(
          indexVDDY_norm_DIFF_U3(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_C4:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexVDDY_norm_DIFF_C4(v, f), outloc);
    } else {
      return interp_to(
          indexVDDY_norm_DIFF_C4(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexVDDY_non_stag unknown method %d.\n"
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
const Field3D indexFDDX_non_stag(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                 DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_x_FluxDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_U1:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexFDDX_norm_DIFF_U1(v, f), outloc);
    } else {
      return interp_to(
          indexFDDX_norm_DIFF_U1(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_C2:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexFDDX_norm_DIFF_C2(v, f), outloc);
    } else {
      return interp_to(
          indexFDDX_norm_DIFF_C2(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_C4:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexFDDX_norm_DIFF_C4(v, f), outloc);
    } else {
      return interp_to(
          indexFDDX_norm_DIFF_C4(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDX_non_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexFDDY_non_stag(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                 DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_y_FluxDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_U1:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexFDDY_norm_DIFF_U1(v, f), outloc);
    } else {
      return interp_to(
          indexFDDY_norm_DIFF_U1(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_C2:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexFDDY_norm_DIFF_C2(v, f), outloc);
    } else {
      return interp_to(
          indexFDDY_norm_DIFF_C2(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_C4:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexFDDY_norm_DIFF_C4(v, f), outloc);
    } else {
      return interp_to(
          indexFDDY_norm_DIFF_C4(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDY_non_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexFDDZ_non_stag(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                 DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_z_FluxDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_U1:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexFDDZ_norm_DIFF_U1(v, f), outloc);
    } else {
      return interp_to(
          indexFDDZ_norm_DIFF_U1(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_C2:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexFDDZ_norm_DIFF_C2(v, f), outloc);
    } else {
      return interp_to(
          indexFDDZ_norm_DIFF_C2(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_C4:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexFDDZ_norm_DIFF_C4(v, f), outloc);
    } else {
      return interp_to(
          indexFDDZ_norm_DIFF_C4(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDZ_non_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexFDDX_non_stag(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                                 DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_x_FluxDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_U1:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexFDDX_norm_DIFF_U1(v, f), outloc);
    } else {
      return interp_to(
          indexFDDX_norm_DIFF_U1(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_C2:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexFDDX_norm_DIFF_C2(v, f), outloc);
    } else {
      return interp_to(
          indexFDDX_norm_DIFF_C2(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_C4:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexFDDX_norm_DIFF_C4(v, f), outloc);
    } else {
      return interp_to(
          indexFDDX_norm_DIFF_C4(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexFDDX_non_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexFDDY_non_stag(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                                 DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_y_FluxDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_U1:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexFDDY_norm_DIFF_U1(v, f), outloc);
    } else {
      return interp_to(
          indexFDDY_norm_DIFF_U1(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_C2:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexFDDY_norm_DIFF_C2(v, f), outloc);
    } else {
      return interp_to(
          indexFDDY_norm_DIFF_C2(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  case DIFF_C4:
    if (v.getLocation() == f.getLocation()) {
      return interp_to(indexFDDY_norm_DIFF_C4(v, f), outloc);
    } else {
      return interp_to(
          indexFDDY_norm_DIFF_C4(interp_to(v, CELL_CENTRE), interp_to(f, CELL_CENTRE)),
          outloc);
    }
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexFDDY_non_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexDDX_stag(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_x_FirstStagDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_C2:
    if (outloc == CELL_XLOW) {
      return indexDDX_on_DIFF_C2(interp_to(f, CELL_CENTRE));
    } else {
      return interp_to(indexDDX_off_DIFF_C2(f), outloc);
    }
    break;
  case DIFF_C4:
    if (outloc == CELL_XLOW) {
      return indexDDX_on_DIFF_C4(interp_to(f, CELL_CENTRE));
    } else {
      return interp_to(indexDDX_off_DIFF_C4(f), outloc);
    }
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDX_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexDDY_stag(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_y_FirstStagDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_C2:
    if (outloc == CELL_YLOW) {
      return indexDDY_on_DIFF_C2(interp_to(f, CELL_CENTRE));
    } else {
      return interp_to(indexDDY_off_DIFF_C2(f), outloc);
    }
    break;
  case DIFF_C4:
    if (outloc == CELL_YLOW) {
      return indexDDY_on_DIFF_C4(interp_to(f, CELL_CENTRE));
    } else {
      return interp_to(indexDDY_off_DIFF_C4(f), outloc);
    }
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDY_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexDDZ_stag(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_z_FirstStagDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_C2:
    if (outloc == CELL_ZLOW) {
      return indexDDZ_on_DIFF_C2(interp_to(f, CELL_CENTRE));
    } else {
      return interp_to(indexDDZ_off_DIFF_C2(f), outloc);
    }
    break;
  case DIFF_C4:
    if (outloc == CELL_ZLOW) {
      return indexDDZ_on_DIFF_C4(interp_to(f, CELL_CENTRE));
    } else {
      return interp_to(indexDDZ_off_DIFF_C4(f), outloc);
    }
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexDDZ_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexDDX_stag(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_x_FirstStagDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_C2:
    if (outloc == CELL_XLOW) {
      return indexDDX_on_DIFF_C2(interp_to(f, CELL_CENTRE));
    } else {
      return interp_to(indexDDX_off_DIFF_C2(f), outloc);
    }
    break;
  case DIFF_C4:
    if (outloc == CELL_XLOW) {
      return indexDDX_on_DIFF_C4(interp_to(f, CELL_CENTRE));
    } else {
      return interp_to(indexDDX_off_DIFF_C4(f), outloc);
    }
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexDDX_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexDDY_stag(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_y_FirstStagDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_C2:
    if (outloc == CELL_YLOW) {
      return indexDDY_on_DIFF_C2(interp_to(f, CELL_CENTRE));
    } else {
      return interp_to(indexDDY_off_DIFF_C2(f), outloc);
    }
    break;
  case DIFF_C4:
    if (outloc == CELL_YLOW) {
      return indexDDY_on_DIFF_C4(interp_to(f, CELL_CENTRE));
    } else {
      return interp_to(indexDDY_off_DIFF_C4(f), outloc);
    }
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexDDY_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        " * DIFF_C4"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexD2DX2_stag(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_x_SecondStagDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_C2:
    if (outloc == CELL_XLOW) {
      return indexD2DX2_on_DIFF_C2(interp_to(f, CELL_CENTRE));
    } else {
      return interp_to(indexD2DX2_off_DIFF_C2(f), outloc);
    }
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DX2_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexD2DY2_stag(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_y_SecondStagDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_C2:
    if (outloc == CELL_YLOW) {
      return indexD2DY2_on_DIFF_C2(interp_to(f, CELL_CENTRE));
    } else {
      return interp_to(indexD2DY2_off_DIFF_C2(f), outloc);
    }
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DY2_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexD2DZ2_stag(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_z_SecondStagDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_C2:
    if (outloc == CELL_ZLOW) {
      return indexD2DZ2_on_DIFF_C2(interp_to(f, CELL_CENTRE));
    } else {
      return interp_to(indexD2DZ2_off_DIFF_C2(f), outloc);
    }
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexD2DZ2_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexD2DX2_stag(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_x_SecondStagDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_C2:
    if (outloc == CELL_XLOW) {
      return indexD2DX2_on_DIFF_C2(interp_to(f, CELL_CENTRE));
    } else {
      return interp_to(indexD2DX2_off_DIFF_C2(f), outloc);
    }
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexD2DX2_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexD2DY2_stag(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_y_SecondStagDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_C2:
    if (outloc == CELL_YLOW) {
      return indexD2DY2_on_DIFF_C2(interp_to(f, CELL_CENTRE));
    } else {
      return interp_to(indexD2DY2_off_DIFF_C2(f), outloc);
    }
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexD2DY2_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_C2"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexVDDX_stag(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                             DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_x_UpwindStagDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_U1:
    if (outloc == CELL_XLOW) {
      return indexVDDX_on_DIFF_U1(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexVDDX_off_DIFF_U1(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  case DIFF_U2:
    if (outloc == CELL_XLOW) {
      return indexVDDX_on_DIFF_U2(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexVDDX_off_DIFF_U2(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  case DIFF_C2:
    if (outloc == CELL_XLOW) {
      return indexVDDX_on_DIFF_C2(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexVDDX_off_DIFF_C2(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  case DIFF_C4:
    if (outloc == CELL_XLOW) {
      return indexVDDX_on_DIFF_C4(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexVDDX_off_DIFF_C4(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDX_stag unknown method %d.\n"
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
const Field3D indexVDDY_stag(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                             DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_y_UpwindStagDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_U1:
    if (outloc == CELL_YLOW) {
      return indexVDDY_on_DIFF_U1(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexVDDY_off_DIFF_U1(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  case DIFF_U2:
    if (outloc == CELL_YLOW) {
      return indexVDDY_on_DIFF_U2(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexVDDY_off_DIFF_U2(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  case DIFF_C2:
    if (outloc == CELL_YLOW) {
      return indexVDDY_on_DIFF_C2(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexVDDY_off_DIFF_C2(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  case DIFF_C4:
    if (outloc == CELL_YLOW) {
      return indexVDDY_on_DIFF_C4(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexVDDY_off_DIFF_C4(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDY_stag unknown method %d.\n"
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
const Field3D indexVDDZ_stag(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                             DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_z_UpwindStagDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_U1:
    if (outloc == CELL_ZLOW) {
      return indexVDDZ_on_DIFF_U1(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexVDDZ_off_DIFF_U1(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  case DIFF_U2:
    if (outloc == CELL_ZLOW) {
      return indexVDDZ_on_DIFF_U2(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexVDDZ_off_DIFF_U2(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  case DIFF_C2:
    if (outloc == CELL_ZLOW) {
      return indexVDDZ_on_DIFF_C2(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexVDDZ_off_DIFF_C2(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  case DIFF_C4:
    if (outloc == CELL_ZLOW) {
      return indexVDDZ_on_DIFF_C4(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexVDDZ_off_DIFF_C4(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexVDDZ_stag unknown method %d.\n"
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
const Field2D indexVDDX_stag(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                             DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_x_UpwindStagDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_U1:
    if (outloc == CELL_XLOW) {
      return indexVDDX_on_DIFF_U1(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexVDDX_off_DIFF_U1(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  case DIFF_U2:
    if (outloc == CELL_XLOW) {
      return indexVDDX_on_DIFF_U2(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexVDDX_off_DIFF_U2(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  case DIFF_C2:
    if (outloc == CELL_XLOW) {
      return indexVDDX_on_DIFF_C2(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexVDDX_off_DIFF_C2(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  case DIFF_C4:
    if (outloc == CELL_XLOW) {
      return indexVDDX_on_DIFF_C4(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexVDDX_off_DIFF_C4(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexVDDX_stag unknown method %d.\n"
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
const Field2D indexVDDY_stag(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                             DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_y_UpwindStagDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_U1:
    if (outloc == CELL_YLOW) {
      return indexVDDY_on_DIFF_U1(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexVDDY_off_DIFF_U1(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  case DIFF_U2:
    if (outloc == CELL_YLOW) {
      return indexVDDY_on_DIFF_U2(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexVDDY_off_DIFF_U2(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  case DIFF_C2:
    if (outloc == CELL_YLOW) {
      return indexVDDY_on_DIFF_C2(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexVDDY_off_DIFF_C2(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  case DIFF_C4:
    if (outloc == CELL_YLOW) {
      return indexVDDY_on_DIFF_C4(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexVDDY_off_DIFF_C4(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexVDDY_stag unknown method %d.\n"
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
const Field3D indexFDDX_stag(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                             DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_x_FluxStagDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_U1:
    if (outloc == CELL_XLOW) {
      return indexFDDX_on_DIFF_U1(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexFDDX_off_DIFF_U1(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDX_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexFDDY_stag(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                             DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_y_FluxStagDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_U1:
    if (outloc == CELL_YLOW) {
      return indexFDDY_on_DIFF_U1(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexFDDY_off_DIFF_U1(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDY_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field3D indexFDDZ_stag(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                             DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_z_FluxStagDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_U1:
    if (outloc == CELL_ZLOW) {
      return indexFDDZ_on_DIFF_U1(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexFDDZ_off_DIFF_U1(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  default:
    throw BoutException("Field3D AiolosMesh::indexFDDZ_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexFDDX_stag(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                             DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_x_FluxStagDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_U1:
    if (outloc == CELL_XLOW) {
      return indexFDDX_on_DIFF_U1(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexFDDX_off_DIFF_U1(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexFDDX_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
const Field2D indexFDDY_stag(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                             DIFF_METHOD method) {
  if (method == DIFF_DEFAULT) {
    method = default_y_FluxStagDeriv;
  }
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  switch (method) {
  case DIFF_U1:
    if (outloc == CELL_YLOW) {
      return indexFDDY_on_DIFF_U1(interp_to(v, CELL_CENTRE), f);
    } else {
      return interp_to(indexFDDY_off_DIFF_U1(v, interp_to(f, CELL_CENTRE)), outloc);
    }
    break;
  default:
    throw BoutException("Field2D AiolosMesh::indexFDDY_stag unknown method %d.\n"
                        "Supported methods are"
                        " * DIFF_U1"
                        "\nNote FFTs are not (yet) supported.",
                        method);
  }; // end switch
}

// This file is auto-generated - do not edit!
// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative
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
  if ((outloc == CELL_XLOW) != (f.getLocation() == CELL_XLOW)) {
    // we are going onto a staggered grid or coming from one
    return indexDDX_stag(f, outloc, method);
  } else {
    return indexDDX_non_stag(f, outloc, method);
  }
}

// This file is auto-generated - do not edit!
// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative
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
  if ((outloc == CELL_YLOW) != (f.getLocation() == CELL_YLOW)) {
    // we are going onto a staggered grid or coming from one
    return indexDDY_stag(f, outloc, method);
  } else {
    return indexDDY_non_stag(f, outloc, method);
  }
}

// This file is auto-generated - do not edit!
// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative
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
  if ((outloc == CELL_ZLOW) != (f.getLocation() == CELL_ZLOW)) {
    // we are going onto a staggered grid or coming from one
    return indexDDZ_stag(f, outloc, method);
  } else {
    return indexDDZ_non_stag(f, outloc, method);
  }
}

// This file is auto-generated - do not edit!
// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative
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
  if ((outloc == CELL_XLOW) != (f.getLocation() == CELL_XLOW)) {
    // we are going onto a staggered grid or coming from one
    return indexDDX_stag(f, outloc, method);
  } else {
    return indexDDX_non_stag(f, outloc, method);
  }
}

// This file is auto-generated - do not edit!
// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative
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
  if ((outloc == CELL_YLOW) != (f.getLocation() == CELL_YLOW)) {
    // we are going onto a staggered grid or coming from one
    return indexDDY_stag(f, outloc, method);
  } else {
    return indexDDY_non_stag(f, outloc, method);
  }
}

// This file is auto-generated - do not edit!
// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative
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
  if ((outloc == CELL_XLOW) != (f.getLocation() == CELL_XLOW)) {
    // we are going onto a staggered grid or coming from one
    return indexD2DX2_stag(f, outloc, method);
  } else {
    return indexD2DX2_non_stag(f, outloc, method);
  }
}

// This file is auto-generated - do not edit!
// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative
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
  if ((outloc == CELL_YLOW) != (f.getLocation() == CELL_YLOW)) {
    // we are going onto a staggered grid or coming from one
    return indexD2DY2_stag(f, outloc, method);
  } else {
    return indexD2DY2_non_stag(f, outloc, method);
  }
}

// This file is auto-generated - do not edit!
// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative
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
  if ((outloc == CELL_ZLOW) != (f.getLocation() == CELL_ZLOW)) {
    // we are going onto a staggered grid or coming from one
    return indexD2DZ2_stag(f, outloc, method);
  } else {
    return indexD2DZ2_non_stag(f, outloc, method);
  }
}

// This file is auto-generated - do not edit!
// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative
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
  if ((outloc == CELL_XLOW) != (f.getLocation() == CELL_XLOW)) {
    // we are going onto a staggered grid or coming from one
    return indexD2DX2_stag(f, outloc, method);
  } else {
    return indexD2DX2_non_stag(f, outloc, method);
  }
}

// This file is auto-generated - do not edit!
// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative
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
  if ((outloc == CELL_YLOW) != (f.getLocation() == CELL_YLOW)) {
    // we are going onto a staggered grid or coming from one
    return indexD2DY2_stag(f, outloc, method);
  } else {
    return indexD2DY2_non_stag(f, outloc, method);
  }
}

// This file is auto-generated - do not edit!
// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative
const Field3D AiolosMesh::indexVDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                    DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (outloc != f.getLocation()) {
    throw BoutException("AiolosMesh::index?DDX: Unhandled case for "
                        "shifting.\nf.getLocation()==outloc is required!");
  }
  if (this->LocalNx == 1) {
    Field3D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if ((outloc == CELL_XLOW) != (v.getLocation() == CELL_XLOW)) {
    // we are going onto a staggered grid or coming from one
    return indexVDDX_stag(v, f, outloc, method);
  } else {
    return indexVDDX_non_stag(v, f, outloc, method);
  }
}

// This file is auto-generated - do not edit!
// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative
const Field3D AiolosMesh::indexVDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                    DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (outloc != f.getLocation()) {
    throw BoutException("AiolosMesh::index?DDX: Unhandled case for "
                        "shifting.\nf.getLocation()==outloc is required!");
  }
  if (this->LocalNy == 1) {
    Field3D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if ((outloc == CELL_YLOW) != (v.getLocation() == CELL_YLOW)) {
    // we are going onto a staggered grid or coming from one
    return indexVDDY_stag(v, f, outloc, method);
  } else {
    return indexVDDY_non_stag(v, f, outloc, method);
  }
}

// This file is auto-generated - do not edit!
// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative
const Field3D AiolosMesh::indexVDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                    DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (outloc != f.getLocation()) {
    throw BoutException("AiolosMesh::index?DDX: Unhandled case for "
                        "shifting.\nf.getLocation()==outloc is required!");
  }
  if (this->LocalNz == 1) {
    Field3D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if ((outloc == CELL_ZLOW) != (v.getLocation() == CELL_ZLOW)) {
    // we are going onto a staggered grid or coming from one
    return indexVDDZ_stag(v, f, outloc, method);
  } else {
    return indexVDDZ_non_stag(v, f, outloc, method);
  }
}

// This file is auto-generated - do not edit!
// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative
const Field2D AiolosMesh::indexVDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                                    DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (outloc != f.getLocation()) {
    throw BoutException("AiolosMesh::index?DDX: Unhandled case for "
                        "shifting.\nf.getLocation()==outloc is required!");
  }
  if (this->LocalNx == 1) {
    Field2D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if ((outloc == CELL_XLOW) != (v.getLocation() == CELL_XLOW)) {
    // we are going onto a staggered grid or coming from one
    return indexVDDX_stag(v, f, outloc, method);
  } else {
    return indexVDDX_non_stag(v, f, outloc, method);
  }
}

// This file is auto-generated - do not edit!
// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative
const Field2D AiolosMesh::indexVDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                                    DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (outloc != f.getLocation()) {
    throw BoutException("AiolosMesh::index?DDX: Unhandled case for "
                        "shifting.\nf.getLocation()==outloc is required!");
  }
  if (this->LocalNy == 1) {
    Field2D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if ((outloc == CELL_YLOW) != (v.getLocation() == CELL_YLOW)) {
    // we are going onto a staggered grid or coming from one
    return indexVDDY_stag(v, f, outloc, method);
  } else {
    return indexVDDY_non_stag(v, f, outloc, method);
  }
}

// This file is auto-generated - do not edit!
// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative
const Field3D AiolosMesh::indexFDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                    DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (outloc != f.getLocation()) {
    throw BoutException("AiolosMesh::index?DDX: Unhandled case for "
                        "shifting.\nf.getLocation()==outloc is required!");
  }
  if (this->LocalNx == 1) {
    Field3D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if ((outloc == CELL_XLOW) != (v.getLocation() == CELL_XLOW)) {
    // we are going onto a staggered grid or coming from one
    return indexFDDX_stag(v, f, outloc, method);
  } else {
    return indexFDDX_non_stag(v, f, outloc, method);
  }
}

// This file is auto-generated - do not edit!
// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative
const Field3D AiolosMesh::indexFDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                    DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (outloc != f.getLocation()) {
    throw BoutException("AiolosMesh::index?DDX: Unhandled case for "
                        "shifting.\nf.getLocation()==outloc is required!");
  }
  if (this->LocalNy == 1) {
    Field3D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if ((outloc == CELL_YLOW) != (v.getLocation() == CELL_YLOW)) {
    // we are going onto a staggered grid or coming from one
    return indexFDDY_stag(v, f, outloc, method);
  } else {
    return indexFDDY_non_stag(v, f, outloc, method);
  }
}

// This file is auto-generated - do not edit!
// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative
const Field3D AiolosMesh::indexFDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                    DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (outloc != f.getLocation()) {
    throw BoutException("AiolosMesh::index?DDX: Unhandled case for "
                        "shifting.\nf.getLocation()==outloc is required!");
  }
  if (this->LocalNz == 1) {
    Field3D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if ((outloc == CELL_ZLOW) != (v.getLocation() == CELL_ZLOW)) {
    // we are going onto a staggered grid or coming from one
    return indexFDDZ_stag(v, f, outloc, method);
  } else {
    return indexFDDZ_non_stag(v, f, outloc, method);
  }
}

// This file is auto-generated - do not edit!
// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative
const Field2D AiolosMesh::indexFDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                                    DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (outloc != f.getLocation()) {
    throw BoutException("AiolosMesh::index?DDX: Unhandled case for "
                        "shifting.\nf.getLocation()==outloc is required!");
  }
  if (this->LocalNx == 1) {
    Field2D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if ((outloc == CELL_XLOW) != (v.getLocation() == CELL_XLOW)) {
    // we are going onto a staggered grid or coming from one
    return indexFDDX_stag(v, f, outloc, method);
  } else {
    return indexFDDX_non_stag(v, f, outloc, method);
  }
}

// This file is auto-generated - do not edit!
// Do check the input parameters. Further decide on whether or not we are doing a
// staggered derivative or a non-staaggered derivative
const Field2D AiolosMesh::indexFDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                                    DIFF_METHOD method, REGION ignored) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  if (outloc != f.getLocation()) {
    throw BoutException("AiolosMesh::index?DDX: Unhandled case for "
                        "shifting.\nf.getLocation()==outloc is required!");
  }
  if (this->LocalNy == 1) {
    Field2D result{0., this};
    result.setLocation(outloc);
    return result;
  }
  if ((outloc == CELL_YLOW) != (v.getLocation() == CELL_YLOW)) {
    // we are going onto a staggered grid or coming from one
    return indexFDDY_stag(v, f, outloc, method);
  } else {
    return indexFDDY_non_stag(v, f, outloc, method);
  }
}
