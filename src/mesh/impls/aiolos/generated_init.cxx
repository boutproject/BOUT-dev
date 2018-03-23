
// This file is auto-generated - do not edit!
DIFF_METHOD default_x_FirstDeriv;
DIFF_METHOD default_x_FirstStagDeriv;
DIFF_METHOD default_x_SecondDeriv;
DIFF_METHOD default_x_SecondStagDeriv;
DIFF_METHOD default_x_UpwindDeriv;
DIFF_METHOD default_x_UpwindStagDeriv;
DIFF_METHOD default_x_FluxDeriv;
DIFF_METHOD default_x_FluxStagDeriv;

// This file is auto-generated - do not edit!
DIFF_METHOD default_y_FirstDeriv;
DIFF_METHOD default_y_FirstStagDeriv;
DIFF_METHOD default_y_SecondDeriv;
DIFF_METHOD default_y_SecondStagDeriv;
DIFF_METHOD default_y_UpwindDeriv;
DIFF_METHOD default_y_UpwindStagDeriv;
DIFF_METHOD default_y_FluxDeriv;
DIFF_METHOD default_y_FluxStagDeriv;

// This file is auto-generated - do not edit!
DIFF_METHOD default_z_FirstDeriv;
DIFF_METHOD default_z_FirstStagDeriv;
DIFF_METHOD default_z_SecondDeriv;
DIFF_METHOD default_z_SecondStagDeriv;
DIFF_METHOD default_z_UpwindDeriv;
DIFF_METHOD default_z_UpwindStagDeriv;
DIFF_METHOD default_z_FluxDeriv;
DIFF_METHOD default_z_FluxStagDeriv;

// This file is auto-generated - do not edit!
void AiolosMesh::derivs_init(Options *option) {
  std::string name;
  Options *dirOption;
  output.write("\tSetting derivatives for direction x:\n");
  dirOption = option->getSection("ddx");

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddx and First
  if (dirOption->isSet("First")) {
    dirOption->get("First", name, "C2");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "C2");
  } else {
    name = "C2";
  }
  if (strcasecmp(name.c_str(), "C2") == 0) {
    default_x_FirstDeriv = DIFF_C2;
    output.write("	          First : Second order central\n");
  } else if (strcasecmp(name.c_str(), "W2") == 0) {
    default_x_FirstDeriv = DIFF_W2;
    output.write("	          First : Second order WENO\n");
  } else if (strcasecmp(name.c_str(), "C4") == 0) {
    default_x_FirstDeriv = DIFF_C4;
    output.write("	          First : Fourth order central\n");
  } else if (strcasecmp(name.c_str(), "S2") == 0) {
    default_x_FirstDeriv = DIFF_S2;
    output.write("	          First : Smoothing 2nd order\n");
  } else {
    throw BoutException("Dont't know what diff method to use for First (direction x, "
                        "tried to use %s)!\nOptions are:\n * C2: Second order central\n "
                        "* W2: Second order WENO\n * C4: Fourth order central\n * S2: "
                        "Smoothing 2nd order",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddx and FirstStag
  if (dirOption->isSet("FirstStag")) {
    dirOption->get("FirstStag", name, "C2");
  } else if (dirOption->isSet("First")) {
    dirOption->get("First", name, "C2");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "C2");
  } else {
    name = "C2";
  }
  if (strcasecmp(name.c_str(), "C2") == 0) {
    default_x_FirstStagDeriv = DIFF_C2;
    output.write("	      FirstStag : Second order central\n");
  } else if (strcasecmp(name.c_str(), "C4") == 0) {
    default_x_FirstStagDeriv = DIFF_C4;
    output.write("	      FirstStag : Fourth order central\n");
  } else {
    throw BoutException("Dont't know what diff method to use for FirstStag (direction x, "
                        "tried to use %s)!\nOptions are:\n * C2: Second order central\n "
                        "* C4: Fourth order central",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddx and Second
  if (dirOption->isSet("Second")) {
    dirOption->get("Second", name, "C2");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "C2");
  } else {
    name = "C2";
  }
  if (strcasecmp(name.c_str(), "C2") == 0) {
    default_x_SecondDeriv = DIFF_C2;
    output.write("	         Second : Second order central\n");
  } else if (strcasecmp(name.c_str(), "C4") == 0) {
    default_x_SecondDeriv = DIFF_C4;
    output.write("	         Second : Fourth order central\n");
  } else {
    throw BoutException("Dont't know what diff method to use for Second (direction x, "
                        "tried to use %s)!\nOptions are:\n * C2: Second order central\n "
                        "* C4: Fourth order central",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddx and SecondStag
  if (dirOption->isSet("SecondStag")) {
    dirOption->get("SecondStag", name, "C2");
  } else if (dirOption->isSet("Second")) {
    dirOption->get("Second", name, "C2");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "C2");
  } else {
    name = "C2";
  }
  if (strcasecmp(name.c_str(), "C2") == 0) {
    default_x_SecondStagDeriv = DIFF_C2;
    output.write("	     SecondStag : Second order central\n");
  } else {
    throw BoutException("Dont't know what diff method to use for SecondStag (direction "
                        "x, tried to use %s)!\nOptions are:\n * C2: Second order central",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddx and Upwind
  if (dirOption->isSet("Upwind")) {
    dirOption->get("Upwind", name, "U1");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "U1");
  } else {
    name = "U1";
  }
  if (strcasecmp(name.c_str(), "U1") == 0) {
    default_x_UpwindDeriv = DIFF_U1;
    output.write("	         Upwind : First order upwinding\n");
  } else if (strcasecmp(name.c_str(), "U2") == 0) {
    default_x_UpwindDeriv = DIFF_U2;
    output.write("	         Upwind : Second order upwinding\n");
  } else if (strcasecmp(name.c_str(), "C2") == 0) {
    default_x_UpwindDeriv = DIFF_C2;
    output.write("	         Upwind : Second order central\n");
  } else if (strcasecmp(name.c_str(), "U3") == 0) {
    default_x_UpwindDeriv = DIFF_U3;
    output.write("	         Upwind : Third order upwinding\n");
  } else if (strcasecmp(name.c_str(), "U4") == 0) {
    default_x_UpwindDeriv = DIFF_U3;
    output.write(
        "	         Upwind : Third order upwinding (Can't do 4th order yet).\n");
  } else if (strcasecmp(name.c_str(), "C4") == 0) {
    default_x_UpwindDeriv = DIFF_C4;
    output.write("	         Upwind : Fourth order central\n");
  } else {
    throw BoutException("Dont't know what diff method to use for Upwind (direction x, "
                        "tried to use %s)!\nOptions are:\n * U1: First order upwinding\n "
                        "* U2: Second order upwinding\n * C2: Second order central\n * "
                        "U3: Third order upwinding\n * U4: Third order upwinding (Can't "
                        "do 4th order yet).\n * C4: Fourth order central",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddx and UpwindStag
  if (dirOption->isSet("UpwindStag")) {
    dirOption->get("UpwindStag", name, "U1");
  } else if (dirOption->isSet("Upwind")) {
    dirOption->get("Upwind", name, "U1");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "U1");
  } else {
    name = "U1";
  }
  if (strcasecmp(name.c_str(), "U1") == 0) {
    default_x_UpwindStagDeriv = DIFF_U1;
    output.write("	     UpwindStag : First order upwinding\n");
  } else if (strcasecmp(name.c_str(), "U2") == 0) {
    default_x_UpwindStagDeriv = DIFF_U2;
    output.write("	     UpwindStag : Second order upwinding\n");
  } else if (strcasecmp(name.c_str(), "C2") == 0) {
    default_x_UpwindStagDeriv = DIFF_C2;
    output.write("	     UpwindStag : Second order central\n");
  } else if (strcasecmp(name.c_str(), "C4") == 0) {
    default_x_UpwindStagDeriv = DIFF_C4;
    output.write("	     UpwindStag : Fourth order central\n");
  } else {
    throw BoutException("Dont't know what diff method to use for UpwindStag (direction "
                        "x, tried to use %s)!\nOptions are:\n * U1: First order "
                        "upwinding\n * U2: Second order upwinding\n * C2: Second order "
                        "central\n * C4: Fourth order central",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddx and Flux
  if (dirOption->isSet("Flux")) {
    dirOption->get("Flux", name, "U1");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "U1");
  } else {
    name = "U1";
  }
  if (strcasecmp(name.c_str(), "U1") == 0) {
    default_x_FluxDeriv = DIFF_U1;
    output.write("	           Flux : First order upwinding\n");
  } else if (strcasecmp(name.c_str(), "C2") == 0) {
    default_x_FluxDeriv = DIFF_C2;
    output.write("	           Flux : Second order central\n");
  } else if (strcasecmp(name.c_str(), "C4") == 0) {
    default_x_FluxDeriv = DIFF_C4;
    output.write("	           Flux : Fourth order central\n");
  } else {
    throw BoutException("Dont't know what diff method to use for Flux (direction x, "
                        "tried to use %s)!\nOptions are:\n * U1: First order upwinding\n "
                        "* C2: Second order central\n * C4: Fourth order central",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddx and FluxStag
  if (dirOption->isSet("FluxStag")) {
    dirOption->get("FluxStag", name, "U1");
  } else if (dirOption->isSet("Flux")) {
    dirOption->get("Flux", name, "U1");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "U1");
  } else {
    name = "U1";
  }
  if (strcasecmp(name.c_str(), "U1") == 0) {
    default_x_FluxStagDeriv = DIFF_U1;
    output.write("	       FluxStag : First order upwinding\n");
  } else {
    throw BoutException("Dont't know what diff method to use for FluxStag (direction x, "
                        "tried to use %s)!\nOptions are:\n * U1: First order upwinding",
                        name.c_str());
  }
  output.write("\tSetting derivatives for direction y:\n");
  dirOption = option->getSection("ddy");

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddy and First
  if (dirOption->isSet("First")) {
    dirOption->get("First", name, "C2");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "C2");
  } else {
    name = "C2";
  }
  if (strcasecmp(name.c_str(), "C2") == 0) {
    default_y_FirstDeriv = DIFF_C2;
    output.write("	          First : Second order central\n");
  } else if (strcasecmp(name.c_str(), "W2") == 0) {
    default_y_FirstDeriv = DIFF_W2;
    output.write("	          First : Second order WENO\n");
  } else if (strcasecmp(name.c_str(), "C4") == 0) {
    default_y_FirstDeriv = DIFF_C4;
    output.write("	          First : Fourth order central\n");
  } else if (strcasecmp(name.c_str(), "S2") == 0) {
    default_y_FirstDeriv = DIFF_S2;
    output.write("	          First : Smoothing 2nd order\n");
  } else {
    throw BoutException("Dont't know what diff method to use for First (direction y, "
                        "tried to use %s)!\nOptions are:\n * C2: Second order central\n "
                        "* W2: Second order WENO\n * C4: Fourth order central\n * S2: "
                        "Smoothing 2nd order",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddy and FirstStag
  if (dirOption->isSet("FirstStag")) {
    dirOption->get("FirstStag", name, "C2");
  } else if (dirOption->isSet("First")) {
    dirOption->get("First", name, "C2");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "C2");
  } else {
    name = "C2";
  }
  if (strcasecmp(name.c_str(), "C2") == 0) {
    default_y_FirstStagDeriv = DIFF_C2;
    output.write("	      FirstStag : Second order central\n");
  } else if (strcasecmp(name.c_str(), "C4") == 0) {
    default_y_FirstStagDeriv = DIFF_C4;
    output.write("	      FirstStag : Fourth order central\n");
  } else {
    throw BoutException("Dont't know what diff method to use for FirstStag (direction y, "
                        "tried to use %s)!\nOptions are:\n * C2: Second order central\n "
                        "* C4: Fourth order central",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddy and Second
  if (dirOption->isSet("Second")) {
    dirOption->get("Second", name, "C2");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "C2");
  } else {
    name = "C2";
  }
  if (strcasecmp(name.c_str(), "C2") == 0) {
    default_y_SecondDeriv = DIFF_C2;
    output.write("	         Second : Second order central\n");
  } else if (strcasecmp(name.c_str(), "C4") == 0) {
    default_y_SecondDeriv = DIFF_C4;
    output.write("	         Second : Fourth order central\n");
  } else {
    throw BoutException("Dont't know what diff method to use for Second (direction y, "
                        "tried to use %s)!\nOptions are:\n * C2: Second order central\n "
                        "* C4: Fourth order central",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddy and SecondStag
  if (dirOption->isSet("SecondStag")) {
    dirOption->get("SecondStag", name, "C2");
  } else if (dirOption->isSet("Second")) {
    dirOption->get("Second", name, "C2");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "C2");
  } else {
    name = "C2";
  }
  if (strcasecmp(name.c_str(), "C2") == 0) {
    default_y_SecondStagDeriv = DIFF_C2;
    output.write("	     SecondStag : Second order central\n");
  } else {
    throw BoutException("Dont't know what diff method to use for SecondStag (direction "
                        "y, tried to use %s)!\nOptions are:\n * C2: Second order central",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddy and Upwind
  if (dirOption->isSet("Upwind")) {
    dirOption->get("Upwind", name, "U1");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "U1");
  } else {
    name = "U1";
  }
  if (strcasecmp(name.c_str(), "U1") == 0) {
    default_y_UpwindDeriv = DIFF_U1;
    output.write("	         Upwind : First order upwinding\n");
  } else if (strcasecmp(name.c_str(), "U2") == 0) {
    default_y_UpwindDeriv = DIFF_U2;
    output.write("	         Upwind : Second order upwinding\n");
  } else if (strcasecmp(name.c_str(), "C2") == 0) {
    default_y_UpwindDeriv = DIFF_C2;
    output.write("	         Upwind : Second order central\n");
  } else if (strcasecmp(name.c_str(), "U3") == 0) {
    default_y_UpwindDeriv = DIFF_U3;
    output.write("	         Upwind : Third order upwinding\n");
  } else if (strcasecmp(name.c_str(), "U4") == 0) {
    default_y_UpwindDeriv = DIFF_U3;
    output.write(
        "	         Upwind : Third order upwinding (Can't do 4th order yet).\n");
  } else if (strcasecmp(name.c_str(), "C4") == 0) {
    default_y_UpwindDeriv = DIFF_C4;
    output.write("	         Upwind : Fourth order central\n");
  } else {
    throw BoutException("Dont't know what diff method to use for Upwind (direction y, "
                        "tried to use %s)!\nOptions are:\n * U1: First order upwinding\n "
                        "* U2: Second order upwinding\n * C2: Second order central\n * "
                        "U3: Third order upwinding\n * U4: Third order upwinding (Can't "
                        "do 4th order yet).\n * C4: Fourth order central",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddy and UpwindStag
  if (dirOption->isSet("UpwindStag")) {
    dirOption->get("UpwindStag", name, "U1");
  } else if (dirOption->isSet("Upwind")) {
    dirOption->get("Upwind", name, "U1");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "U1");
  } else {
    name = "U1";
  }
  if (strcasecmp(name.c_str(), "U1") == 0) {
    default_y_UpwindStagDeriv = DIFF_U1;
    output.write("	     UpwindStag : First order upwinding\n");
  } else if (strcasecmp(name.c_str(), "U2") == 0) {
    default_y_UpwindStagDeriv = DIFF_U2;
    output.write("	     UpwindStag : Second order upwinding\n");
  } else if (strcasecmp(name.c_str(), "C2") == 0) {
    default_y_UpwindStagDeriv = DIFF_C2;
    output.write("	     UpwindStag : Second order central\n");
  } else if (strcasecmp(name.c_str(), "C4") == 0) {
    default_y_UpwindStagDeriv = DIFF_C4;
    output.write("	     UpwindStag : Fourth order central\n");
  } else {
    throw BoutException("Dont't know what diff method to use for UpwindStag (direction "
                        "y, tried to use %s)!\nOptions are:\n * U1: First order "
                        "upwinding\n * U2: Second order upwinding\n * C2: Second order "
                        "central\n * C4: Fourth order central",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddy and Flux
  if (dirOption->isSet("Flux")) {
    dirOption->get("Flux", name, "U1");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "U1");
  } else {
    name = "U1";
  }
  if (strcasecmp(name.c_str(), "U1") == 0) {
    default_y_FluxDeriv = DIFF_U1;
    output.write("	           Flux : First order upwinding\n");
  } else if (strcasecmp(name.c_str(), "C2") == 0) {
    default_y_FluxDeriv = DIFF_C2;
    output.write("	           Flux : Second order central\n");
  } else if (strcasecmp(name.c_str(), "C4") == 0) {
    default_y_FluxDeriv = DIFF_C4;
    output.write("	           Flux : Fourth order central\n");
  } else {
    throw BoutException("Dont't know what diff method to use for Flux (direction y, "
                        "tried to use %s)!\nOptions are:\n * U1: First order upwinding\n "
                        "* C2: Second order central\n * C4: Fourth order central",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddy and FluxStag
  if (dirOption->isSet("FluxStag")) {
    dirOption->get("FluxStag", name, "U1");
  } else if (dirOption->isSet("Flux")) {
    dirOption->get("Flux", name, "U1");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "U1");
  } else {
    name = "U1";
  }
  if (strcasecmp(name.c_str(), "U1") == 0) {
    default_y_FluxStagDeriv = DIFF_U1;
    output.write("	       FluxStag : First order upwinding\n");
  } else {
    throw BoutException("Dont't know what diff method to use for FluxStag (direction y, "
                        "tried to use %s)!\nOptions are:\n * U1: First order upwinding",
                        name.c_str());
  }
  output.write("\tSetting derivatives for direction z:\n");
  dirOption = option->getSection("ddz");

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddz and First
  if (dirOption->isSet("First")) {
    dirOption->get("First", name, "C2");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "C2");
  } else {
    name = "C2";
  }
  if (strcasecmp(name.c_str(), "C2") == 0) {
    default_z_FirstDeriv = DIFF_C2;
    output.write("	          First : Second order central\n");
  } else if (strcasecmp(name.c_str(), "W2") == 0) {
    default_z_FirstDeriv = DIFF_W2;
    output.write("	          First : Second order WENO\n");
  } else if (strcasecmp(name.c_str(), "C4") == 0) {
    default_z_FirstDeriv = DIFF_C4;
    output.write("	          First : Fourth order central\n");
  } else if (strcasecmp(name.c_str(), "S2") == 0) {
    default_z_FirstDeriv = DIFF_S2;
    output.write("	          First : Smoothing 2nd order\n");
  } else {
    throw BoutException("Dont't know what diff method to use for First (direction z, "
                        "tried to use %s)!\nOptions are:\n * C2: Second order central\n "
                        "* W2: Second order WENO\n * C4: Fourth order central\n * S2: "
                        "Smoothing 2nd order",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddz and FirstStag
  if (dirOption->isSet("FirstStag")) {
    dirOption->get("FirstStag", name, "C2");
  } else if (dirOption->isSet("First")) {
    dirOption->get("First", name, "C2");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "C2");
  } else {
    name = "C2";
  }
  if (strcasecmp(name.c_str(), "C2") == 0) {
    default_z_FirstStagDeriv = DIFF_C2;
    output.write("	      FirstStag : Second order central\n");
  } else if (strcasecmp(name.c_str(), "C4") == 0) {
    default_z_FirstStagDeriv = DIFF_C4;
    output.write("	      FirstStag : Fourth order central\n");
  } else {
    throw BoutException("Dont't know what diff method to use for FirstStag (direction z, "
                        "tried to use %s)!\nOptions are:\n * C2: Second order central\n "
                        "* C4: Fourth order central",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddz and Second
  if (dirOption->isSet("Second")) {
    dirOption->get("Second", name, "C2");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "C2");
  } else {
    name = "C2";
  }
  if (strcasecmp(name.c_str(), "C2") == 0) {
    default_z_SecondDeriv = DIFF_C2;
    output.write("	         Second : Second order central\n");
  } else if (strcasecmp(name.c_str(), "C4") == 0) {
    default_z_SecondDeriv = DIFF_C4;
    output.write("	         Second : Fourth order central\n");
  } else {
    throw BoutException("Dont't know what diff method to use for Second (direction z, "
                        "tried to use %s)!\nOptions are:\n * C2: Second order central\n "
                        "* C4: Fourth order central",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddz and SecondStag
  if (dirOption->isSet("SecondStag")) {
    dirOption->get("SecondStag", name, "C2");
  } else if (dirOption->isSet("Second")) {
    dirOption->get("Second", name, "C2");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "C2");
  } else {
    name = "C2";
  }
  if (strcasecmp(name.c_str(), "C2") == 0) {
    default_z_SecondStagDeriv = DIFF_C2;
    output.write("	     SecondStag : Second order central\n");
  } else {
    throw BoutException("Dont't know what diff method to use for SecondStag (direction "
                        "z, tried to use %s)!\nOptions are:\n * C2: Second order central",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddz and Upwind
  if (dirOption->isSet("Upwind")) {
    dirOption->get("Upwind", name, "U1");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "U1");
  } else {
    name = "U1";
  }
  if (strcasecmp(name.c_str(), "U1") == 0) {
    default_z_UpwindDeriv = DIFF_U1;
    output.write("	         Upwind : First order upwinding\n");
  } else if (strcasecmp(name.c_str(), "U2") == 0) {
    default_z_UpwindDeriv = DIFF_U2;
    output.write("	         Upwind : Second order upwinding\n");
  } else if (strcasecmp(name.c_str(), "C2") == 0) {
    default_z_UpwindDeriv = DIFF_C2;
    output.write("	         Upwind : Second order central\n");
  } else if (strcasecmp(name.c_str(), "U3") == 0) {
    default_z_UpwindDeriv = DIFF_U3;
    output.write("	         Upwind : Third order upwinding\n");
  } else if (strcasecmp(name.c_str(), "U4") == 0) {
    default_z_UpwindDeriv = DIFF_U3;
    output.write(
        "	         Upwind : Third order upwinding (Can't do 4th order yet).\n");
  } else if (strcasecmp(name.c_str(), "C4") == 0) {
    default_z_UpwindDeriv = DIFF_C4;
    output.write("	         Upwind : Fourth order central\n");
  } else {
    throw BoutException("Dont't know what diff method to use for Upwind (direction z, "
                        "tried to use %s)!\nOptions are:\n * U1: First order upwinding\n "
                        "* U2: Second order upwinding\n * C2: Second order central\n * "
                        "U3: Third order upwinding\n * U4: Third order upwinding (Can't "
                        "do 4th order yet).\n * C4: Fourth order central",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddz and UpwindStag
  if (dirOption->isSet("UpwindStag")) {
    dirOption->get("UpwindStag", name, "U1");
  } else if (dirOption->isSet("Upwind")) {
    dirOption->get("Upwind", name, "U1");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "U1");
  } else {
    name = "U1";
  }
  if (strcasecmp(name.c_str(), "U1") == 0) {
    default_z_UpwindStagDeriv = DIFF_U1;
    output.write("	     UpwindStag : First order upwinding\n");
  } else if (strcasecmp(name.c_str(), "U2") == 0) {
    default_z_UpwindStagDeriv = DIFF_U2;
    output.write("	     UpwindStag : Second order upwinding\n");
  } else if (strcasecmp(name.c_str(), "C2") == 0) {
    default_z_UpwindStagDeriv = DIFF_C2;
    output.write("	     UpwindStag : Second order central\n");
  } else if (strcasecmp(name.c_str(), "C4") == 0) {
    default_z_UpwindStagDeriv = DIFF_C4;
    output.write("	     UpwindStag : Fourth order central\n");
  } else {
    throw BoutException("Dont't know what diff method to use for UpwindStag (direction "
                        "z, tried to use %s)!\nOptions are:\n * U1: First order "
                        "upwinding\n * U2: Second order upwinding\n * C2: Second order "
                        "central\n * C4: Fourth order central",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddz and Flux
  if (dirOption->isSet("Flux")) {
    dirOption->get("Flux", name, "U1");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "U1");
  } else {
    name = "U1";
  }
  if (strcasecmp(name.c_str(), "U1") == 0) {
    default_z_FluxDeriv = DIFF_U1;
    output.write("	           Flux : First order upwinding\n");
  } else if (strcasecmp(name.c_str(), "C2") == 0) {
    default_z_FluxDeriv = DIFF_C2;
    output.write("	           Flux : Second order central\n");
  } else if (strcasecmp(name.c_str(), "C4") == 0) {
    default_z_FluxDeriv = DIFF_C4;
    output.write("	           Flux : Fourth order central\n");
  } else {
    throw BoutException("Dont't know what diff method to use for Flux (direction z, "
                        "tried to use %s)!\nOptions are:\n * U1: First order upwinding\n "
                        "* C2: Second order central\n * C4: Fourth order central",
                        name.c_str());
  }

  // This file is auto-generated - do not edit!
  // Setting derivatives for ddz and FluxStag
  if (dirOption->isSet("FluxStag")) {
    dirOption->get("FluxStag", name, "U1");
  } else if (dirOption->isSet("Flux")) {
    dirOption->get("Flux", name, "U1");
  } else if (dirOption->isSet("all")) {
    dirOption->get("all", name, "U1");
  } else {
    name = "U1";
  }
  if (strcasecmp(name.c_str(), "U1") == 0) {
    default_z_FluxStagDeriv = DIFF_U1;
    output.write("	       FluxStag : First order upwinding\n");
  } else {
    throw BoutException("Dont't know what diff method to use for FluxStag (direction z, "
                        "tried to use %s)!\nOptions are:\n * U1: First order upwinding",
                        name.c_str());
  }
}
