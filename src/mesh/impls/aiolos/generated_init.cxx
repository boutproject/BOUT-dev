
// This file is auto-generated - do not edit!
DIFF_METHOD default_stencil[8][3];
enum AIOLOS_DIFF_TYPE {
  AIOLOS_First = 0,
  AIOLOS_FirstStag = 1,
  AIOLOS_Second = 2,
  AIOLOS_SecondStag = 3,
  AIOLOS_Upwind = 4,
  AIOLOS_UpwindStag = 5,
  AIOLOS_Flux = 6,
  AIOLOS_FluxStag = 7,
};

// This file is auto-generated - do not edit!

struct available_stencils {
  const char *key;
  const char *desc;
  DIFF_METHOD method;
};

struct stencils_to_check {
  int id;
  const char *desc;
  const char *default_;
  std::vector<available_stencils> available;
  const char *error;
  std::vector<const char *> option_names;
};

void AiolosMesh::derivs_init(Options *option) {
  std::string name;
  Options *dirOption;
  Options *defOption = option->getSection("diff");
  for (int di : {0, 1, 2}) {
    const char *dds, *d_str;
    bool found;

    if (di == 0) {
      dds = "ddx";
      d_str = "x";
    }
    if (di == 1) {
      dds = "ddy";
      d_str = "y";
    }
    if (di == 2) {
      dds = "ddz";
      d_str = "z";
    }
    output_info.write("\tSetting derivatives for direction %s:\n", d_str);
    dirOption = option->getSection(dds);

    std::vector<stencils_to_check> diff_types{
        // First
        {
            // id
            0,
            // desc
            "First",
            // default
            "C2",
            // list of all available stencils
            {
                {"C2", "Second order central", DIFF_C2},
                {"W2", "Second order WENO", DIFF_W2},
                {"C4", "Fourth order central", DIFF_C4},
                {"S2", "Smoothing 2nd order", DIFF_S2},
            },
            // string for error
            "\n * C2: Second order central\n * W2: Second order WENO\n * C4: Fourth "
            "order central\n * S2: Smoothing 2nd order",
            // list of names to check
            {"First", "all"},
        },
        // FirstStag
        {
            // id
            1,
            // desc
            "FirstStag",
            // default
            "C2",
            // list of all available stencils
            {
                {"C2", "Second order central", DIFF_C2},
                {"C4", "Fourth order central", DIFF_C4},
            },
            // string for error
            "\n * C2: Second order central\n * C4: Fourth order central",
            // list of names to check
            {"FirstStag", "First", "all"},
        },
        // Second
        {
            // id
            2,
            // desc
            "Second",
            // default
            "C2",
            // list of all available stencils
            {
                {"C2", "Second order central", DIFF_C2},
                {"C4", "Fourth order central", DIFF_C4},
            },
            // string for error
            "\n * C2: Second order central\n * C4: Fourth order central",
            // list of names to check
            {"Second", "all"},
        },
        // SecondStag
        {
            // id
            3,
            // desc
            "SecondStag",
            // default
            "C2",
            // list of all available stencils
            {
                {"C2", "Second order central", DIFF_C2},
            },
            // string for error
            "\n * C2: Second order central",
            // list of names to check
            {"SecondStag", "Second", "all"},
        },
        // Upwind
        {
            // id
            4,
            // desc
            "Upwind",
            // default
            "U1",
            // list of all available stencils
            {
                {"U1", "First order upwinding", DIFF_U1},
                {"U2", "Second order upwinding", DIFF_U2},
                {"C2", "Second order central", DIFF_C2},
                {"U3", "Third order upwinding", DIFF_U3},
                {"U4", "Third order upwinding (Can't do 4th order yet).", DIFF_U3},
                {"C4", "Fourth order central", DIFF_C4},
            },
            // string for error
            "\n * U1: First order upwinding\n * U2: Second order upwinding\n * C2: "
            "Second order central\n * U3: Third order upwinding\n * U4: Third order "
            "upwinding (Can't do 4th order yet).\n * C4: Fourth order central",
            // list of names to check
            {"Upwind", "all"},
        },
        // UpwindStag
        {
            // id
            5,
            // desc
            "UpwindStag",
            // default
            "U1",
            // list of all available stencils
            {
                {"U1", "First order upwinding", DIFF_U1},
                {"U2", "Second order upwinding", DIFF_U2},
                {"C2", "Second order central", DIFF_C2},
                {"C4", "Fourth order central", DIFF_C4},
            },
            // string for error
            "\n * U1: First order upwinding\n * U2: Second order upwinding\n * C2: "
            "Second order central\n * C4: Fourth order central",
            // list of names to check
            {"UpwindStag", "Upwind", "all"},
        },
        // Flux
        {
            // id
            6,
            // desc
            "Flux",
            // default
            "U1",
            // list of all available stencils
            {
                {"U1", "First order upwinding", DIFF_U1},
                {"C2", "Second order central", DIFF_C2},
                {"C4", "Fourth order central", DIFF_C4},
            },
            // string for error
            "\n * U1: First order upwinding\n * C2: Second order central\n * C4: Fourth "
            "order central",
            // list of names to check
            {"Flux", "all"},
        },
        // FluxStag
        {
            // id
            7,
            // desc
            "FluxStag",
            // default
            "U1",
            // list of all available stencils
            {
                {"U1", "First order upwinding", DIFF_U1},
            },
            // string for error
            "\n * U1: First order upwinding",
            // list of names to check
            {"FluxStag", "Flux", "all"},
        },
    };

    // This file is auto-generated - do not edit!
    for (const auto &diff_type : diff_types) {
      output_debug.write("Setting derivative for %s and %s", dds, diff_type.desc);
      name = diff_type.default_;
      found = false;
      for (auto opt : {dirOption, defOption}) {
        for (const auto &strf : diff_type.option_names) {
          if (opt->isSet(strf)) {
            opt->get(strf, name, diff_type.default_);
            found = true;
            break;
          }
        }
        if (found)
          break;
      }
      found = false;
      for (const auto &stencil : diff_type.available) {
        if (!found) {
          if (strcasecmp(name.c_str(), stencil.key) == 0) {
            default_stencil[diff_type.id][di] = stencil.method;
            output_info.write("\t%15s : %s\n", stencil.key, stencil.desc);
            found = true;
          }
        }
      }
      if (!found) {
        throw BoutException("Dont't know what diff method to use for %s (direction %s, "
                            "tried to use %s)!\nOptions are:%s",
                            diff_type.desc, d_str, name.c_str(), diff_type.error);
      }
    }
  }
}
