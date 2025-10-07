#include "bout/bout.hxx"
#include "bout/field_factory.hxx"

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  using bout::globals::mesh;

  Field3D input{FieldFactory::get()->create3D("input_field", Options::getRoot(), mesh)};

  // Communicate to calculate parallel transform
  mesh->communicate(input);

  Options dump;
  // Add mesh geometry variables
  mesh->outputVars(dump);

  auto* factory = FieldFactory::get();
  {
    Field3D solution{factory->create3D("grad_par_solution", Options::getRoot(), mesh)};
    Field3D result{Grad_par(input)};
    Field3D error{result - solution};

    dump["grad_par_l_2"] = sqrt(mean(SQ(error), true, "RGN_NOBNDRY"));
    dump["grad_par_l_inf"] = max(abs(error), true, "RGN_NOBNDRY");

    dump["grad_par_result"] = result;
    dump["grad_par_error"] = error;
    dump["grad_par_input"] = input;
    dump["grad_par_solution"] = solution;

    for (int slice = 1; slice < mesh->ystart; ++slice) {
      dump[fmt::format("grad_par_input.ynext(-{})", slice)] = input.ynext(-slice);
      dump[fmt::format("grad_par_input.ynext({})", slice)] = input.ynext(slice);
    }
  }
  {
    Field3D solution{factory->create3D("grad2_par2_solution", Options::getRoot(), mesh)};
    Field3D result{Grad2_par2(input)};
    Field3D error{result - solution};

    dump["grad2_par2_l_2"] = sqrt(mean(SQ(error), true, "RGN_NOBNDRY"));
    dump["grad2_par2_l_inf"] = max(abs(error), true, "RGN_NOBNDRY");

    dump["grad2_par2_result"] = result;
    dump["grad2_par2_error"] = error;
    dump["grad2_par2_input"] = input;
    dump["grad2_par2_solution"] = solution;

    for (int slice = 1; slice < mesh->ystart; ++slice) {
      dump[fmt::format("grad2_par2_input.ynext(-{})", slice)] = input.ynext(-slice);
      dump[fmt::format("grad2_par2_input.ynext({})", slice)] = input.ynext(slice);
    }
  }

  bout::writeDefaultOutputFile(dump);

  BoutFinalise();
}
