#include "bout/bout.hxx"
#include "bout/field.hxx"
#include "bout/field3d.hxx"
#include "bout/field_factory.hxx"
#include "bout/fv_ops.hxx"
#include "bout/globals.hxx"
#include "bout/options.hxx"
#include "bout/options_io.hxx"
#include "bout/utils.hxx"

#include <fmt/format.h>

#include <cmath>
#include <string>

namespace {
auto fv_op_test(const std::string& name, Options& dump, const Field3D& input,
                const Field3D& result, std::string suffix = "") {
  auto* mesh = input.getMesh();
  const Field3D solution{FieldFactory::get()->create3D(fmt::format("{}_solution", name),
                                                       Options::getRoot(), mesh)};
  const Field3D error{result - solution};

  dump[fmt::format("{}{}_l_2", name, suffix)] =
      sqrt(mean(SQ(error), true, "RGN_NOBNDRY"));
  dump[fmt::format("{}{}_l_inf", name, suffix)] = max(abs(error), true, "RGN_NOBNDRY");

  dump[fmt::format("{}{}_result", name, suffix)] = result;
  dump[fmt::format("{}{}_error", name, suffix)] = error;
  dump[fmt::format("{}{}_input", name, suffix)] = input;
  dump[fmt::format("{}{}_solution", name, suffix)] = solution;
}
} // namespace

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  using bout::globals::mesh;

  Field3D input{FieldFactory::get()->create3D("input_field", Options::getRoot(), mesh)};
  Field3D v{FieldFactory::get()->create3D("v", Options::getRoot(), mesh)};

  // Communicate to calculate parallel transform.
  mesh->communicate(input, v);

  Options dump;
  // Add mesh geometry variables
  mesh->outputVars(dump);
  dump["v"] = v;

  // Dummy variable for *_mod overloads
  Field3D flow_ylow;

  fv_op_test("FV_Div_par", dump, input, FV::Div_par<FV::MC>(input, v, v), "_MC");
  fv_op_test("FV_Div_par_mod", dump, input,
             FV::Div_par_mod<FV::MC>(input, v, v, flow_ylow), "_MC");
  fv_op_test("FV_Div_par_fvv", dump, input, FV::Div_par_fvv<FV::MC>(input, v, v), "_MC");

  fv_op_test("FV_Div_par", dump, input, FV::Div_par<FV::Upwind>(input, v, v), "_Upwind");
  fv_op_test("FV_Div_par_mod", dump, input,
             FV::Div_par_mod<FV::Upwind>(input, v, v, flow_ylow), "_Upwind");
  fv_op_test("FV_Div_par_fvv", dump, input, FV::Div_par_fvv<FV::Upwind>(input, v, v),
             "_Upwind");

  fv_op_test("FV_Div_par", dump, input, FV::Div_par<FV::Fromm>(input, v, v), "_Fromm");
  fv_op_test("FV_Div_par_mod", dump, input,
             FV::Div_par_mod<FV::Fromm>(input, v, v, flow_ylow), "_Fromm");
  fv_op_test("FV_Div_par_fvv", dump, input, FV::Div_par_fvv<FV::Fromm>(input, v, v),
             "_Fromm");

  fv_op_test("FV_Div_par", dump, input, FV::Div_par<FV::MinMod>(input, v, v), "_MinMod");
  fv_op_test("FV_Div_par_mod", dump, input,
             FV::Div_par_mod<FV::MinMod>(input, v, v, flow_ylow), "_MinMod");
  fv_op_test("FV_Div_par_fvv", dump, input, FV::Div_par_fvv<FV::MinMod>(input, v, v),
             "_MinMod");

  fv_op_test("FV_Div_par", dump, input, FV::Div_par<FV::Superbee>(input, v, v),
             "_Superbee");
  fv_op_test("FV_Div_par_mod", dump, input,
             FV::Div_par_mod<FV::Superbee>(input, v, v, flow_ylow), "_Superbee");
  fv_op_test("FV_Div_par_fvv", dump, input, FV::Div_par_fvv<FV::Superbee>(input, v, v),
             "_Superbee");

  fv_op_test("FV_Div_par_K_Grad_par", dump, input, FV::Div_par_K_Grad_par(v, input));
  fv_op_test("FV_Div_par_K_Grad_par_mod", dump, input,
             Div_par_K_Grad_par_mod(v, input, flow_ylow));

  bout::writeDefaultOutputFile(dump);

  BoutFinalise();
}
