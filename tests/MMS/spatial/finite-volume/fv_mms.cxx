#include "bout/bout.hxx"
#include "bout/difops.hxx"
#include "bout/field.hxx"
#include "bout/field3d.hxx"
#include "bout/field_factory.hxx"
#include "bout/fv_ops.hxx"
#include "bout/globals.hxx"
#include "bout/options.hxx"
#include "bout/options_io.hxx"

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
  // Small amplitude profile
  Field3D input_small{
      FieldFactory::get()->create3D("input_field_small", Options::getRoot(), mesh)};
  Field3D v_small{FieldFactory::get()->create3D("v_small", Options::getRoot(), mesh)};
  // Velocity profile with stagnation points at 0 and pi
  Field3D v_stagnation{
      FieldFactory::get()->create3D("v_stagnation", Options::getRoot(), mesh)};

  Field3D wave{FieldFactory::get()->create3D("wave_speed", Options::getRoot(), mesh)};

  // Communicate to calculate parallel transform.
  mesh->communicate(input, v, v_stagnation, wave, input_small, v_small);

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

  fv_op_test("FV_Div_par", dump, input, FV::Div_par<FV::VanAlbada>(input, v, v),
             "_VanAlbada");
  fv_op_test("FV_Div_par_mod", dump, input,
             FV::Div_par_mod<FV::VanAlbada>(input, v, v, flow_ylow), "_VanAlbada");
  fv_op_test("FV_Div_par_fvv", dump, input, FV::Div_par_fvv<FV::VanAlbada>(input, v, v),
             "_VanAlbada");
  fv_op_test("FV_Div_par_mod_stagnation", dump, input,
             FV::Div_par_mod<FV::VanAlbada>(input, v_stagnation, wave, flow_ylow),
             "_VanAlbada");
  fv_op_test("FV_Div_par_fvv_stagnation", dump, input,
             FV::Div_par_fvv<FV::VanAlbada>(input, v_stagnation, wave), "_VanAlbada");

  fv_op_test("FV_Div_par", dump, input, FV::Div_par<FV::VanAlbadaSmooth>(input, v, v),
             "_VanAlbadaSmooth");
  fv_op_test("FV_Div_par_mod", dump, input,
             FV::Div_par_mod<FV::VanAlbadaSmooth>(input, v, v, flow_ylow),
             "_VanAlbadaSmooth");
  fv_op_test("FV_Div_par_fvv", dump, input,
             FV::Div_par_fvv<FV::VanAlbadaSmooth>(input, v, v), "_VanAlbadaSmooth");
  fv_op_test("FV_Div_par_mod_stagnation", dump, input,
             FV::Div_par_mod<FV::VanAlbadaSmooth>(input, v_stagnation, wave, flow_ylow),
             "_VanAlbadaSmooth");
  fv_op_test("FV_Div_par_fvv_stagnation", dump, input,
             FV::Div_par_fvv<FV::VanAlbadaSmooth>(input, v_stagnation, wave),
             "_VanAlbadaSmooth");

  // Small amplitude: VanAlbadaSmooth should match the central slope (Fromm),
  // while VanAlbada stays limited at extrema
  fv_op_test("FV_Div_par_small", dump, input_small,
             FV::Div_par<FV::Fromm>(input_small, v_small, v_small), "_Fromm");
  fv_op_test("FV_Div_par_mod_small", dump, input_small,
             FV::Div_par_mod<FV::Fromm>(input_small, v_small, v_small, flow_ylow),
             "_Fromm");
  fv_op_test("FV_Div_par_fvv_small", dump, input_small,
             FV::Div_par_fvv<FV::Fromm>(input_small, v_small, v_small), "_Fromm");

  fv_op_test("FV_Div_par_small", dump, input_small,
             FV::Div_par<FV::VanAlbada>(input_small, v_small, v_small), "_VanAlbada");
  fv_op_test("FV_Div_par_mod_small", dump, input_small,
             FV::Div_par_mod<FV::VanAlbada>(input_small, v_small, v_small, flow_ylow),
             "_VanAlbada");
  fv_op_test("FV_Div_par_fvv_small", dump, input_small,
             FV::Div_par_fvv<FV::VanAlbada>(input_small, v_small, v_small), "_VanAlbada");

  fv_op_test("FV_Div_par_small", dump, input_small,
             FV::Div_par<FV::VanAlbadaSmooth>(input_small, v_small, v_small),
             "_VanAlbadaSmooth");
  fv_op_test(
      "FV_Div_par_mod_small", dump, input_small,
      FV::Div_par_mod<FV::VanAlbadaSmooth>(input_small, v_small, v_small, flow_ylow),
      "_VanAlbadaSmooth");
  fv_op_test("FV_Div_par_fvv_small", dump, input_small,
             FV::Div_par_fvv<FV::VanAlbadaSmooth>(input_small, v_small, v_small),
             "_VanAlbadaSmooth");

  fv_op_test("FV_Div_par_K_Grad_par", dump, input, FV::Div_par_K_Grad_par(v, input));
  fv_op_test("FV_Div_par_K_Grad_par_mod", dump, input,
             Div_par_K_Grad_par_mod(v, input, flow_ylow));

  bout::writeDefaultOutputFile(dump);

  BoutFinalise();
}
