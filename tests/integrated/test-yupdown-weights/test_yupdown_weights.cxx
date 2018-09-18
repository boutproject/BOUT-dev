#include <bout.hxx>
#include <derivs.hxx>
#include "../../../src/mesh/parallel/shiftedmetricinterp.hxx"

// Y derivative using yup() and ydown() fields
const Field3D DDY_yud(const Field3D &f) {
  Field3D result;
  result.allocate();

  result = 0.0;
  
  for(int i=0;i<mesh->LocalNx;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->LocalNz;k++)
        result(i,j,k) = 0.5*(f.yup()(i,j+1,k) - f.ydown()(i,j-1,k));

  return result;
}

// Y derivative constructed from interpolation weights
const Field3D DDY_weights(const Field3D &f) {

  ShiftedMetricInterp* sm = new ShiftedMetricInterp(*mesh);
  Field3D result;
  result.allocate();

  result = 0.0;
  
  for(int i=0;i<mesh->LocalNx;i++){
    for(int j=mesh->ystart;j<=mesh->yend;j++){
      for(int k=0;k<mesh->LocalNz;k++){
	std::vector<Interpolation::positionsAndWeights> pw_up = sm->interp_yup->getWeightsForYApproximation(i,j,k,1);
	std::vector<Interpolation::positionsAndWeights> pw_down = sm->interp_ydown->getWeightsForYApproximation(i,j,k,-1);

	for (const auto &p : pw_up){
	  result(i,j,k) += 0.5*p.weight*f(p.i,p.j,p.k);
	}
	for (const auto &p : pw_down){
	  result(i,j,k) -= 0.5*p.weight*f(p.i,p.j,p.k);
	}
      }
    }
  }

  return result;
}

int main(int argc, char** argv) {

  BoutInitialise(argc, argv);

  // Read variable from mesh
  Field3D var;
  mesh->get(var, "var");
  
  // Var starts in orthogonal X-Z coordinates

  // Calculate yup and ydown
  mesh->communicate(var);

  // Calculate d/dy using yup() and ydown() fields
  Field3D ddy = DDY_yud(var);
  // This is also equal to:
  //Field3D ddy = DDY(var);

  // Calculate d/dy using Hermite spline weights
  Field3D ddy2 = DDY_weights(var);
  
  SAVE_ONCE2(ddy, ddy2);
  dump.write();

  BoutFinalise();

  return 0;
}
