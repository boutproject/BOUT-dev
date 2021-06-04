#include <bout.hxx>
#include <derivs.hxx>
#include "../../../src/mesh/parallel/shiftedmetricinterp.hxx"

// Y derivative using yup() and ydown() fields
Field3D DDY_yud(const Field3D &f) {
  Field3D result{0.0};
  const auto* mesh = f.getMesh();

  for(int i=0;i<mesh->LocalNx;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->LocalNz;k++)
        result(i,j,k) = 0.5*(f.yup()(i,j+1,k) - f.ydown()(i,j-1,k));

  return result;
}

// Y derivative constructed from interpolation weights
Field3D DDY_weights(const Field3D &f) {

  ParallelTransform& pt = f.getCoordinates()->getParallelTransform();
  Field3D result{0.0};
  const auto* mesh = f.getMesh();

  for(int i=mesh->xstart;i<=mesh->xend;i++){
    for(int j=mesh->ystart;j<=mesh->yend;j++){
      for(int k=0;k<mesh->LocalNz;k++){
        std::vector<ParallelTransform::PositionsAndWeights> pw_up = pt.getWeightsForYUpApproximation(i,j,k);
        std::vector<ParallelTransform::PositionsAndWeights> pw_down = pt.getWeightsForYDownApproximation(i,j,k);

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
  using bout::globals::mesh;

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
  
  Options::root()["ddy"] = ddy;
  Options::root()["ddy2"] = ddy2;
  bout::writeDefaultOutputFile();

  BoutFinalise();

  return 0;
}
