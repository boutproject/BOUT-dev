#include <bout.hxx>
#include <boutmain.hxx>

#include <bout/region.hxx>
#include <bout/assert.hxx>

Field3D input, expectedDdz, expectedD2dz2, actualDdz, actualD2dz2;
BoutReal errorDdz, errorD2dz2;

int physics_init(bool restarting) {
  output.enable(false);
  SOLVE_FOR(input);
  SOLVE_FOR(expectedDdz);
  SOLVE_FOR(expectedD2dz2);
  SAVE_REPEAT2(errorDdz,actualDdz);
  SAVE_REPEAT2(errorD2dz2,actualD2dz2);
  actualDdz = 0.;
  actualD2dz2 = 0.;
  return 0;
}

int physics_run(BoutReal t) {
  auto region = mesh->getRegion("RGN_ALL");
  auto regionZp = region; regionZp.periodicShift(1, mesh->LocalNz);
  auto zpInd = regionZp.getIndices();
  auto regionZm = region; regionZm.periodicShift(-1, mesh->LocalNz);
  auto zmInd = regionZm.getIndices();

  auto dzInv = 1.0/mesh->coordinates()->dz;
  BLOCK_REGION_LOOP_COUNTER(region, index, counter,
			    actualD2dz2[index] =
			    dzInv*dzInv*(input[zmInd[counter]] - 2*input[index] + input[zpInd[counter]]);
			    actualDdz[index] =
			    0.5*dzInv*(input[zpInd[counter]] - input[zmInd[counter]]);
			    );

  errorDdz = max(abs(expectedDdz-actualDdz), true);
  errorD2dz2 = max(abs(expectedD2dz2-actualD2dz2), true);
  output.enable(true);
  output<<errorDdz<<" "<<errorD2dz2<<endl;
  output.enable(false);

  ddt(input) = 0.;
  ddt(expectedDdz) = 0.;
  ddt(expectedD2dz2) = 0.;
  return 0;
}
