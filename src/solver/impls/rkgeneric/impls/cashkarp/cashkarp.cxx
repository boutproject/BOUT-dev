#include "cashkarp.hxx"

CASHKARPScheme::CASHKARPScheme(Options *options):RKScheme(options){
  //Set characteristics of scheme
  numStages = 6;
  numOrders = 2;
  order = 4;
  label = "cashkarp";
  followHighOrder = false;

  OPTION(options, followHighOrder, followHighOrder);

  //Allocate coefficient arrays
  stageCoeffs.reallocate(numStages, numStages);
  resultCoeffs.reallocate(numStages, numOrders);
  timeCoeffs.reallocate(numStages);

  //Zero out arrays (shouldn't be needed, but do for testing)
  for(int i=0;i<numStages;i++){
    timeCoeffs[i]=0.;
    for(int j=0;j<numStages;j++){
      stageCoeffs(i, j) = 0.;
    }
    for(int j=0;j<numOrders;j++){
      resultCoeffs(i, j) = 0.;
    }
  }

  //////////////////////////////////
  //Set coefficients : stageCoeffs
  //////////////////////////////////
  //Level 0
  stageCoeffs(0, 0) = 0.0;
  //Level 1
  stageCoeffs(1, 0) = 1.0 / 5.0;
  //Level 2
  stageCoeffs(2, 0) = 3.0 / 40.0;
  stageCoeffs(2, 1) = 9.0 / 40.0;
  //Level 3
  stageCoeffs(3, 0) = 3.0 / 10.0;
  stageCoeffs(3, 1) = -9.0 / 10.0;
  stageCoeffs(3, 2) = 6.0 / 5.0;
  //Level 4
  stageCoeffs(4, 0) = -11.0 / 54.0;
  stageCoeffs(4, 1) = 5.0 / 2.0;
  stageCoeffs(4, 2) = -70.0 / 27.0;
  stageCoeffs(4, 3) = 35.0 / 27.0;
  //Level 5
  stageCoeffs(5, 0) = 1631.0 / 55296.0;
  stageCoeffs(5, 1) = 175.0 / 512.0;
  stageCoeffs(5, 2) = 575.0 / 13824.0;
  stageCoeffs(5, 3) = 44275.0 / 110592.0;
  stageCoeffs(5, 4) = 253.0 / 4096.0;

  //////////////////////////////////
  //Set coefficients : resultCoeffs
  //////////////////////////////////
  //Level 0
  resultCoeffs(0, 0) = 37.0 / 378.0;
  resultCoeffs(0, 1) = 2825.0 / 27648.0;
  //Level 1
  resultCoeffs(1, 0) = 0.0;
  resultCoeffs(1, 1) = 0.0;
  //Level 2
  resultCoeffs(2, 0) = 250.0 / 621.0;
  resultCoeffs(2, 1) = 18575.0 / 48384.0;
  //Level 3
  resultCoeffs(3, 0) = 125.0 / 594.0;
  resultCoeffs(3, 1) = 13525.0 / 55296.0;
  //Level 4
  resultCoeffs(4, 0) = 0.0;
  resultCoeffs(4, 1) = 277.0 / 14336.0;
  //Level 5
  resultCoeffs(5, 0) = 512.0 / 1771.0;
  resultCoeffs(5, 1) = 1.0 / 4.0;

  //////////////////////////////////
  //Set coefficients : timeCoeffs
  //////////////////////////////////
  //Level 0
  timeCoeffs[0] = 0.0;
  //Level 1
  timeCoeffs[1] = 1.0/5.0;
  //Level 2
  timeCoeffs[2] = 3.0/10.0;
  //Level 3
  timeCoeffs[3] = 3.0/5.0;
  //Level 4
  timeCoeffs[4] = 1.0;
  //Level 5
  timeCoeffs[5] = 7.0/8.0;

}

CASHKARPScheme::~CASHKARPScheme(){
  //Do my cleanup
  
}

