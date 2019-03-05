#include "rkf34.hxx"

RKF34Scheme::RKF34Scheme(Options *options):RKScheme(options){
  //Set characteristics of scheme
  numStages = 5;
  numOrders = 2;
  order = 3;
  label = "rkf34";
  followHighOrder = false;//true;

  OPTION(options, followHighOrder, followHighOrder);
  if(followHighOrder){
    dtfac = 0.9;
  }else{
    dtfac = 0.9; //Could make this different if needed
  }

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
  stageCoeffs(1, 0) = 1.0 / 4.0;
  //Level 2
  stageCoeffs(2, 0) = 4.0 / 81.0;
  stageCoeffs(2, 1) = 32.0 / 81.0;
  //Level 3
  stageCoeffs(3, 0) = 57.0 / 98.0;
  stageCoeffs(3, 1) = -432.0 / 343.0;
  stageCoeffs(3, 2) = 1053.0 / 686.0;
  //Level 4
  stageCoeffs(4, 0) = 1.0 / 6.0;
  stageCoeffs(4, 1) = 0.0;
  stageCoeffs(4, 2) = 27.0 / 52.0;
  stageCoeffs(4, 3) = 49.0 / 156.0;

  //////////////////////////////////
  //Set coefficients : resultCoeffs
  //////////////////////////////////
  //Level 0
  resultCoeffs(0, 0) = 43.0 / 288.0;
  resultCoeffs(0, 1) = 1.0 / 6.0;
  //Level 1
  resultCoeffs(1, 0) = 0.0;
  resultCoeffs(1, 1) = 0.0;
  //Level 2
  resultCoeffs(2, 0) = 243.0 / 416.0;
  resultCoeffs(2, 1) = 27.0 / 52.0;
  //Level 3
  resultCoeffs(3, 0) = 343.0 / 1872.0;
  resultCoeffs(3, 1) = 49.0 / 156.0;
  //Level 4
  resultCoeffs(4, 0) = 1.0 / 12.0;
  resultCoeffs(4, 1) = 0.0;

  //////////////////////////////////
  //Set coefficients : timeCoeffs
  //////////////////////////////////
  //Level 0
  timeCoeffs[0] = 0.0;
  //Level 1
  timeCoeffs[1] = 1.0/4.0;
  //Level 2
  timeCoeffs[2] = 4.0/9.0;
  //Level 3
  timeCoeffs[3] = 6.0/7.0;
  //Level 4
  timeCoeffs[4] = 1.0;

}

RKF34Scheme::~RKF34Scheme(){
  //Do my cleanup
  
}

