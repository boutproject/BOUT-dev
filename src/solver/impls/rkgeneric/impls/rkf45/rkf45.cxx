#include "rkf45.hxx"

RKF45Scheme::RKF45Scheme(Options *options):RKScheme(options){
  //Set characteristics of scheme
  numStages = 6;
  numOrders = 2;
  order = 4;
  label = "rkf45";
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
  stageCoeffs(1, 0) = 1.0 / 4.0;
  //Level 2
  stageCoeffs(2, 0) = 3.0 / 32.0;
  stageCoeffs(2, 1) = 9.0 / 32.0;
  //Level 3
  stageCoeffs(3, 0) = 1932.0 / 2197.0;
  stageCoeffs(3, 1) = -7200.0 / 2197.0;
  stageCoeffs(3, 2) = 7296.0 / 2197.0;
  //Level 4
  stageCoeffs(4, 0) = 439.0 / 216.0;
  stageCoeffs(4, 1) = -8.0;
  stageCoeffs(4, 2) = 3680.0 / 513.0;
  stageCoeffs(4, 3) = -845.0 / 4104.0;
  //Level 5
  stageCoeffs(5, 0) = -8.0 / 27.0;
  stageCoeffs(5, 1) = 2.0;
  stageCoeffs(5, 2) = -3544.0 / 2565.0;
  stageCoeffs(5, 3) = 1859.0 / 4104.0;
  stageCoeffs(5, 4) = -11.0 / 40.0;

  //////////////////////////////////
  //Set coefficients : resultCoeffs
  //////////////////////////////////
  //Level 0
  resultCoeffs(0, 0) = 16.0 / 135.0;
  resultCoeffs(0, 1) = 25.0 / 216.0;
  //Level 1
  resultCoeffs(1, 0) = 0.0;
  resultCoeffs(1, 1) = 0.0;
  //Level 2
  resultCoeffs(2, 0) = 6656.0 / 12825.0;
  resultCoeffs(2, 1) = 1408.0 / 2565.0;
  //Level 3
  resultCoeffs(3, 0) = 28561.0 / 56430.0;
  resultCoeffs(3, 1) = 2197.0 / 4104.0;
  //Level 4
  resultCoeffs(4, 0) = -9.0 / 50.0;
  resultCoeffs(4, 1) = -1.0 / 5.0;
  //Level 5
  resultCoeffs(5, 0) = 2.0 / 55.0;
  resultCoeffs(5, 1) = 0.0;

  //////////////////////////////////
  //Set coefficients : timeCoeffs
  //////////////////////////////////
  //Level 0
  timeCoeffs[0] = 0.0;
  //Level 1
  timeCoeffs[1] = 1.0/4.0;
  //Level 2
  timeCoeffs[2] = 3.0/8.0;
  //Level 3
  timeCoeffs[3] = 12.0/13.0;
  //Level 4
  timeCoeffs[4] = 1.0;
  //Level 5
  timeCoeffs[5] = 1.0/2.0;

}

RKF45Scheme::~RKF45Scheme(){
  //Do my cleanup
  
}

