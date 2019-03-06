#include "rk4simple.hxx"

RK4SIMPLEScheme::RK4SIMPLEScheme(Options *options):RKScheme(options){
  //Set characteristics of scheme
  numStages = 11;
  numOrders = 2;
  order = 4;
  label = "rk4simple";
  followHighOrder = true;

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
  //Single large step
  stageCoeffs(0, 0) = 0.0;
  stageCoeffs(1, 0) = 1.0 / 2.0;
  stageCoeffs(2, 1) = 1.0 / 2.0;
  stageCoeffs(3, 2) = 1.0;

  //First half step
  stageCoeffs(4, 0) = 1.0 / 4.0;
  stageCoeffs(5, 4) = 1.0 / 4.0;
  stageCoeffs(6, 5) = 1.0 / 2.0;

  //Second half step -- standard
  stageCoeffs(7, 6) = 0.0;
  stageCoeffs(8, 7) = 1.0 / 4.0;
  stageCoeffs(9, 8) = 1.0 / 4.0;
  stageCoeffs(10, 9) = 1.0 / 2.0;

  //Second half step -- result from first half
  stageCoeffs(7, 0) = 1.0 / 12.0;
  stageCoeffs(7, 4) = 1.0 / 6.0;
  stageCoeffs(7, 5) = 1.0 / 6.0;
  stageCoeffs(7, 6) = 1.0 / 12.0;

  stageCoeffs(8, 0) = 1.0 / 12.0;
  stageCoeffs(8, 4) = 1.0 / 6.0;
  stageCoeffs(8, 5) = 1.0 / 6.0;
  stageCoeffs(8, 6) = 1.0 / 12.0;

  stageCoeffs(9, 0) = 1.0 / 12.0;
  stageCoeffs(9, 4) = 1.0 / 6.0;
  stageCoeffs(9, 5) = 1.0 / 6.0;
  stageCoeffs(9, 6) = 1.0 / 12.0;

  stageCoeffs(10, 0) = 1.0 / 12.0;
  stageCoeffs(10, 4) = 1.0 / 6.0;
  stageCoeffs(10, 5) = 1.0 / 6.0;
  stageCoeffs(10, 6) = 1.0 / 12.0;

  //////////////////////////////////
  //Set coefficients : resultCoeffs
  //////////////////////////////////
  //Large time step
  resultCoeffs(0, 1) = 1.0 / 6.0;
  resultCoeffs(1, 1) = 1.0 / 3.0;
  resultCoeffs(2, 1) = 1.0 / 3.0;
  resultCoeffs(3, 1) = 1.0 / 6.0;

  //small time step
  resultCoeffs(0, 0) = 1.0 / 12.0;
  resultCoeffs(4, 0) = 1.0 / 6.0;
  resultCoeffs(5, 0) = 1.0 / 6.0;
  resultCoeffs(6, 0) = 1.0 / 12.0;
  //Second half
  resultCoeffs(7, 0) = 1.0 / 12.0;
  resultCoeffs(8, 0) = 1.0 / 6.0;
  resultCoeffs(9, 0) = 1.0 / 6.0;
  resultCoeffs(10, 0) = 1.0 / 12.0;

  //////////////////////////////////
  //Set coefficients : timeCoeffs
  //////////////////////////////////
  //Large step
  timeCoeffs[0] = 0.0;
  timeCoeffs[1] = 1.0/2.0;
  timeCoeffs[2] = 1.0/2.0;
  timeCoeffs[3] = 1.0;

  //First small step
  timeCoeffs[4] = 1.0/4.0;
  timeCoeffs[5] = 1.0/4.0;
  timeCoeffs[6] = 1.0/2.0;

  //Second small step
  timeCoeffs[7] = 1.0/2.0+0.0;
  timeCoeffs[8] = 1.0/2.0+1.0/4.0;
  timeCoeffs[9] = 1.0/2.0+1.0/4.0;
  timeCoeffs[10] = 1.0/2.0+1.0/2.0;

}

RK4SIMPLEScheme::~RK4SIMPLEScheme(){
  //Do my cleanup
  
}

BoutReal RK4SIMPLEScheme::setOutputStates(const Array<BoutReal> &start, const BoutReal dt,
                                          Array<BoutReal> &resultFollow) {
  //return RKScheme::setOutputStates(start,dt,resultFollow);
  if(followHighOrder){
    for(int i=0;i<nlocal;i++){
      if(adaptive){
        resultAlt[i] =
            start[i] +
            dt * (resultCoeffs(0, 1) * steps(0, i) + resultCoeffs(1, 1) * steps(1, i) +
                  resultCoeffs(2, 1) * steps(2, i) + resultCoeffs(3, 1) * steps(3, i));
      }

      resultFollow[i] =
          start[i] +
          dt * (resultCoeffs(0, 0) * steps(0, i) + resultCoeffs(4, 0) * steps(4, i) +
                resultCoeffs(5, 0) * steps(5, i) + resultCoeffs(6, 0) * steps(6, i));

      resultFollow[i] =
          resultFollow[i] +
          dt * (resultCoeffs(7, 0) * steps(7, i) + resultCoeffs(8, 0) * steps(8, i) +
                resultCoeffs(9, 0) * steps(9, i) + resultCoeffs(10, 0) * steps(10, i));
    }
  }else{
    for(int i=0;i<nlocal;i++){

      resultFollow[i] =
          start[i] +
          dt * (resultCoeffs(0, 1) * steps(0, i) + resultCoeffs(1, 1) * steps(1, i) +
                resultCoeffs(2, 1) * steps(2, i) + resultCoeffs(3, 1) * steps(3, i));

      if(adaptive) {
        resultAlt[i] =
            start[i] +
            dt * (resultCoeffs(0, 0) * steps(0, i) + resultCoeffs(4, 0) * steps(4, i) +
                  resultCoeffs(5, 0) * steps(5, i) + resultCoeffs(6, 0) * steps(6, i));

        resultAlt[i] =
            resultAlt[i] +
            dt * (resultCoeffs(7, 0) * steps(7, i) + resultCoeffs(8, 0) * steps(8, i) +
                  resultCoeffs(9, 0) * steps(9, i) + resultCoeffs(10, 0) * steps(10, i));
      }
    }
  }

  return getErr(resultFollow, resultAlt);
}
