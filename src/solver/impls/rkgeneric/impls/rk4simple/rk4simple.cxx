#include "rk4simple.hxx"

RK4SIMPLEScheme::RK4SIMPLEScheme(Options *options):RKScheme(options){
  //Set characteristics of scheme
  numStages = 12;
  numOrders = 2;
  order = 4;
  label = "rk4simple";
  followHighOrder = true;//false;

  //Allocate coefficient arrays
  stageCoeffs = rmatrix(numStages,numStages);
  resultCoeffs = rmatrix(numStages,numOrders);
  timeCoeffs = new BoutReal[numStages];

  //Zero out arrays (shouldn't be needed, but do for testing)
  for(int i=0;i<numStages;i++){
    timeCoeffs[i]=0.;
    for(int j=0;j<numStages;j++){
      stageCoeffs[i][j]=0.;
    }
    for(int j=0;j<numOrders;j++){
      resultCoeffs[i][j]=0.;
    }
  }

  //////////////////////////////////
  //Set coefficients : stageCoeffs
  //////////////////////////////////
  //Single large step
  stageCoeffs[0][0] = 0.0;
  stageCoeffs[1][0] = 1.0/2.0;
  stageCoeffs[2][1] = 1.0/2.0;
  stageCoeffs[3][2] = 1.0;

  //First half step
  stageCoeffs[4][3] = 0.0;
  stageCoeffs[5][4] = 1.0/4.0;
  stageCoeffs[6][5] = 1.0/4.0;
  stageCoeffs[7][6] = 1.0/2.0;

  //Second half step -- standard
  stageCoeffs[8][7] = 0.0;
  stageCoeffs[9][8] = 1.0/4.0;
  stageCoeffs[10][9] = 1.0/4.0;
  stageCoeffs[11][10] = 1.0/2.0;
  //Second half step -- result from first half
  stageCoeffs[8][4] = 1.0/12.0; stageCoeffs[8][5] = 1.0/6.0;
  stageCoeffs[8][6] = 1.0/6.0; stageCoeffs[8][7] = 1.0/12.0;

  stageCoeffs[9][4] = 1.0/12.0; stageCoeffs[9][5] = 1.0/6.0;
  stageCoeffs[9][6] = 1.0/6.0; stageCoeffs[9][7] = 1.0/12.0;

  stageCoeffs[10][4] = 1.0/12.0; stageCoeffs[10][5] = 1.0/6.0;
  stageCoeffs[10][6] = 1.0/6.0; stageCoeffs[10][7] = 1.0/12.0;

  stageCoeffs[11][4] = 1.0/12.0; stageCoeffs[11][5] = 1.0/6.0;
  stageCoeffs[11][6] = 1.0/6.0; stageCoeffs[11][7] = 1.0/12.0;

  // stageCoeffs[8][3] = 1.0/12.0; stageCoeffs[8][4] = 1.0/6.0;
  // stageCoeffs[8][5] = 1.0/6.0; stageCoeffs[8][6] = 1.0/12.0;

  // stageCoeffs[9][3] = 1.0/12.0; stageCoeffs[9][4] = 1.0/6.0;
  // stageCoeffs[9][5] = 1.0/6.0; stageCoeffs[9][6] = 1.0/12.0;

  // stageCoeffs[10][3] = 1.0/12.0; stageCoeffs[10][4] = 1.0/6.0;
  // stageCoeffs[10][5] = 1.0/6.0; stageCoeffs[10][6] = 1.0/12.0;

  // stageCoeffs[11][3] = 1.0/12.0; stageCoeffs[11][4] = 1.0/6.0;
  // stageCoeffs[11][5] = 1.0/6.0; stageCoeffs[11][6] = 1.0/12.0;


  //////////////////////////////////
  //Set coefficients : resultCoeffs
  //////////////////////////////////
  //Large time step
  resultCoeffs[0][1] = 1.0/6.0;
  resultCoeffs[1][1] = 1.0/3.0;
  resultCoeffs[2][1] = 1.0/3.0;
  resultCoeffs[3][1] = 1.0/6.0;

  //small time step
  resultCoeffs[4][0] = 1.0/12.0;
  resultCoeffs[5][0] = 1.0/6.0;
  resultCoeffs[6][0] = 1.0/6.0;
  resultCoeffs[7][0] = 1.0/12.0;
  //Second half
  resultCoeffs[8][0] = 1.0/12.0;
  resultCoeffs[9][0] = 1.0/6.0;
  resultCoeffs[10][0] = 1.0/6.0;
  resultCoeffs[11][0] = 1.0/12.0;

  //////////////////////////////////
  //Set coefficients : timeCoeffs
  //////////////////////////////////
  //Large step
  timeCoeffs[0] = 0.0;
  timeCoeffs[1] = 1.0/2.0;
  timeCoeffs[2] = 1.0/2.0;
  timeCoeffs[3] = 1.0;

  //First small step
  timeCoeffs[4] = 0.0;
  timeCoeffs[5] = 1.0/4.0;
  timeCoeffs[6] = 1.0/4.0;
  timeCoeffs[7] = 1.0/2.0;

  //Second small step
  timeCoeffs[8] = 1.0/2.0+0.0;
  timeCoeffs[9] = 1.0/2.0+1.0/4.0;
  timeCoeffs[10] = 1.0/2.0+1.0/4.0;
  timeCoeffs[11] = 1.0/2.0+1.0/2.0;

};

RK4SIMPLEScheme::~RK4SIMPLEScheme(){
  //Do my cleanup
  
};

