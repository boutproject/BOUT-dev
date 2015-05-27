#include "rkschemefactory.hxx"
#include <bout/rkscheme.hxx>
#include <output.hxx>
#include <cmath>
#include <boutcomm.hxx>

////////////////////
// PUBLIC
////////////////////

//Initialise
RKScheme::RKScheme(Options *opts){
  //Currently not reading anything from the options here

  //Init the pointer arrays to null
  stageCoeffs = (BoutReal**)NULL;
  resultCoeffs = (BoutReal**)NULL;
  timeCoeffs = (BoutReal*)NULL;
  steps = (BoutReal**)NULL;
};

//Cleanup
RKScheme::~RKScheme(){
  //stageCoeffs
  free_rmatrix(stageCoeffs);

  //resultCoeffs
  free_rmatrix(resultCoeffs);

  //steps
  free_rmatrix(steps);

  //timeCoeffs
  delete[] timeCoeffs;
};

//Finish generic initialisation
void RKScheme::init(const int nlocal){
  //Allocate storage for stages
  steps = rmatrix(numStages,nlocal);
};

BoutReal RKScheme::setCurTime(const BoutReal timeIn, const BoutReal dt, const int curStage){
  return timeIn+dt*timeCoeffs[curStage];
};

void RKScheme::setCurState(const BoutReal *start, BoutReal *out, const int nlocal, const int curStage, const BoutReal dt){
  for(int i=0;i<nlocal;i++){
    out[i] = start[i];
    if(curStage==1) continue;
    for(int j=0;j<curStage;j++){
      out[i] = out[i] + stageCoeffs[curStage][j]*dt*steps[j][i];
    };
  }
};

void RKScheme::setOutputStates(const BoutReal *start, BoutReal *resultFollow, BoutReal *resultAlt, const int nlocal, const BoutReal dt){

  int followInd, altInd;
  if(followHighOrder){
    followInd=0; altInd=1;
  }else{
    followInd=1; altInd=0;
  }

  //Initialise the return data
  for(int i=0;i<nlocal;i++){
    resultFollow[i]=start[i];
    resultAlt[i]=start[i];
  }

  //Now construct the two solutions
  for(int curStage=0;curStage<numStages;curStage++){
    for(int i=0;i<nlocal;i++){
      resultFollow[i]=resultFollow[i]+dt*resultCoeffs[curStage][followInd]*steps[curStage][i];
      resultAlt[i]=resultAlt[i]+dt*resultCoeffs[curStage][altInd]*steps[curStage][i];
    }
  }
}

////////////////////
// PRIVATE
////////////////////

//Print out details of the scheme
void RKScheme::debugPrint(){
};

void RKScheme::setStageCoeffs(int stage, const vector<BoutReal> coeffs){
};

void RKScheme::setResultCoeffs(int stage, const BoutReal highOrder, const BoutReal lowOrder){
};

void RKScheme::setTimeCoeffs(int stage, const BoutReal coeff){
};
