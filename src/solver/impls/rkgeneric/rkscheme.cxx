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
  ///These arrays are allocated in the derived class, should
  ///we really free them there as well?
  
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
  steps = rmatrix(getStageCount(),nlocal);
};

//Get the time at given stage
BoutReal RKScheme::setCurTime(const BoutReal timeIn, const BoutReal dt, const int curStage){
  return timeIn+dt*timeCoeffs[curStage];
};

//Get the state vector at given stage
void RKScheme::setCurState(const BoutReal *start, BoutReal *out, const int nlocal, 
			   const int curStage, const BoutReal dt){

  //Set the initial stage
  for(int i=0;i<nlocal;i++){
    out[i] = start[i];
  }

  //If on the first stage we don't need to modify the state
  if(curStage==0) return;//Don't realy need this as below loop won't execute
  if(curStage==1) return;//This shouldn't be here but solution blows up without it at the moment.
  
  //Construct the current state from previous results -- This is expensive
  for(int j=0;j<curStage;j++){
    BoutReal fac=stageCoeffs[curStage][j]*dt;
    for(int i=0;i<nlocal;i++){
      out[i] = out[i] + fac*steps[j][i];
    };
  };
};

//Construct the system state at the next time
void RKScheme::setOutputStates(const BoutReal *start, BoutReal *resultFollow, 
			       BoutReal *resultAlt, const int nlocal, const BoutReal dt){
  //Only really need resultAlt in order to calculate the error, so if not adaptive could avoid it
  //*and* techinically we can write resultFollow-resultAlt in terms of resultCoeffs and steps.

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
  for(int curStage=0;curStage<getStageCount();curStage++){
    for(int i=0;i<nlocal;i++){
      resultFollow[i]=resultFollow[i]+dt*resultCoeffs[curStage][followInd]*steps[curStage][i];
      resultAlt[i]=resultAlt[i]+dt*resultCoeffs[curStage][altInd]*steps[curStage][i];
    }
  }
}

BoutReal RKScheme::updateTimestep(const BoutReal dt, const BoutReal rtol, const BoutReal err){
  BoutReal dtNew;
  dtNew = dt*pow(rtol/(2.0*err),1.0/(order+1.0));
  return dtNew;
}

////////////////////
// PRIVATE
////////////////////
