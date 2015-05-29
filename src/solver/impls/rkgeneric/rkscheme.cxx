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
void RKScheme::init(const int nlocalIn, const int neqIn, const BoutReal atolIn, const BoutReal rtolIn){
  //Store configuration data
  nlocal = nlocalIn;
  neq = neqIn;
  atol = atolIn;
  rtol = rtolIn;
  
  //Allocate storage for stages
  steps = rmatrix(getStageCount(),nlocal);
  zeroSteps();

  //Will probably only want the following when debugging, but leave it on for now
  verifyCoeffs();
  printButcherTableau();
};

//Get the time at given stage
BoutReal RKScheme::setCurTime(const BoutReal timeIn, const BoutReal dt, const int curStage){
  return timeIn+dt*timeCoeffs[curStage];
};

//Get the state vector at given stage
void RKScheme::setCurState(const BoutReal *start, BoutReal *out, const int curStage, 
			   const BoutReal dt){

  //Set the initial stage
  for(int i=0;i<nlocal;i++){
    out[i] = start[i];
  }

  //If on the first stage we don't need to modify the state
  if(curStage==0) return;//Don't realy need this as below loop won't execute
  //if(curStage==1) return;//This shouldn't be here but solution blows up without it at the moment.
  
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
			       BoutReal *resultAlt, const BoutReal dt){
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

BoutReal RKScheme::updateTimestep(const BoutReal dt, const BoutReal err){
  BoutReal dtNew;
  // BoutReal stepFac=1.0;
  // if(err<rtol){
  //   dtNew = stepFac*dt*pow(rtol/(2.0*err),1.0/(order+1.0));
  // }else{
  //   dtNew = stepFac*dt*pow(rtol/(2.0*err),1.0/(order));
  // };

  //dtNew = dt*pow(rtol/(2.0*err),1.0/(order+1.0));
  dtNew = dt/pow(2.0*err/rtol,1.0/(order+1.0));

  return dtNew;
}

////////////////////
// PRIVATE
////////////////////

//Check that the coefficients are consistent
void RKScheme::verifyCoeffs(){
  output<<endl;

  //Header
  output<<string(50,'-')<<endl;
  output<<"RK Coefficients consistency check"<<endl;
  output<<string(50,'-')<<endl;

  //Check time and stage coefficients
  output<<setw(10)<<"TimeCoeff"<<" | "<<setw(10)<<"SumStageCoeff"<<endl;
  bool warn=false;
  for(int i=0;i<getStageCount();i++){
    BoutReal tmp=0;
    for(int j=0;j<i;j++){
      tmp+=stageCoeffs[i][j];
    };
    output<<setw(10)<<timeCoeffs[i]<<" | "<<setw(10)<<tmp<<endl;
    if(fabs(timeCoeffs[i]-tmp)>atol) warn=true;
  };

  //Optional warning
  if(warn){
    output<<string(50,'=')<<endl;
    output<<"WARNING: Stage/Time coefficients not consistent"<<endl;
    output<<string(50,'=')<<endl;
    warn=false;
  };

  //Check results coefficients
  output<<string(50,'-')<<endl;
  output<<"Results coefficients (should be 1)"<<endl;
  for(int j=0;j<getNumOrders();j++){
    BoutReal tmp=0;
    for(int i=0;i<getStageCount();i++){
      tmp+=resultCoeffs[i][j];
    };
    output<<"Order : "<<j<<" = "<<tmp<<endl;
    if(fabs(1.0-tmp)>atol) warn=true;
  }

  //Optional warning
  if(warn){
    output<<string(50,'=')<<endl;
    output<<"WARNING: Result coefficients not consistent"<<endl;
    output<<string(50,'=')<<endl;
    warn=false;
  };

  //Footer
  output<<string(50,'-')<<endl;
  output<<endl;
}

void RKScheme::printButcherTableau(){
  int width=15;
  int totalWidth=width*(1+getStageCount())+3;
  output<<endl;
  
  //Header
  output<<string(totalWidth,'-')<<endl;
  output<<"Butcher Tableau"<<endl;
  output<<string(totalWidth,'-')<<endl;

  //Time and stage coeffs
  for(int i=0;i<getStageCount();i++){
    output<<setw(width)<<timeCoeffs[i]<<" | ";
    for(int j=0;j<getStageCount();j++){
      output<<setw(width)<<stageCoeffs[i][j];
    }
    output<<endl;
  };

  //Divider
  output<<string(width,'-')<<" | "<<string(width*getStageCount(),'-')<<endl;

  //Result coeffs for each order
  for(int i=0;i<getNumOrders();i++){
    output<<setw(width)<<i<<" | ";
    for(int j=0;j<getStageCount();j++){
      output<<setw(width)<<resultCoeffs[j][i];
    }
    output<<endl;
  };

  //Footer
  output<<string(totalWidth,'-')<<endl;
  output<<endl;
};

void RKScheme::zeroSteps(){
  for(int i=0;i<getStageCount();i++){
    for(int j=0;j<nlocal;j++){
      steps[i][j]=0.;
    }
  }
};

