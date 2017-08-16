#include "rkschemefactory.hxx"
#include <bout/rkscheme.hxx>
#include <output.hxx>
#include <cmath>
#include <boutcomm.hxx>
#include "unused.hxx"

////////////////////
// PUBLIC
////////////////////

//Initialise
RKScheme::RKScheme(Options *UNUSED(opts))
    : steps(nullptr), stageCoeffs(nullptr), resultCoeffs(nullptr), timeCoeffs(nullptr),
      resultAlt(nullptr) {
  // Currently not reading anything from the options here

  // Initialise internals
  dtfac = 1.0; // Time step factor
}

//Cleanup
RKScheme::~RKScheme(){
  ///These arrays are allocated in the derived class, should
  ///we really free them there as well?
  
  //stageCoeffs
  if(stageCoeffs != nullptr) free_matrix(stageCoeffs);

  //resultCoeffs
  if(stageCoeffs != nullptr) free_matrix(resultCoeffs);

  //steps
  if(stageCoeffs != nullptr) free_matrix(steps);

  //timeCoeffs
  if(stageCoeffs != nullptr) delete[] timeCoeffs;
  
  if(stageCoeffs != nullptr) delete[] resultAlt;
}

//Finish generic initialisation
void RKScheme::init(const int nlocalIn, const int neqIn, const bool adaptiveIn, const BoutReal atolIn, 
		    const BoutReal rtolIn, Options *options){

  bool diagnose;
  OPTION(options, dtfac, dtfac); //Time step adjustment factor
  OPTION(options, diagnose, false); //Diagnostics enabled?

  //Store configuration data
  nlocal = nlocalIn;
  neq = neqIn;
  atol = atolIn;
  rtol = rtolIn;
  adaptive = adaptiveIn;

  //Allocate storage for stages
  steps = matrix<BoutReal>(getStageCount(),nlocal);
  zeroSteps();

  //Allocate array for storing alternative order result
  if(adaptive) resultAlt = new BoutReal[nlocal]; //Result--alternative order

  //Will probably only want the following when debugging, but leave it on for now
  if(diagnose){
    verifyCoeffs();
    printButcherTableau();
  }
}

//Get the time at given stage
BoutReal RKScheme::setCurTime(const BoutReal timeIn, const BoutReal dt, const int curStage){
  return timeIn+dt*timeCoeffs[curStage];
}

//Get the state vector at given stage
void RKScheme::setCurState(const BoutReal *start, BoutReal *out, const int curStage, 
			   const BoutReal dt){

  //Set the initial stage
  for(int i=0;i<nlocal;i++){
    out[i] = start[i];
  }

  //If on the first stage we don't need to modify the state
  if(curStage==0) return;//Don't realy need this as below loop won't execute
  
  //Construct the current state from previous results -- This is expensive
  for(int j=0;j<curStage;j++){
    if(abs(stageCoeffs[curStage][j]) < atol) continue;
    BoutReal fac=stageCoeffs[curStage][j]*dt;
    for(int i=0;i<nlocal;i++){
      out[i] = out[i] + fac*steps[j][i];
    }
  }
}

//Construct the system state at the next time
BoutReal RKScheme::setOutputStates(const BoutReal *start, const BoutReal dt, BoutReal *resultFollow){
  //Only really need resultAlt in order to calculate the error, so if not adaptive could avoid it
  //*and* technically we can write resultFollow-resultAlt in terms of resultCoeffs and steps.

  int followInd, altInd;
  if(followHighOrder){
    followInd=0; altInd=1;
  }else{
    followInd=1; altInd=0;
  }

  //NOTE: It's slightly slower to split construction of output
  //      into two separate loops as we're doing below *when*
  //      we want both solutions, but it's faster when not adaptive.
  //      To get the best of both worlds we could probably do something
  //      like:
  /*      
	  if(adapative){
	    constructOutputs;
	  }else{
            constructOutput;
	  }
   */

  //Get the result
  constructOutput(start,dt,followInd,resultFollow);

  //If adaptive get the second state
  if(adaptive){
    constructOutput(start,dt,altInd,resultAlt);
  }

  //Get the error coefficient
  return getErr(resultFollow,resultAlt);
}

BoutReal RKScheme::updateTimestep(const BoutReal dt, const BoutReal err){
  return dtfac*dt*pow(rtol/(2.0*err),1.0/(order+1.0));
}

////////////////////
// PRIVATE
////////////////////

//Estimate the error, given two solutions
BoutReal RKScheme::getErr(BoutReal *solA, BoutReal *solB){
  BoutReal err=0.;

  //If not adaptive don't care about the error
  if(!adaptive){return err;}

  //Get local part of relative error
  BoutReal local_err = 0.;
  for(int i=0;i<nlocal;i++) {
    local_err += fabs(solA[i] - solB[i]) / ( fabs(solA[i]) + fabs(solB[i]) + atol );
  }
  //Reduce over procs
  if(MPI_Allreduce(&local_err, &err, 1, MPI_DOUBLE, MPI_SUM, BoutComm::get())) {
    throw BoutException("MPI_Allreduce failed");
  }
  //Normalise by number of values
  err /= (BoutReal) neq;
  
  return err;
}

void RKScheme::constructOutput(const BoutReal *start, const BoutReal dt, 
			       const int index, BoutReal *sol){
  //Initialise the return data
  for(int i=0;i<nlocal;i++){
    sol[i]=start[i];
  }

  //Construct the solution
  for(int curStage=0;curStage<getStageCount();curStage++){
    if(resultCoeffs[curStage][index] == 0.) continue; //Real comparison not great
    BoutReal fac=dt*resultCoeffs[curStage][index];
    for(int i=0;i<nlocal;i++){
      sol[i]=sol[i]+fac*steps[curStage][i];
    }
  }
  
}

void RKScheme::constructOutputs(const BoutReal *start, const BoutReal dt, 
				const int indexFollow, const int indexAlt, 
				BoutReal *solFollow, BoutReal *solAlt){
  //Initialise the return data
  for(int i=0;i<nlocal;i++){
    solFollow[i]=start[i];
    solAlt[i]=start[i];
  }

  //Construct the solution
  for(int curStage=0;curStage<getStageCount();curStage++){
    BoutReal facFol=dt*resultCoeffs[curStage][indexFollow];
    BoutReal facAlt=dt*resultCoeffs[curStage][indexAlt];

    for(int i=0;i<nlocal;i++){
      solFollow[i]=solFollow[i]+facFol*steps[curStage][i];
      solAlt[i]=solAlt[i]+facAlt*steps[curStage][i];
    }
  }
  
}

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
    }
    output<<setw(10)<<timeCoeffs[i]<<" | "<<setw(10)<<tmp<<endl;
    if(fabs(timeCoeffs[i]-tmp)>atol) warn=true;
  }

  //Optional warning
  if(warn){
    output<<string(50,'=')<<endl;
    output<<"WARNING: Stage/Time coefficients not consistent"<<endl;
    output<<string(50,'=')<<endl;
    warn=false;
  }

  //Check results coefficients
  output<<string(50,'-')<<endl;
  output<<"Results coefficients (should be 1)"<<endl;
  for(int j=0;j<getNumOrders();j++){
    BoutReal tmp=0;
    for(int i=0;i<getStageCount();i++){
      tmp+=resultCoeffs[i][j];
    }
    output<<"Order : "<<j<<" = "<<tmp<<endl;
    if(fabs(1.0-tmp)>atol) warn=true;
  }

  //Optional warning
  if(warn){
    output<<string(50,'=')<<endl;
    output<<"WARNING: Result coefficients not consistent"<<endl;
    output<<string(50,'=')<<endl;
  }

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
  }

  //Divider
  output<<string(width,'-')<<" | "<<string(width*getStageCount(),'-')<<endl;

  //Result coeffs for each order
  for(int i=0;i<getNumOrders();i++){
    output<<setw(width)<<i<<" | ";
    for(int j=0;j<getStageCount();j++){
      output<<setw(width)<<resultCoeffs[j][i];
    }
    output<<endl;
  }

  //Footer
  output<<string(totalWidth,'-')<<endl;
  output<<endl;
}

void RKScheme::zeroSteps(){
  for(int i=0;i<getStageCount();i++){
    for(int j=0;j<nlocal;j++){
      steps[i][j]=0.;
    }
  }
}

