#include "rkschemefactory.hxx"
#include "unused.hxx"
#include <bout/rkscheme.hxx>
#include <boutcomm.hxx>
#include <cmath>
#include <output.hxx>

////////////////////
// PUBLIC
////////////////////

//Initialise
RKScheme::RKScheme(Options *UNUSED(opts)) {
  // Currently not reading anything from the options here

  // Initialise internals
  dtfac = 1.0; // Time step factor
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
  steps.reallocate(getStageCount(), nlocal);
  zeroSteps();

  //Allocate array for storing alternative order result
  if (adaptive)
    resultAlt.reallocate(nlocal); // Result--alternative order

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
void RKScheme::setCurState(const Array<BoutReal> &start, Array<BoutReal> &out,
                           const int curStage, const BoutReal dt) {

  //Set the initial stage
  BOUT_OMP(parallel for)
  for(int i=0;i<nlocal;i++){
    out[i] = start[i];
  }

  //If on the first stage we don't need to modify the state
  if(curStage==0) return;//Don't realy need this as below loop won't execute
  
  //Construct the current state from previous results -- This is expensive
  for(int j=0;j<curStage;j++){
    if (std::abs(stageCoeffs(curStage, j)) < atol)
      continue;
    BoutReal fac = stageCoeffs(curStage, j) * dt;

    BOUT_OMP(parallel for)
    for(int i=0;i<nlocal;i++){
      out[i] = out[i] + fac * steps(j, i);
    }
  }
}

//Construct the system state at the next time
BoutReal RKScheme::setOutputStates(const Array<BoutReal> &start, const BoutReal dt,
                                   Array<BoutReal> &resultFollow) {
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
    constructOutput(start, dt, altInd, resultAlt);
  }

  //Get the error coefficient
  return getErr(resultFollow, resultAlt);
}

BoutReal RKScheme::updateTimestep(const BoutReal dt, const BoutReal err){
  return dtfac*dt*pow(rtol/(2.0*err),1.0/(order+1.0));
}

////////////////////
// PRIVATE
////////////////////

//Estimate the error, given two solutions
BoutReal RKScheme::getErr(Array<BoutReal> &solA, Array<BoutReal> &solB) {
  BoutReal err=0.;

  //If not adaptive don't care about the error
  if(!adaptive){return err;}

  //Get local part of relative error
  BoutReal local_err = 0.;

  // Note because the order of operation is not deterministic
  // we expect slightly different round-off error each time this
  // is called and hence the nrhs may no longer be exactly
  // repeatable with this parallelisation.
  BOUT_OMP(parallel for reduction(+:local_err))
  for(int i=0;i<nlocal;i++) {
    local_err +=
        std::abs(solA[i] - solB[i]) / (std::abs(solA[i]) + std::abs(solB[i]) + atol);
  }
  //Reduce over procs
  if(MPI_Allreduce(&local_err, &err, 1, MPI_DOUBLE, MPI_SUM, BoutComm::get())) {
    throw BoutException("MPI_Allreduce failed");
  }
  //Normalise by number of values
  err /= static_cast<BoutReal>(neq);

  return err;
}

void RKScheme::constructOutput(const Array<BoutReal> &start, const BoutReal dt,
                               const int index, Array<BoutReal> &sol) {
  //Initialise the return data
  BOUT_OMP(parallel for)
  for(int i=0;i<nlocal;i++){
    sol[i]=start[i];
  }

  //Construct the solution
  for(int curStage=0;curStage<getStageCount();curStage++){
    if (resultCoeffs(curStage, index) == 0.)
      continue; // Real comparison not great
    BoutReal fac = dt * resultCoeffs(curStage, index);
    BOUT_OMP(parallel for)
    for(int i=0;i<nlocal;i++){
      sol[i] = sol[i] + fac * steps(curStage, i);
    }
  }
  
}

void RKScheme::constructOutputs(const Array<BoutReal> &start, const BoutReal dt,
                                const int indexFollow, const int indexAlt,
                                Array<BoutReal> &solFollow, Array<BoutReal> &solAlt) {
  //Initialise the return data
  BOUT_OMP(parallel for)
  for(int i=0;i<nlocal;i++){
    solFollow[i]=start[i];
    solAlt[i]=start[i];
  }

  //Construct the solution
  for(int curStage=0;curStage<getStageCount();curStage++){
    BoutReal facFol = dt * resultCoeffs(curStage, indexFollow);
    BoutReal facAlt = dt * resultCoeffs(curStage, indexAlt);
    BOUT_OMP(parallel for)
    for(int i=0;i<nlocal;i++){
      solFollow[i] = solFollow[i] + facFol * steps(curStage, i);
      solAlt[i] = solAlt[i] + facAlt * steps(curStage, i);
    }
  }
  
}

//Check that the coefficients are consistent
void RKScheme::verifyCoeffs(){
  output<<endl;

  //Header
  output<<std::string(50,'-')<<endl;
  output<<"RK Coefficients consistency check"<<endl;
  output<<std::string(50,'-')<<endl;

  //Check time and stage coefficients
  output<<std::setw(10)<<"TimeCoeff"<<" | "<<std::setw(10)<<"SumStageCoeff"<<endl;
  bool warn=false;
  for(int i=0;i<getStageCount();i++){
    BoutReal tmp=0;
    for(int j=0;j<i;j++){
      tmp += stageCoeffs(i, j);
    }
    output<<std::setw(10)<<timeCoeffs[i]<<" | "<<std::setw(10)<<tmp<<endl;
    if(std::abs(timeCoeffs[i]-tmp)>atol) warn=true;
  }

  //Optional warning
  if(warn){
    output<<std::string(50,'=')<<endl;
    output<<"WARNING: Stage/Time coefficients not consistent"<<endl;
    output<<std::string(50,'=')<<endl;
    warn=false;
  }

  //Check results coefficients
  output<<std::string(50,'-')<<endl;
  output<<"Results coefficients (should be 1)"<<endl;
  for(int j=0;j<getNumOrders();j++){
    BoutReal tmp=0;
    for(int i=0;i<getStageCount();i++){
      tmp += resultCoeffs(i, j);
    }
    output<<"Order : "<<j<<" = "<<tmp<<endl;
    if(std::abs(1.0-tmp)>atol) warn=true;
  }

  //Optional warning
  if(warn){
    output<<std::string(50,'=')<<endl;
    output<<"WARNING: Result coefficients not consistent"<<endl;
    output<<std::string(50,'=')<<endl;
  }

  //Footer
  output<<std::string(50,'-')<<endl;
  output<<endl;
}

void RKScheme::printButcherTableau(){
  int width=15;
  int totalWidth=width*(1+getStageCount())+3;
  output<<endl;
  
  //Header
  output<<std::string(totalWidth,'-')<<endl;
  output<<"Butcher Tableau"<<endl;
  output<<std::string(totalWidth,'-')<<endl;

  //Time and stage coeffs
  for(int i=0;i<getStageCount();i++){
    output<<std::setw(width)<<timeCoeffs[i]<<" | ";
    for(int j=0;j<getStageCount();j++){
      output << std::setw(width) << stageCoeffs(i, j);
    }
    output<<endl;
  }

  //Divider
  output<<std::string(width,'-')<<" | "<<std::string(width*getStageCount(),'-')<<endl;

  //Result coeffs for each order
  for(int i=0;i<getNumOrders();i++){
    output<<std::setw(width)<<i<<" | ";
    for(int j=0;j<getStageCount();j++){
      output << std::setw(width) << resultCoeffs(j, i);
    }
    output<<endl;
  }

  //Footer
  output<<std::string(totalWidth,'-')<<endl;
  output<<endl;
}

void RKScheme::zeroSteps(){
  for(int i=0;i<getStageCount();i++){
    for(int j=0;j<nlocal;j++){
      steps(i, j) = 0.;
    }
  }
}

