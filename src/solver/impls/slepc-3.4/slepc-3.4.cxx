/**************************************************************************
 * Interface to SLEPc solver
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/


#ifdef BOUT_HAS_SLEPC_3_4
//Hacked together by <DD>

#include "slepc-3.4.hxx"

#include <globals.hxx>

#include <stdlib.h>

#include <interpolation.hxx> // Cell interpolation
#include <msg_stack.hxx>
#include <output.hxx>
#include <boutcomm.hxx>

static char help[] = "BOUT++: Uses finite difference methods to solve plasma fluid problems in curvilinear coordinates";

//The callback function for the shell matrix-multiply operation
//A simple wrapper around the SlepcSolver advanceStep routine
PetscErrorCode advanceStepWrapper(Mat matOperator, Vec inData, Vec outData){
  PetscFunctionBegin;
  SlepcSolver* ctx;
  MatShellGetContext(matOperator,(void**) &ctx); //Here we set the ctx pointer to the solver instance
  PetscFunctionReturn(ctx->advanceStep(matOperator,inData,outData)); //Actually advance
}

//The callback function for the eigenvalue comparison
//A simple wrapper around the SlepcSolver compareEigs routine
PetscErrorCode compareEigsWrapper(PetscScalar ar, PetscScalar ai, PetscScalar br, PetscScalar bi,
				  PetscInt *res, void *ctx){
  PetscFunctionBegin;
  //Cast context as SlepcSolver and call the actual compare routine
  //  tmpCompare=(static_cast<SlepcSolver *>(ctx))->compareEigs(ar,ai,br,bi);
  SlepcSolver* myCtx;
  myCtx=(SlepcSolver*)ctx;
  myCtx->compareState=myCtx->compareEigs(ar,ai,br,bi);
 
  res = &myCtx->compareState;
  PetscFunctionReturn(0);
}


//The callback function for the monitor
//A simple wrapper around the SlepcSolver compareEigs routine
PetscErrorCode monitorWrapper(EPS eps, PetscInt its, PetscInt nconv,
				  PetscScalar *eigr, PetscScalar *eigi,
				  PetscReal* errest, PetscInt nest, void *mctx){
  PetscFunctionBegin;
  //Cast context as SlepcSolver and call the actual compare routine
  SlepcSolver* myCtx;
  myCtx=(SlepcSolver*)mctx;
  myCtx->monitor(its,nconv,eigr,eigi,errest,nest);
  PetscFunctionReturn(0);
}

SlepcSolver::SlepcSolver(){
  has_constraints = false;
  initialised = false;
  
  //Setup actual advance solver -- We have to do this now to ensure that any 
  //calls that are passed through to the actual solver are allowed
  advanceSolver=NULL;
  f0=NULL;
  f1=NULL;
  Options *globalOptions = Options::getRoot();
  Options *options=globalOptions->getSection("slepc");
  setAdvanceSolver(options);
}

SlepcSolver::~SlepcSolver(){
  if(initialised){
    //Free memory
    if(eps){EPSDestroy(&eps);};
    if(shellMat){MatDestroy(&shellMat);};
    if(advanceSolver){delete advanceSolver;};
    if(f0){delete[] f0;};
    if(f1){delete[] f1;};
    initialised = false;
  }
}

int SlepcSolver::init(bool restarting, int NOUT, BoutReal TIMESTEP) {

#ifdef CHECK
  int msg_point = msg_stack.push("Initialising SLEPc-3.4 solver");
#endif

  //Report initialisation
  output.write("Initialising SLEPc-3.4 solver\n");  
  if(selfSolve){
    Solver::init(restarting,NOUT,TIMESTEP);
  }

  //If no advanceSolver then can only advance one step at a time
  if(selfSolve){
    NOUT=1;
  }

  //Save for use later
  nout = NOUT;
  tstep = TIMESTEP;

  //Read options
  readOptions();
  comm=PETSC_COMM_WORLD;

  //Initialise advanceSolver if not self
  if(!selfSolve){
    //--> First set solver parameters which are usually done automatically
    copySettingsToAdvanceSolver();
    //--> Now call actual init
    advanceSolver->init(restarting,NOUT,TIMESTEP);
  }

  //Calculate grid sizes
  localSize=getLocalN();

  //Also create vector for derivs etc. if SLEPc in charge of solving
  if(selfSolve){    
    // Allocate memory
    f0 = new BoutReal[localSize];
    f1 = new BoutReal[localSize];
  }
  
  // Get total problem size
  int neq;
  if(MPI_Allreduce(&localSize, &neq, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    throw BoutException("MPI_Allreduce failed in SlepcSolver::init");
  }
  
  output.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n",
	       n3Dvars(), n2Dvars(), neq, localSize);
  
  //Create EPS solver
  createEPS();

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  //Return ok
  return 0;
}

//Read in the slepc options and set various things
void SlepcSolver::readOptions(){
  Options *globalOptions = Options::getRoot();
  Options *options=globalOptions->getSection("slepc");

  //Slepc settings
  options->get("nEig",nEig,0);//0 means keep the current value, i.e. autoset
  options->get("tol",tol,1.0e-6);//Tolerance --> Note this is on SLEPc eig not BOUT
  options->get("maxIt",maxIt,PETSC_DECIDE);
  options->get("targRe",targRe,0.0); //Target frequency when using user eig comparison
  options->get("targIm",targIm,0.0); //Target growth rate when using user eig comparison
  
  //Convert bout targs to slepc
  if(targRe==0.0 && targIm==0.0){
    target=999.0;
  }else{
    PetscScalar slepcRe,slepcIm;
    boutToSlepc(targRe,targIm,slepcRe,slepcIm);
    dcomplex tmp(slepcRe,slepcIm);
    target=abs(tmp);
  }
  options->get("target",target,target); //If 999 we don't set the target. This is SLEPc eig target
  
  options->get("userWhich",userWhich,false);

  //Generic settings
  options->get("useInitial",useInitial,true);
  options->get("debugMonitor",debugMonitor,false);
}

void SlepcSolver::copySettingsToAdvanceSolver(){
  //Exit if don't have an advance solver
  if(selfSolve){return;}

  advanceSolver->setRHS(phys_run);
  advanceSolver->setPrecon(prefunc);

  //Could attach a monitor here, if so should we remove it from
  //the SlepcSolver?
  if(debugMonitor){advanceSolver->addMonitor(bout_monitor);}

  return;
}

int SlepcSolver::run() {
  //Now the basic idea with slepc is that:
  //Whilst n_eig_converged<n_eig_desired do
  //1. Let slepc set the initial fields
  //2. Use the advanceSolver to evolve fields over a certain time
  //3. Package up the evolved fields for slepc to analyse
  //Once converged data found:
  //1. Get converged eigenvalues/eigenvectors
  //2. Write to file
  //3. Clean up
  //--> The first section is handled by calling EPSSolve(eps) with appropriate shellMat
  //--> The second section has to be handled by this solver

  //Find the eigenvalues
  EPSSolve(eps);

  //The following prints the used solver settings
  EPSView(eps,PETSC_VIEWER_STDOUT_WORLD);

  //Analyse and dump to file
  analyseResults();

  //Return ok
  return 0;
}

//Create a solver instance based on options
void SlepcSolver::setAdvanceSolver(Options *options){
  if(options==NULL){
    options=Options::getRoot()->getSection("slepc");
  }

  string typeTmp;
  //  SolverType defType=SolverFactory::getInstance()->getDefaultSolverType();
  //The above would be nice but doesn't currently seem to work so do:
  SolverType defType=SOLVERRK4;
  SolverType selfType=SOLVERSLEPCSELF;

  options->get("advancesolver",typeTmp,string(defType));
  SolverType type=typeTmp.c_str();

  if(typeTmp == selfType){
    selfSolve=true;
    //Could we set advanceSolver=this in order to avoid branches on selfSolve later on?
  }else{
    selfSolve=false;
    //throw BoutException("Don't currently support use of solvers other than 'self' with SlepcSolver.");
    //Guard to protect against using anything that can't reset. Note to query the solver we have
    //to create an instance which means we must destroy it if we want to make a new instance.
    advanceSolver=SolverFactory::getInstance()->createSolver(type,NULL);
    if(!advanceSolver->canReset){
      output<<"WARNING:: CURRENTLY DON'T SUPPORT advanceSolver="<<type<<" --> Overriding to use "<<defType<<endl;
      delete advanceSolver;
      type=defType;
      advanceSolver=SolverFactory::getInstance()->createSolver(type,NULL);
    }
    //advanceSolver=SolverFactory::getInstance()->createSolver(type,NULL);
  }
}

//This routine takes a Vec type object of length localSize and 
//unpacks it into the local fields
void SlepcSolver::vecToFields(Vec &inVec){
  /*
    Whilst this routine does indeed populate the field variables
    most (/all?) solvers overwrite this data on call to run() as
    they actually evolve a single vector which they unpack (not
    unlike this routine). The typical run work flow is something
    like:
      1. Unpack internal vector A into fields.
      2. Use phys_run to calculate ddt of fields.
      3. Pack ddt into an internal vector B.
      4. Use B to evolve A.
    Unfortunately each solver implementation uses a different 
    internal vector approach which means it's a implementation
    specific thing which makes it hard/impossible to actually
    set from here. Possible solutions are:
      1. Add a virtual function to solver, which each implementation
         overrides, which just copies field data into appropriate
	 internal vector or copies out of internal vector into 
	 fields (for fieldsToVec).
      2. ?
   */

  //Get pointer to data
  PetscScalar *point;
  VecGetArray(inVec,&point);

  //Copy data from point into fields
  if(selfSolve){
    load_vars(point);
  }else{
    advanceSolver->load_vars(point);
    advanceSolver->resetInternalFields();
  }

  //Restore array
  VecRestoreArray(inVec,&point);
}

//This routine packs the local fields into a vector
void SlepcSolver::fieldsToVec(Vec &outVec){
  //Get pointer to data
  PetscScalar *point;
  VecGetArray(outVec,&point);

  //Copy fields into point
  if(selfSolve){
    save_vars(point);
  }else{
    advanceSolver->save_vars(point);
  }

  //Restore array
  VecRestoreArray(outVec,&point);
}

//Create a shell matrix operator
void SlepcSolver::createShellMat(){
  output<<"Creating shellMat with local size : "<<localSize<<endl;//<<" ("<<localSize2D<<"+"<<localSize3D<<")"<<endl;

  //Create the shell matrix
  MatCreateShell(comm,localSize,localSize,PETSC_DETERMINE,
  		 PETSC_DETERMINE,this,&shellMat); //Note we pass the this reference as the matrix context

  //Define the mat_mult operation --> Define what routine return M.x, where M
  //is the time advance operator and x are the initial conditions
  MatShellSetOperation(shellMat,MATOP_MULT,(void(*)())&advanceStepWrapper);

  //The above function callback can cause issues as member functions have a hidden "this" 
  //argument which means when Slepc calls advanceStep(Mat,Vec,Vec) this is redefined to Mat,
  //and the two Vecs are mangled, this messes up memory and the code crashes.
  //Alternatives include:
  //  1. Making advanceStep a static function
  //  2. Making advanceStep a non-member function
  //These alternatives generally divorce the advanceStep from the SlepcSolver class
  //which might make it difficult to access required data.
  //We've therefore gone with a third option where we define the callback to be a
  //non-member function which then gets the context pointer attached to shellMat
  //which points to the SlepcSolver instance. This can then be used to call the
  //advanceStep member function as if it were called within this.
}

//Create an EPS Solver
void SlepcSolver::createEPS(){
  //First need to create shell matrix
  createShellMat();

  //Now construct EPS
  EPSCreate(comm,&eps);
  EPSSetOperators(eps,shellMat,NULL);
  EPSSetProblemType(eps,EPS_NHEP);//Non-hermitian
  
  //Probably want to read options and set EPS properties
  //at this point.
  EPSSetDimensions(eps,nEig,PETSC_DECIDE,PETSC_DECIDE);
  EPSSetTolerances(eps,tol,maxIt);
  if(! (target==999)){
    EPSSetTarget(eps,target);
  }
  
  //Set the user comparison function
  if(userWhich){
    EPSSetEigenvalueComparison(eps,&compareEigsWrapper,this);
    EPSSetWhichEigenpairs(eps,EPS_WHICH_USER);
  }

  //Update options from command line
  EPSSetFromOptions(eps);

  //Register a monitor
  EPSMonitorSet(eps,&monitorWrapper,this,NULL);

  //Should probably call a routine here which interrogates eps
  //to determine the important settings that have been used and dump
  //the settings to screen/file/dmp?
  //I think there may be a Slepc flag which will do this (to screen)
  //but not sure if we can force this is the code (without messing with argv).

  //Set initial space i.e. first guess
  if(useInitial){
    Vec initVec, rightVec;
    MatGetVecs(shellMat,&rightVec,&initVec);
    fieldsToVec(initVec);
    EPSSetInitialSpace(eps,1,&initVec);
    VecDestroy(&initVec); 
    VecDestroy(&rightVec);
  };
}

//This routine takes initial conditions provided by SLEPc, uses this to set the fields,
//advances them with the attached solver and then returns the evolved fields in a slepc
//structure.
int SlepcSolver::advanceStep(Mat &matOperator, Vec &inData, Vec &outData){
  //First unpack input into fields
  vecToFields(inData);

  //Now advance
  int retVal;
  if(selfSolve){
    retVal=run_rhs(0.0);
    //Here we add dt*ddt(Fields) to fields to advance solution (cf. Euler)
    save_vars(f0);
    save_derivs(f1);
    for(int iVec=0;iVec<localSize;iVec++){ //THIS WILL IGNORE THE X-BOUNDARIES -- IS THAT OK?
      f0[iVec]+=f1[iVec]*tstep;
    }
    load_vars(f0);
  }else{
    // advanceSolver->copyFieldsToInternal(); //This resets the solver and sets initial fields -- Not implemented anywhere yet!
    retVal=advanceSolver->run();
  }

  //Now pack evolved fields into output
  fieldsToVec(outData);

  //Return
  return retVal;
}

int SlepcSolver::compareEigs(PetscScalar ar, PetscScalar ai, PetscScalar br, PetscScalar bi){
  BoutReal arBout, aiBout, brBout, biBout;
  
  //First convert to BOUT values
  slepcToBout(ar,ai,arBout,aiBout);
  slepcToBout(br,bi,brBout,biBout);

  //Now we calculate the distance between eigenvalues and target.
  BoutReal da, db;
  da=sqrt(pow(arBout-targRe,2)+pow(aiBout-targIm,2));
  db=sqrt(pow(brBout-targRe,2)+pow(biBout-targIm,2));

  // output<<"Picking between "<<ar<<"+i"<<ai<<" and "<<br<<"+i"<<bi<<endl;
  // output<<"(targRe="<<targRe<<", targIm="<<targIm<<")"<<endl;
  // output<<"   da : "<<da<<endl;
  // output<<"   db : "<<db<<endl;
  //Now we decide which eigenvalue is preferred.
  int retVal;

  //Largest growth rate
  // if(aiBout>biBout){
  //   retVal=-1;
  // }else if(biBout>aiBout){
  //   retVal=1;
  // }else{
  //   retVal=0;
  // }
  // return retVal;


  //Smallest distance from complex target
  //If prefer B we return +ve
  if(da>db){
    retVal=1;
  //If prefer A we return -ve
  }else if(db>da){
    retVal=-1;
  //If we don't prefer either we return 0
  }else{
    retVal=0;
  };
  //output<<"--> Return "<<retVal<<endl;
  return retVal;
}

void SlepcSolver::monitor(PetscInt its, PetscInt nconv, PetscScalar eigr[], PetscScalar eigi[], PetscReal errest[], PetscInt nest){
  static int nConvPrev=0;

  //No output until after first iteration
  if(its<1){return;}

  BoutReal reEigBout, imEigBout;
  string joinNum, joinNumSlepc;
  slepcToBout(eigr[nconv],eigi[nconv],reEigBout,imEigBout);
  if(imEigBout<0){
    joinNum="";
  }else{
    joinNum="+";
  }

  //This line more or less replicates the normal slepc output (when using -eps_monitor)
  //but reports Bout eigenvalues rather than the Slepc values. Note we haven't changed error estimate.
  output<<" "<<its<<" nconv="<<nconv<<"\t first unconverged value (error) "<<reEigBout<<joinNum<<imEigBout<<"i\t ("<<errest[nconv]<<")"<<endl;

  //The following can be quite noisy so may want to add a flag to disable/enable.
  int newConv=nconv-nConvPrev;
  if(newConv>0){
    output<<"Found "<<newConv<<" new converged eigenvalues:"<<endl;

    for (PetscInt i=nConvPrev;i<nconv;i++){
      slepcToBout(eigr[i],eigi[i],reEigBout,imEigBout);
      if(imEigBout<0){
	joinNum="";
      }else{
	joinNum="+";
      }
      if(eigi[i]<0){
	joinNumSlepc="";
      }else{
	joinNumSlepc="+";
      }
      output<<"\t"<<i<<"\t: "<<eigr[i]<<joinNumSlepc<<eigi[i]<<"i --> "<<reEigBout<<joinNum<<imEigBout<<"i"<<endl;
    }
  }

  //Update the number of converged modes already investigated.
  nConvPrev=nconv;
};

//Convert a slepc eigenvalue to a BOUT one
void SlepcSolver::slepcToBout(PetscScalar &reEigIn, PetscScalar &imEigIn,
			      BoutReal &reEigOut, BoutReal &imEigOut){

  //The slepc eigenvalue is actually
  //Exp(-i*Eig_Bout*tstep)
  //where Eig_BOUT is the actual eigenvalue and tstep is the time step
  //the solution is evolved over.
  dcomplex slepcEig(reEigIn,imEigIn), ci(0.0,1.0);
  dcomplex boutEig;
  boutEig=ci*log(slepcEig)/(tstep*nout);
  
  //Set return values
  reEigOut=boutEig.Real();
  imEigOut=boutEig.Imag();
}

//Convert a BOUT++ eigenvalue to a Slepc one
void SlepcSolver::boutToSlepc(BoutReal &reEigIn, BoutReal &imEigIn,
			      PetscScalar &reEigOut, PetscScalar &imEigOut){

  dcomplex boutEig(reEigIn,imEigIn), ci(0.0,1.0);
  dcomplex slepcEig;

  slepcEig=exp(-ci*boutEig*tstep*nout);

  //Set return values
  reEigOut=slepcEig.Real();
  imEigOut=slepcEig.Imag();
}


//Interrogate eps to find out how many eigenvalues we've found etc.
void SlepcSolver::analyseResults(){
  PetscInt nEigFound;

  //Find how many eigenvalues have been found
  EPSGetConverged(eps,&nEigFound);

  //Now loop over each converged eigenpair and output eigenvalue
  if(nEigFound>0){
    output<<"Converged eigenvalues :"<<endl;
    output<<"\tIndex\tSlepc eig (mag.)\t\t\tBOUT eig (mag.)"<<endl;

    iteration=0;

    //Declare and create vectors to store eigenfunctions
    Vec vecReal, vecImag;
    MatGetVecs(shellMat,&vecReal,&vecImag);

    //This allows us to set the simtime in bout++.cxx directly
    //rather than calling the monitors which are noisy |--> Not very nice way to do this
    extern BoutReal simtime;

    for(PetscInt iEig=0; iEig<nEigFound; iEig++){
      PetscScalar reEig, imEig;
      BoutReal reEigBout, imEigBout;
      EPSGetEigenvalue(eps,iEig,&reEig,&imEig);
      dcomplex slepcEig(reEig,imEig);
      output<<"\t"<<iEig<<"\t"<<reEig<<" "<<imEig<<"i\t("<<abs(slepcEig)<<")";
      slepcToBout(reEig,imEig,reEigBout,imEigBout);
      dcomplex boutEig(reEigBout,imEigBout);
      output<<"\t"<<reEigBout<<" "<<imEigBout<<"i\t("<<abs(boutEig)<<")"<<endl;

      //Get eigenvector
      EPSGetEigenvector(eps,iEig,vecReal,vecImag);

      //Write real part of eigen data
      //First dump real part to fields
      vecToFields(vecReal);
      //Set the simtime to omega
      simtime=reEigBout;
      //call_monitors(simtime,iteration,nEigFound*2);
      //Write to file
      dump.write();
      iteration++;

      //Now write imaginary part of eigen data
      //First dump imag part to fields
      vecToFields(vecImag);
      //Set the simtime to gamma
      simtime=imEigBout;
      //call_monitors(simtime,iteration,nEigFound*2);
      //Write to file
      dump.write();
      iteration++;
    }

    //Destroy vectors
    VecDestroy(&vecReal); VecDestroy(&vecImag);
  }
  else{
    output<<"Warning : No converged eigenvalues found!"<<endl;
  }
}

#endif // BOUT_HAS_SLEPC_3_4
