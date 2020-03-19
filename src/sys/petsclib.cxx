
#ifdef BOUT_HAS_PETSC

#include "boutcomm.hxx"
#include "bout/openmpwrap.hxx"
#include <bout/petsclib.hxx>

#include <output.hxx>

// Define all the static member variables
int PetscLib::count = 0;
char PetscLib::help[] = "BOUT++: Uses finite difference methods to solve plasma fluid problems in curvilinear coordinates";
int *PetscLib::pargc = nullptr;
char ***PetscLib::pargv = nullptr;
PetscLogEvent PetscLib::USER_EVENT = 0;

PetscLib::PetscLib() {
  BOUT_OMP(critical(PetscLib))
  {
    if(count == 0) {
      // Initialise PETSc

      output << "Initialising PETSc\n";
      PETSC_COMM_WORLD = BoutComm::getInstance()->getComm();
      PetscInitialize(pargc,pargv,PETSC_NULL,help);
      PetscLogEventRegister("Total BOUT++",0,&USER_EVENT);
      PetscLogEventBegin(USER_EVENT,0,0,0,0);
    }
    count++;
  }
}

PetscLib::~PetscLib() {
  BOUT_OMP(critical(PetscLib))
  {
    count--;
    if(count == 0) {
      // Finalise PETSc
      output << "Finalising PETSc\n";
      PetscLogEventEnd(USER_EVENT,0,0,0,0);
      PetscFinalize();
    }
  }
}

void PetscLib::cleanup() {
  BOUT_OMP(critical(PetscLib))
  {
    if(count > 0) {
      output << "Finalising PETSc. Warning: Instances of PetscLib still exist.\n";
      PetscLogEventEnd(USER_EVENT,0,0,0,0);
      PetscFinalize();

      count = 0; // ensure that finalise is not called again later
    }
  }
}

#endif // BOUT_HAS_PETSC

