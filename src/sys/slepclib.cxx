#include "bout/build_config.hxx"

#if BOUT_HAS_SLEPC

#include <bout/slepclib.hxx>
#include <output.hxx>

// Define all the static member variables
int SlepcLib::count = 0;
char SlepcLib::help[] = "BOUT++: Uses finite difference methods to solve plasma fluid problems in curvilinear coordinates";
int* SlepcLib::pargc = nullptr;
char*** SlepcLib::pargv = nullptr;
PetscLogEvent SlepcLib::USER_EVENT = 0;

SlepcLib::SlepcLib() {
  if(count == 0) {
    // Initialise SLEPc

    output << "Initialising SLEPc\n";
    SlepcInitialize(pargc,pargv,PETSC_NULL,help);
    PetscLogEventRegister("Total BOUT++",0,&USER_EVENT);
    PetscLogEventBegin(USER_EVENT,0,0,0,0);
  }
  count++;
}

SlepcLib::~SlepcLib() {
  count--;
  if(count == 0) {
    // Finalise Slepc
    output << "Finalising SLEPc\n";
    PetscLogEventEnd(USER_EVENT,0,0,0,0);
    SlepcFinalize();
  }
}

void SlepcLib::cleanup() {
  if(count == 0)
    return; // Either never initialised, or already cleaned up

  output << "Finalising SLEPCc. Warning: Instances of SlepcLib still exist.\n";
  PetscLogEventEnd(USER_EVENT,0,0,0,0);
  SlepcFinalize();

  count = 0; // ensure that finalise is not called again later
}

#endif // BOUT_HAS_SLEPC
