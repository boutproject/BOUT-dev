
#ifdef BOUT_HAS_PETSC

#include "boutcomm.hxx"
#include <bout/petsclib.hxx>

#include <output.hxx>

// Define all the static member variables
int PetscLib::count = 0;
char PetscLib::help[] = "BOUT++: Uses finite difference methods to solve plasma fluid problems in curvilinear coordinates";
int *PetscLib::pargc = nullptr;
char ***PetscLib::pargv = nullptr;
PetscLogEvent PetscLib::USER_EVENT = 0;

PetscLib::PetscLib(Options* opt) {
  if(count == 0) {
    // Initialise PETSc
    
    output << "Initialising PETSc\n";
    PETSC_COMM_WORLD = BoutComm::getInstance()->getComm();
    PetscInitialize(pargc,pargv,PETSC_NULL,help);
    PetscLogEventRegister("Total BOUT++",0,&USER_EVENT);
    PetscLogEventBegin(USER_EVENT,0,0,0,0);
  }

  if (count == 0 or opt != nullptr) {
    // Pass options to Petsc's global options database.
    // (PetscOptions type exists for non-global options, but its use is not discussed in
    // the Petsc manual, so ignoring the possibility here.)

    if (opt == nullptr) {
      // Options read by default from the [petsc] section of the input
      opt = Options::getRoot()->getSection("petsc");
    }

    Options& options = *opt;

    // Pass all options in the section to Petsc
    for (auto& i : options.getChildren()) {
      if (not i.second.isValue()) {
        throw BoutException("Found subsection %s in %s when reading Petsc options - only "
            "values are allowed in the Petsc options, not subsections",
            i.first.c_str(), options.str().c_str());
      }
      // Note, option names in the input file don't start with "-", but need to be passed
      // to Petsc with "-" prepended
      PetscErrorCode ierr;
      if (lowercase(i.second) == "true") {
        // Petsc flag with no value
        ierr = PetscOptionsSetValue(nullptr, ("-"+i.first).c_str(), nullptr);
      } else {
        // Option with actual value to pass
        ierr = PetscOptionsSetValue(nullptr, ("-"+i.first).c_str(),
            i.second.as<std::string>().c_str());
      }
      if (ierr) {
        throw BoutException("PetscOptionsSetValue returned error code %i", ierr);
      }
    }
  }
  count++;
}

PetscLib::~PetscLib() {
  count--;
  if(count == 0) {
    // Finalise PETSc
    output << "Finalising PETSc\n";
    PetscLogEventEnd(USER_EVENT,0,0,0,0);
    PetscFinalize();
  }
}

void PetscLib::cleanup() {
  if(count == 0)
    return; // Either never initialised, or already cleaned up

  output << "Finalising PETSc. Warning: Instances of PetscLib still exist.\n";
  PetscLogEventEnd(USER_EVENT,0,0,0,0);
  PetscFinalize();
  
  count = 0; // ensure that finalise is not called again later
}

#endif // BOUT_HAS_PETSC

