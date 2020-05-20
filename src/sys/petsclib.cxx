
#ifdef BOUT_HAS_PETSC

#include "boutcomm.hxx"
#include "options.hxx"
#include "bout/openmpwrap.hxx"
#include <bout/petsclib.hxx>

#include <output.hxx>

// Define all the static member variables
int PetscLib::count = 0;
int PetscLib::unique_id = 0;
char PetscLib::help[] = "BOUT++: Uses finite difference methods to solve plasma fluid problems in curvilinear coordinates";
int *PetscLib::pargc = nullptr;
char ***PetscLib::pargv = nullptr;
PetscLogEvent PetscLib::USER_EVENT = 0;

PetscLib::PetscLib(Options* opt) {
  BOUT_OMP(critical(PetscLib))
  {
    if(count == 0) {
      // Initialise PETSc

      output << "Initialising PETSc\n";
      PETSC_COMM_WORLD = BoutComm::getInstance()->getComm();
      PetscInitialize(pargc,pargv,PETSC_NULL,help);
      PetscLogEventRegister("Total BOUT++",0,&USER_EVENT);
      PetscLogEventBegin(USER_EVENT,0,0,0,0);

    // Load global PETSc options from the [petsc] section of the input
    setPetscOptions(Options::root()["petsc"], "");
  }

  if (opt != nullptr and opt->isSection()) {
    // Use options specific to this PetscLib
    // Pass options to PETSc's global options database, with a unique prefix, that will be
    // passed to a KSP later.
    // (PetscOptions type exists for non-global options, but apparently is only for user
    // options, and cannot be passed to KSP, etc. Non-global options can be passed by
    // defining a custom prefix for the options string, and then passing that to the KSP.)

    options_prefix = "boutpetsclib" + std::to_string(unique_id) + "_";

    setPetscOptions(*opt, options_prefix);
    }
    count++;
  unique_id++;
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

void PetscLib::setOptionsFromInputFile(KSP& ksp) {
  auto ierr = KSPSetOptionsPrefix(ksp, options_prefix.c_str());
  if (ierr) {
    throw BoutException("KSPSetOptionsPrefix failed with error %i", ierr);
  }

  ierr = KSPSetFromOptions(ksp);
  if (ierr) {
    throw BoutException("KSPSetFromOptions failed with error %i", ierr);
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

void PetscLib::setPetscOptions(Options& options, const std::string& prefix) {
  // Pass all options in the section to PETSc
  for (auto& i : options.getChildren()) {
    if (not i.second.isValue()) {
      throw BoutException("Found subsection %s in %s when reading PETSc options - only "
          "values are allowed in the PETSc options, not subsections",
          i.first.c_str(), options.str().c_str());
    }

    // Note, option names in the input file don't start with "-", but need to be passed
    // to PETSc with "-" prepended
    auto petsc_option_name = "-"+prefix+i.first;

    // "true" is the value given to an option with no value, when read from BOUT.inp. Also
    // when nullptr is passed to PetscOptionsSetValue for a boolean option, it defaults to
    // true so we should always be OK passing nullptr for null or "true".
    const char* value =
      lowercase(i.second) == "true" ? nullptr : i.second.as<std::string>().c_str();

    PetscErrorCode ierr;
#if PETSC_VERSION_GE(3, 7, 0)
    ierr = PetscOptionsSetValue(nullptr, petsc_option_name.c_str(), value);
#else
// no PetscOptions as first argument
    ierr = PetscOptionsSetValue(petsc_option_name.c_str(), value);
#endif
    if (ierr) {
      throw BoutException("PetscOptionsSetValue returned error code %i when setting %s",
                          ierr, petsc_option_name.c_str());
    }
  }
}
#endif // BOUT_HAS_PETSC

