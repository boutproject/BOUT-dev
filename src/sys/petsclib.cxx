#include "bout/build_defines.hxx"

#if BOUT_HAS_PETSC

#include "bout/boutcomm.hxx"
#include "bout/boutexception.hxx"
#include "bout/openmpwrap.hxx"
#include "bout/options.hxx"
#include "bout/output.hxx"
#include "bout/petsclib.hxx"

#include <petscerror.h>
#include <petscksp.h>
#include <petsclog.h>
#include <petscoptions.h>
#include <petscsnes.h>
#include <petscsys.h>
#include <petscversion.h>

#include <string>

namespace {
constexpr const char* PetscLibHelp =
    "BOUT++: Uses finite difference methods to solve plasma fluid "
    "problems in curvilinear coordinates";

void setPetscOptions(Options& options, const std::string& prefix) {
  // Pass all options in the section to PETSc
  for (const auto& child : options.getChildren()) {
    if (not child.second.isValue()) {
      throw BoutException("Found subsection {} in {} when reading PETSc options - only "
                          "values are allowed in the PETSc options, not subsections",
                          child.first, options.str());
    }

    // Note, option names in the input file don't start with "-", but need to be passed
    // to PETSc with "-" prepended
    auto petsc_option_name = "-" + prefix + child.first;

    auto str_value = child.second.as<std::string>();
    // "true" is the value given to an option with no value, when read from BOUT.inp. Also
    // when nullptr is passed to PetscOptionsSetValue for a boolean option, it defaults to
    // true so we should always be OK passing nullptr for null or "true".
    const char* value = str_value == "true" ? nullptr : str_value.c_str();

#if PETSC_VERSION_GE(3, 7, 0)
    const auto ierr = PetscOptionsSetValue(nullptr, petsc_option_name.c_str(), value);
#else
    // no PetscOptions as first argument
    const auto ierr = PetscOptionsSetValue(petsc_option_name.c_str(), value);
#endif
    if (ierr != 0) {
      throw BoutException("PetscOptionsSetValue returned error code {} when setting {}",
                          ierr, petsc_option_name);
    }
  }
}
} // namespace

PetscLib::PetscLib(Options* opt) {
  BOUT_OMP_SAFE(critical(PetscLib))
  {
    if (count == 0) {
      // Initialise PETSc

      // Load global PETSc options from the [petsc] section of the input
      // Note: This should be before PetscInitialize so that some options
      //       can modify initialization e.g. -log_view.
      setPetscOptions(Options::root()["petsc"], "");

      setenv("PETSC_OPTIONS", "-options_left 0", 1);

      output << "Initialising PETSc\n";
      PETSC_COMM_WORLD = BoutComm::getInstance()->getComm();
      PetscInitialize(pargc, pargv, nullptr, PetscLibHelp);
      PetscPopSignalHandler();

      PetscLogEventRegister("Total BOUT++", 0, &USER_EVENT);
      PetscLogEventBegin(USER_EVENT, 0, 0, 0, 0);
    }

    if (opt != nullptr and opt->isSection()) {
      // Use options specific to this PetscLib
      // Pass options to PETSc's global options database, with a unique prefix, that will be
      // passed to a KSP later.
      // (PetscOptions type exists for non-global options, but apparently is only for user
      // options, and cannot be passed to KSP, etc. Non-global options can be passed by
      // defining a custom prefix for the options string, and then passing that to the KSP.)

      options_prefix = "boutpetsclib_" + opt->str();

      setPetscOptions((*opt)["petsc"], options_prefix);
    }

    count++;
  }
}

PetscLib::~PetscLib() {
  BOUT_OMP_SAFE(critical(PetscLib))
  {
    count--;
    if (count == 0) {
      // Finalise PETSc
      output << "Finalising PETSc\n";
      PetscLogEventEnd(USER_EVENT, 0, 0, 0, 0);
      PetscFinalize();
    }
  }
}

void PetscLib::setOptionsFromInputFile(KSP& ksp) {
  assertIerr(KSPSetOptionsPrefix(ksp, options_prefix.c_str()), "KSPSetOptionsPrefix");

  assertIerr(KSPSetFromOptions(ksp), "KSPSetFromOptions");
}

void PetscLib::setOptionsFromInputFile(SNES& snes) {
  BOUT_DO_PETSC(SNESSetOptionsPrefix(snes, options_prefix.c_str()));

  BOUT_DO_PETSC(SNESSetFromOptions(snes));
}

void PetscLib::cleanup() {
  BOUT_OMP_SAFE(critical(PetscLib))
  {
    if (count > 0) {
      output << "Finalising PETSc. Warning: Instances of PetscLib still exist.\n";
      PetscLogEventEnd(USER_EVENT, 0, 0, 0, 0);
      PetscFinalize();

      count = 0; // ensure that finalise is not called again later
    }
  }
}

BoutException PetscLib::SNESFailure(SNES& snes) {
  SNESConvergedReason reason = SNES_CONVERGED_ITERATING;
  BOUT_DO_PETSC(SNESGetConvergedReason(snes, &reason));
#if PETSC_VERSION_GE(3, 15, 0)
  const char* message{nullptr};
  BOUT_DO_PETSC(SNESGetConvergedReasonString(snes, &message));
#else
  const char* message{""};
#endif
  return BoutException("SNES failed to converge. Reason: {} ({:d})", message,
                       static_cast<int>(reason));
}
#endif // BOUT_HAS_PETSC
