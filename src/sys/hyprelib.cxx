#include "bout/build_config.hxx"

#if BOUT_HAS_HYPRE

#include "boutcomm.hxx"
#include "options.hxx"
#include "bout/openmpwrap.hxx"
#include <bout/hyprelib.hxx>

#include <output.hxx>

// Define all the static member variables
int HypreLib::count = 0;

HypreLib::HypreLib() {
  BOUT_OMP(critical(HypreLib))
  {
    if(count == 0) {
      // Initialise Hypre

      output << "Initialising Hypre\n";
      HYPRE_Init();

    }

    count++;
  }
}

HypreLib::~HypreLib() {
  BOUT_OMP(critical(HypreLib))
  {
    count--;
    if(count == 0) {
      output << "Finalising Hypre\n";
      HYPRE_Finalize();
    }
  }
}

void HypreLib::cleanup() {
  BOUT_OMP(critical(HypreLib))
  {
    if(count > 0) {
      output << "Finalising Hypre. Warning: Instances of HypreLib still exist.\n";
      HYPRE_Finalize();

      count = 0; // ensure that finalise is not called again later
    }
  }
}

#endif // BOUT_HAS_HYPRE
