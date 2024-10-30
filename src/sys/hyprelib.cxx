#include "bout/build_defines.hxx"

#if BOUT_HAS_HYPRE

#include <bout/hyprelib.hxx>

#include "bout/boutcomm.hxx"
#include "bout/openmpwrap.hxx"
#include "bout/options.hxx"
#include "bout/output.hxx"
#include "bout/unused.hxx"

#include <HYPRE.h>
#include <HYPRE_utilities.h>
#include <_hypre_utilities.h>

namespace bout {
// Define all the static member variables
int HypreLib::count = 0;

#if BOUT_HAS_CUDA
static constexpr auto BOUT_HYPRE_EXEC = HYPRE_EXEC_DEVICE;
static constexpr auto BOUT_HYPRE_MEMORY = HYPRE_MEMORY_DEVICE;
#else
static constexpr auto BOUT_HYPRE_EXEC = HYPRE_EXEC_HOST;
static constexpr auto BOUT_HYPRE_MEMORY = HYPRE_MEMORY_HOST;
#endif

HypreLib::HypreLib() {
  BOUT_OMP_SAFE(critical(HypreLib))
  {
    if (count == 0) { // Initialise once
      output_progress.write("Initialising Hypre\n");
      HYPRE_Init();
      hypre_HandleDefaultExecPolicy(hypre_handle()) = BOUT_HYPRE_EXEC;
      hypre_HandleMemoryLocation(hypre_handle()) = BOUT_HYPRE_MEMORY;
    }
    count++;
  }
}

HypreLib::HypreLib([[maybe_unused]] const HypreLib& other) noexcept {
  BOUT_OMP_SAFE(critical(HypreLib))
  {
    // No need to initialise Hypre, because it must already be initialised
    count++; // Copying, so increase count
  }
}

HypreLib::HypreLib([[maybe_unused]] HypreLib&& other) noexcept {
  BOUT_OMP_SAFE(critical(HypreLib))
  {
    // No need to initialise Hypre, because it must already be initialised
    count++; // Creating a new Hyprelib object; other will be deleted
  }
}

HypreLib::~HypreLib() {
  BOUT_OMP_SAFE(critical(HypreLib))
  {
    count--;
    if (count == 0) {
      output_progress.write("Finalising Hypre\n");
      HYPRE_Finalize();
    }
  }
}

void HypreLib::cleanup() {
  BOUT_OMP_SAFE(critical(HypreLib))
  {
    if (count > 0) {
      output << "Finalising Hypre. Warning: Instances of HypreLib still exist.\n";
      HYPRE_Finalize();
      count = 0; // ensure that finalize is not called again later
    }
  }
}
} // namespace bout

#endif // BOUT_HAS_HYPRE
