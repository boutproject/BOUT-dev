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
    output << "Initialising Hypre\n";
    HYPRE_Init();  
#ifdef BOUT_USE_CUDA
    hypre_HandleDefaultExecPolicy(hypre_handle()) = HYPRE_EXEC_DEVICE;
    HYPRE_MemoryLocation memory_location = HYPRE_MEMORY_DEVICE;
    hypre_HandleMemoryLocation(hypre_handle()) = memory_location;
#else
    hypre_HandleDefaultExecPolicy(hypre_handle()) = HYPRE_EXEC_HOST;
    HYPRE_MemoryLocation memory_location = HYPRE_MEMORY_HOST;
    hypre_HandleMemoryLocation(hypre_handle()) = memory_location;
#endif
  }
  count++;
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
      count = 0; // ensure that finalize is not called again later
    }
  }
}

#endif // BOUT_HAS_HYPRE
