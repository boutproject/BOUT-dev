#ifndef BOUT_BUILD_OPTIONS_HXX
#define BOUT_BUILD_OPTIONS_HXX

#include "bout/build_defines.hxx"

// Convert macro to a string constant
#define STRINGIFY1(x) #x
#define STRINGIFY(x) STRINGIFY1(x)

namespace bout {
namespace build {
constexpr auto check_level = BOUT_CHECK_LEVEL;
constexpr auto openmp_schedule = STRINGIFY(BOUT_OPENMP_SCHEDULE);

constexpr auto has_fftw = static_cast<bool>(BOUT_HAS_FFTW);
constexpr auto has_gettext = static_cast<bool>(BOUT_HAS_GETTEXT);
constexpr auto has_lapack = static_cast<bool>(BOUT_HAS_LAPACK);
constexpr auto has_legacy_netcdf = static_cast<bool>(BOUT_HAS_LEGACY_NETCDF);
constexpr auto has_netcdf = static_cast<bool>(BOUT_HAS_NETCDF);
constexpr auto has_petsc = static_cast<bool>(BOUT_HAS_PETSC);
constexpr auto has_pretty_function = static_cast<bool>(BOUT_HAS_PRETTY_FUNCTION);
constexpr auto has_pvode = static_cast<bool>(BOUT_HAS_PVODE);
constexpr auto has_scorep = static_cast<bool>(BOUT_HAS_SCOREP);
constexpr auto has_uuid_system_generator =
    static_cast<bool>(BOUT_HAS_UUID_SYSTEM_GENERATOR);
constexpr auto has_slepc = static_cast<bool>(BOUT_HAS_SLEPC);
constexpr auto has_sundials = static_cast<bool>(BOUT_HAS_SUNDIALS);
constexpr auto use_backtrace = static_cast<bool>(BOUT_USE_BACKTRACE);
constexpr auto use_color = static_cast<bool>(BOUT_USE_COLOR);
constexpr auto use_openmp = static_cast<bool>(BOUT_USE_OPENMP);
constexpr auto use_output_debug = static_cast<bool>(BOUT_USE_OUTPUT_DEBUG);
constexpr auto use_sigfpe = static_cast<bool>(BOUT_USE_SIGFPE);
constexpr auto use_signal = static_cast<bool>(BOUT_USE_SIGNAL);
constexpr auto use_track = static_cast<bool>(BOUT_USE_TRACK);
constexpr auto use_metric_3d = static_cast<bool>(BOUT_USE_METRIC_3D);

} // namespace build
} // namespace bout

#undef STRINGIFY1
#undef STRINGIFY

#endif // BOUT_BUILD_OPTIONS_HXX
