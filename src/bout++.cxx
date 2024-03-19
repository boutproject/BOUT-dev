/**************************************************************************
 *
 * Main BOUT++ functions
 * Adapted from the BOUT code by B.Dudson, University of York, Oct 2007
 *
 **************************************************************************
 * Copyright 2010-2023 BOUT++ contributors
 *
 * Contact Ben Dudson, dudson2@llnl.gov
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

#include "bout/build_config.hxx"

const char DEFAULT_DIR[] = "data";

#define GLOBALORIGIN

#include "bout++-time.hxx"
#include "bout/boundary_factory.hxx"
#include "bout/boutcomm.hxx"
#include "bout/boutexception.hxx"
#include "bout/coordinates_accessor.hxx"
#include "bout/hyprelib.hxx"
#include "bout/interpolation_xz.hxx"
#include "bout/interpolation_z.hxx"
#include "bout/invert/laplacexz.hxx"
#include "bout/invert_laplace.hxx"
#include "bout/invert_parderiv.hxx"
#include "bout/mpi_wrapper.hxx"
#include "bout/msg_stack.hxx"
#include "bout/openmpwrap.hxx"
#include "bout/optionsreader.hxx"
#include "bout/output.hxx"
#include "bout/petsclib.hxx"
#include "bout/revision.hxx"
#include "bout/rkscheme.hxx"
#include "bout/slepclib.hxx"
#include "bout/solver.hxx"
#include "bout/sys/timer.hxx"
#include "bout/version.hxx"

#define BOUT_NO_USING_NAMESPACE_BOUTGLOBALS
#include "bout/bout.hxx"
#undef BOUT_NO_USING_NAMESPACE_BOUTGLOBALS

#if BOUT_HAS_ADIOS
#include "bout/adios_object.hxx"
#endif

#include <fmt/format.h>

#include <csignal>
#include <ctime>
#include <string>
#include <vector>

#include <sys/stat.h>

// Value passed at compile time
// Used for MD5SUM, BOUT_LOCALE_PATH, and REVISION
#define BUILDFLAG1_(x) #x
#define BUILDFLAG(x) BUILDFLAG1_(x)

#define INDIRECT1_BOUTMAIN(a) #a
#define INDIRECT0_BOUTMAIN(...) INDIRECT1_BOUTMAIN(#__VA_ARGS__)
#define STRINGIFY(a) INDIRECT0_BOUTMAIN(a)

// Define S_ISDIR if not defined by system headers (that is, MSVC)
// Taken from https://github.com/curl/curl/blob/e59540139a398dc70fde6aec487b19c5085105af/lib/curl_setup.h#L748-L751
#if !defined(S_ISDIR) && defined(S_IFMT) && defined(S_IFDIR)
#define S_ISDIR(m) (((m)&S_IFMT) == S_IFDIR)
#endif

#ifdef _MSC_VER
#include <windows.h>
inline auto getpid() -> int { return GetCurrentProcessId(); }
#else
// POSIX headers
#include <unistd.h>
#endif

#if BOUT_USE_SIGFPE
#include <fenv.h>
#endif

using std::string;

BoutReal simtime{0.0};
int iteration{0};
bool user_requested_exit = false;

void bout_signal_handler(int sig);   // Handles signals
std::string time_to_hms(BoutReal t); // Converts to h:mm:ss.s format
char get_spin();                     // Produces a spinning bar

// Return the string "enabled" or "disabled"
namespace {
constexpr auto is_enabled(bool enabled) { return enabled ? "enabled" : "disabled"; }
} // namespace

/*!
  Initialise BOUT++

  Inputs
  ------

  The command-line arguments argc and argv are passed by
  reference, and pointers to these will be stored in various
  places in BOUT++.

  Outputs
  -------

  Any non-zero return value should halt the simulation. If the return value is
  less than zero, the exit status from BOUT++ is 0, otherwise it is the return
  value of BoutInitialise.

 */
int BoutInitialise(int& argc, char**& argv) {

  using namespace bout::experimental;

  setupSignalHandler(bout_signal_handler);

  setupGetText();

  CommandLineArgs args;
  try {
    args = parseCommandLineArgs(argc, argv);
  } catch (const BoutException& e) {
    output_error << _("Bad command line arguments:\n") << e.what() << std::endl;
    return 1;
  }

  try {
    checkDataDirectoryIsAccessible(args.data_dir);

    // Set the command-line arguments
    SlepcLib::setArgs(argc, argv); // SLEPc initialisation
    PetscLib::setArgs(argc, argv); // PETSc initialisation
    Solver::setArgs(argc, argv);   // Solver initialisation
    BoutComm::setArgs(argc, argv); // MPI initialisation

    const int MYPE = BoutComm::rank();

    setupBoutLogColor(args.color_output, MYPE);

    setupOutput(args.data_dir, args.log_file, args.verbosity, MYPE);

    savePIDtoFile(args.data_dir, MYPE);

#if BOUT_HAS_ADIOS
    bout::ADIOSInit(BoutComm::get());
#endif

    // Print the different parts of the startup info
    printStartupHeader(MYPE, BoutComm::size());
    printCompileTimeOptions();
    printCommandLineArguments(args.original_argv);

    // Load settings file
    OptionsReader* reader = OptionsReader::getInstance();
    // Ideally we'd use the long options for `datadir` and
    // `optionfile` here, but we'd need to call parseCommandLine
    // _first_ in order to do that and set the source, etc., but we
    // need to call that _second_ in order to override the input file
    reader->read(Options::getRoot(), "{}/{}", args.data_dir, args.opt_file);

    // Get options override from command-line
    reader->parseCommandLine(Options::getRoot(), args.argv);

    // Get the variables back out so they count as having been used
    // when checking for unused options. They normally _do_ get used,
    // but it's possible that only happens in BoutFinalise, which is
    // too late for that check.
    const auto datadir = Options::root()["datadir"].withDefault<std::string>(DEFAULT_DIR);
    [[maybe_unused]] const auto optionfile =
        Options::root()["optionfile"].withDefault<std::string>(args.opt_file);
    const auto settingsfile =
        Options::root()["settingsfile"].withDefault<std::string>(args.set_file);

    setRunStartInfo(Options::root());

    if (MYPE == 0) {
      writeSettingsFile(Options::root(), datadir, settingsfile);
    }

    bout::globals::mpi = new MpiWrapper();

    // Create the mesh
    bout::globals::mesh = Mesh::create();
    // Load from sources. Required for Field initialisation
    bout::globals::mesh->load();

    // time_report options are used in BoutFinalise, i.e. after we
    // check for unused options
    Options::root()["time_report"].setConditionallyUsed();

  } catch (const BoutException& e) {
    output_error.write(_("Error encountered during initialisation: {:s}\n"), e.what());
    throw;
  }

  return 0;
}

namespace bout {
namespace experimental {
void setupSignalHandler(SignalHandler signal_handler) {
#if BOUT_USE_SIGNAL
  std::signal(SIGSEGV, signal_handler);
  std::signal(SIGBUS, signal_handler);
#endif
#if BOUT_USE_SIGFPE
  std::signal(SIGFPE, signal_handler);
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

#ifndef _MSC_VER
  /// Trap SIGUSR1 to allow a clean exit after next write
  std::signal(SIGUSR1, signal_handler);
#endif
}

// This is currently just an alias to the existing handler
void defaultSignalHandler(int sig) { bout_signal_handler(sig); }

void setupGetText() {
#if BOUT_HAS_GETTEXT
  // Setting the i18n environment
  //
  // For libraries:
  // https://www.gnu.org/software/gettext/manual/html_node/Libraries.html
  //
  try {
    // Note: Would like to use std::locale::global
    //    std::locale::global(std::locale(""));
    // but the Numeric aspect causes problems parsing input strings
    //
    // Note: Since BOUT++ is a library, it shouldn't really call setlocale;
    //       that should be part of main().
    std::setlocale(LC_ALL, "");
    std::setlocale(LC_NUMERIC, "C");

    bindtextdomain(GETTEXT_PACKAGE, BUILDFLAG(BOUT_LOCALE_PATH));
  } catch (const std::runtime_error& e) {
#if 1
    fmt::print(
        stderr,
        FMT_STRING(
            "WARNING: Could not set locale. Check the LANG environment variable "
            "(get available values by running 'locale -a'). If LANG is correct, there "
            "may be "
            "a problem with the BOUT_LOCALE_PATH={:s} that BOUT++ was compiled with.\n"),
        BUILDFLAG(BOUT_LOCALE_PATH));
#endif
  }
#endif // BOUT_HAS_GETTEXT
}

/// Print all of the available implementations for a given `Factory`
template <class Factory>
[[noreturn]] void printAvailableImplementations() {
  const auto factory = Factory::getInstance();

  for (const auto& implementation : factory.listAvailable()) {
    std::cout << implementation << "\n";
  }
  auto unavailable = factory.listUnavailableReasons();
  if (not unavailable.empty()) {
    std::cout << fmt::format("\nThe following {}s are currently unavailable:\n",
                             Factory::type_name);
    for (const auto& implementation : unavailable) {
      std::cout << implementation << "\n";
    }
  }
  std::exit(EXIT_SUCCESS);
}

/// Print all the Options used in constructing \p type
template <class Factory>
[[noreturn]] void printTypeOptions(const std::string& type) {
  const auto factory = Factory::getInstance();

  // Make sure all the type construction is quiet
  output.disable();
  output_error.disable();
  output_warn.disable();
  output_progress.disable();
  output_info.disable();
  output_verbose.disable();
  output_debug.disable();

  // There are some non-optional arguments to mesh we'll need to
  // supply if we don't have an input file
  Options& mesh_options = Options::root()["mesh"];
  mesh_options["MXG"] = 1;
  mesh_options["MYG"] = 1;
  mesh_options["nx"] = 4;
  mesh_options["ny"] = 4;
  mesh_options["nz"] = 4;

  // We might need a global mesh for some types, so best make one
  bout::globals::mpi = new MpiWrapper();
  bout::globals::mesh = Mesh::create();
  bout::globals::mesh->load();

  // An empty Options that we'll later check for used values
  Options help_options;

  // Most likely failure is typo in type name, so we definitely want
  // to print that
  try {
    factory.create(type, &help_options[Factory::section_name]);
  } catch (const BoutException& error) {
    std::cout << error.what() << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Now we can print all the options used in constructing our
  // type. Note that this does require all the options are used in the
  // constructor, and not in a `init` method or similar
  std::cout << fmt::format("Input options for {} '{}':\n\n", Factory::type_name, type);
  std::cout << fmt::format("{:id}\n", help_options);
  std::exit(EXIT_SUCCESS);
}

/// Handle the command line options for the given `Factory`: listing
/// all types, and used options for a given type
template <class Factory>
void handleFactoryHelp(const std::string& current_arg, int i, int argc, char** argv,
                       bool plural_needs_e = false) {
  const auto name = lowercase(Factory::type_name);

  const auto list_arg = fmt::format("--list-{}{}s", name, plural_needs_e ? "e" : "");
  const auto help_arg = fmt::format("--help-{}", name);

  if (current_arg == list_arg) {
    printAvailableImplementations<Factory>();
  }

  if (current_arg == help_arg) {
    if (i + 1 >= argc) {
      throw BoutException(_("Usage is {} {} <name>\n"), argv[0], help_arg);
    }
    printTypeOptions<Factory>(argv[i + 1]);
  }
}

auto parseCommandLineArgs(int argc, char** argv) -> CommandLineArgs {
  /// NB: "restart" and "append" are now caught by options
  /// Check for help flag separately
  for (int i = 1; i < argc; i++) {
    const std::string current_arg{argv[i]};
    if (current_arg == "-h" || current_arg == "--help") {
      // Print help message -- note this will be displayed once per processor as we've not
      // started MPI yet.
      output.write(_("Usage: {:s} [-d <data directory>] [-f <options filename>] [restart "
                     "[append]] [VAR=VALUE]\n"),
                   argv[0]);
      output.write(
          _("\n"
            "  -d <data directory>\t\tLook in <data directory> for input/output files\n"
            "  -f <options filename>\t\tUse OPTIONS given in <options filename>\n"
            "  -o <settings filename>\tSave used OPTIONS given to <options filename>\n"
            "  -l, --log <log filename>\tPrint log to <log filename>\n"
            "  -v, --verbose\t\t\tIncrease verbosity\n"
            "  -q, --quiet\t\t\tDecrease verbosity\n"));
#if BOUT_USE_COLOR
      output.write(_("  -c, --color\t\t\tColor output using bout-log-color\n"));
#endif
      output.write(
          _("  --print-config\t\tPrint the compile-time configuration\n"
            "  --list-solvers\t\tList the available time solvers\n"
            "  --help-solver <solver>\tPrint help for the given time solver\n"
            "  --list-laplacians\t\tList the available Laplacian inversion solvers\n"
            "  --help-laplacian <laplacian>\tPrint help for the given Laplacian "
            "inversion solver\n"
            "  --list-laplacexz\t\tList the available LaplaceXZ inversion solvers\n"
            "  --help-laplacexz <laplacexz>\tPrint help for the given LaplaceXZ "
            "inversion solver\n"
            "  --list-invertpars\t\tList the available InvertPar solvers\n"
            "  --help-invertpar <invertpar>\tPrint help for the given InvertPar solver\n"
            "  --list-rkschemes\t\tList the available Runge-Kutta schemes\n"
            "  --help-rkscheme <rkscheme>\tPrint help for the given Runge-Kutta scheme\n"
            "  --list-meshes\t\t\tList the available Meshes\n"
            "  --help-mesh <mesh>\t\tPrint help for the given Mesh\n"
            "  --list-xzinterpolations\tList the available XZInterpolations\n"
            "  --help-xzinterpolation <xzinterpolation>\tPrint help for the given "
            "XZInterpolation\n"
            "  --list-zinterpolations\tList the available ZInterpolations\n"
            "  --help-zinterpolation <zinterpolation>\tPrint help for the given "
            "ZInterpolation\n"
            "  -h, --help\t\t\tThis message\n"
            "  restart [append]\t\tRestart the simulation. If append is specified, "
            "append to the existing output files, otherwise overwrite them\n"
            "  VAR=VALUE\t\t\tSpecify a VALUE for input parameter VAR\n"
            "\nFor all possible input parameters, see the user manual and/or the "
            "physics model source (e.g. {:s}.cxx)\n"),
          argv[0]);

      std::exit(EXIT_SUCCESS);
    }
    if (current_arg == "--print-config") {
      printCompileTimeOptions();
      std::exit(EXIT_SUCCESS);
    }
    handleFactoryHelp<SolverFactory>(current_arg, i, argc, argv);
    handleFactoryHelp<LaplaceFactory>(current_arg, i, argc, argv);
    handleFactoryHelp<LaplaceXZFactory>(current_arg, i, argc, argv);
    handleFactoryHelp<InvertParFactory>(current_arg, i, argc, argv);
    handleFactoryHelp<RKSchemeFactory>(current_arg, i, argc, argv);
    handleFactoryHelp<MeshFactory>(current_arg, i, argc, argv, true);
    handleFactoryHelp<XZInterpolationFactory>(current_arg, i, argc, argv);
    handleFactoryHelp<ZInterpolationFactory>(current_arg, i, argc, argv);
  }

  CommandLineArgs args;
  args.original_argv.reserve(argc);
  std::copy(argv, argv + argc, std::back_inserter(args.original_argv));
  args.argv = args.original_argv;

  // Parse single-letter and non-Options arguments, expanding the
  // single-letter arguments into their canonical names, and removing
  // the non-Options arguments
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "-d") {
      // Set data directory
      if (i + 1 >= argc) {
        throw BoutException(_("Usage is {:s} -d <data directory>\n"), argv[0]);
      }

      args.data_dir = argv[++i];
      args.argv[i - 1] = "datadir=";

    } else if (string(argv[i]) == "-f") {
      // Set options file
      if (i + 1 >= argc) {
        throw BoutException(_("Usage is {:s} -f <options filename>\n"), argv[0]);
      }

      args.opt_file = argv[++i];
      args.argv[i - 1] = "optionfile=";

    } else if (string(argv[i]) == "-o") {
      // Set options file
      if (i + 1 >= argc) {
        throw BoutException(_("Usage is {:s} -o <settings filename>\n"), argv[0]);
      }

      args.set_file = argv[++i];
      args.argv[i - 1] = "settingsfile=";

    } else if ((string(argv[i]) == "-l") || (string(argv[i]) == "--log")) {
      // Set log file
      if (i + 1 >= argc) {
        throw BoutException(_("Usage is {:s} -l <log filename>\n"), argv[0]);
      }

      args.log_file = argv[++i];
      args.argv[i - 1] = "logfile=";

    } else if ((string(argv[i]) == "-v") || (string(argv[i]) == "--verbose")) {
      args.verbosity++;
      args.argv[i] = "";

    } else if ((string(argv[i]) == "-q") || (string(argv[i]) == "--quiet")) {
      args.verbosity--;
      args.argv[i] = "";

    } else if ((string(argv[i]) == "-c") || (string(argv[i]) == "--color")) {
      // Add color to the output by piping through bout-log-color
      // This is done after checking all command-line inputs
      // in case -c is set multiple times
      args.color_output = true;
      args.argv[i] = "";
    }
  }

  // Remove empty values from canonicalised arguments
  bout::utils::erase(args.argv, "");

  if (args.set_file == args.opt_file) {
    throw BoutException(
        _("Input and output file for settings must be different.\nProvide -o <settings "
          "file> to avoid this issue.\n"));
  }

  return args;
}

void checkDataDirectoryIsAccessible(const std::string& data_dir) {
  struct stat test;
  if (stat(data_dir.c_str(), &test) == 0) {
    if (!S_ISDIR(test.st_mode)) {
      throw BoutException(_("DataDir \"{:s}\" is not a directory\n"), data_dir);
    }
  } else {
    throw BoutException(_("DataDir \"{:s}\" does not exist or is not accessible\n"),
                        data_dir);
  }
}

void savePIDtoFile(const std::string& data_dir, int MYPE) {
  std::stringstream filename;
  filename << data_dir << "/.BOUT.pid." << MYPE;
  std::ofstream pid_file;
  pid_file.open(filename.str(), std::ios::out | std::ios::trunc);

  if (not pid_file.is_open()) {
    throw BoutException(_("Could not create PID file {:s}"), filename.str());
  }

  pid_file << getpid() << "\n";
  pid_file.close();
}

void printStartupHeader(int MYPE, int NPES) {
  output_progress.write(_("BOUT++ version {:s}\n"), bout::version::full);
  output_progress.write(_("Revision: {:s}\n"), bout::version::revision);
#ifdef MD5SUM
  output_progress.write("MD5 checksum: {:s}\n", BUILDFLAG(MD5SUM));
#endif
  output_progress.write(_("Code compiled on {:s} at {:s}\n\n"), boutcompiledate,
                        boutcompiletime);
  output_info.write("B.Dudson (University of York), M.Umansky (LLNL) 2007\n");
  output_info.write("Based on BOUT by Xueqiao Xu, 1999\n\n");

  output_info.write(_("Processor number: {:d} of {:d}\n\n"), MYPE, NPES);

  output_info.write("pid: {:d}\n\n", getpid());
}

void printCompileTimeOptions() {
  output_info.write(_("Compile-time options:\n"));

  using namespace bout::build;

  output_info.write(_("\tRuntime error checking {}"), is_enabled(check_level > 0));
  if (check_level > 0) {
    output_info.write(_(", level {}"), check_level);
  }
  output_info.write("\n");

#ifdef PNCDF
  output_info.write(_("\tParallel NetCDF support enabled\n"));
#else
  output_info.write(_("\tParallel NetCDF support disabled\n"));
#endif

  output_info.write(_("\tMetrics mode is {}\n"), use_metric_3d ? "3D" : "2D");
  output_info.write(_("\tFFT support {}\n"), is_enabled(has_fftw));
  output_info.write(_("\tNatural language support {}\n"), is_enabled(has_gettext));
  output_info.write(_("\tLAPACK support {}\n"), is_enabled(has_lapack));
  // Horrible nested ternary to set this at compile time
  constexpr auto netcdf_flavour =
      has_netcdf ? (has_legacy_netcdf ? " (Legacy)" : " (NetCDF4)") : "";
  output_info.write(_("\tNetCDF support {}{}\n"), is_enabled(has_netcdf), netcdf_flavour);
  output_info.write(_("\tADIOS support {}\n"), is_enabled(has_adios));
  output_info.write(_("\tPETSc support {}\n"), is_enabled(has_petsc));
  output_info.write(_("\tPretty function name support {}\n"),
                    is_enabled(has_pretty_function));
  output_info.write(_("\tPVODE support {}\n"), is_enabled(has_pvode));
  output_info.write(_("\tScore-P support {}\n"), is_enabled(has_scorep));
  output_info.write(_("\tSLEPc support {}\n"), is_enabled(has_slepc));
  output_info.write(_("\tSUNDIALS support {}\n"), is_enabled(has_sundials));
  output_info.write(_("\tBacktrace in exceptions {}\n"), is_enabled(use_backtrace));
  output_info.write(_("\tColour in logs {}\n"), is_enabled(use_color));
  output_info.write(_("\tOpenMP parallelisation {}, using {} threads\n"),
                    is_enabled(use_openmp), omp_get_max_threads());
  output_info.write(_("\tExtra debug output {}\n"), is_enabled(use_output_debug));
  output_info.write(_("\tFloating-point exceptions {}\n"), is_enabled(use_sigfpe));
  output_info.write(_("\tSignal handling support {}\n"), is_enabled(use_signal));
  output_info.write(_("\tField name tracking {}\n"), is_enabled(use_track));
  output_info.write(_("\tMessage stack {}\n"), is_enabled(use_msgstack));

  // The stringify is needed here as BOUT_FLAGS_STRING may already contain quoted strings
  // which could cause problems (e.g. terminate strings).
  output_info.write(_("\tCompiled with flags : {:s}\n"), STRINGIFY(BOUT_FLAGS_STRING));
}

void printCommandLineArguments(const std::vector<std::string>& original_argv) {
  output_info.write(_("\tCommand line options for this run : "));
  for (auto& arg : original_argv) {
    output_info << arg << " ";
  }
  output_info.write("\n");
}

bool setupBoutLogColor(bool color_output, int MYPE) {
#if BOUT_USE_COLOR
  if (color_output && (MYPE == 0)) {
    // Color stdout by piping through bout-log-color script
    // Only done on processor 0, since this is the only processor which writes to stdout
    // This uses popen, fileno and dup2 functions, which are POSIX
    bool success = false;

    // Run bout-log-color through the shell. This should share stdout with BOUT++,
    // and read stdin from the pipe
    FILE* outpipe = popen("bout-log-color", "w");

    if (outpipe != nullptr) {
      // Valid pipe
      // Get the integer file descriptor
      int fno = fileno(outpipe);
      if (fno != -1) {
        // Valid file descriptor

        // Note: We can get to here if bout-log-color failed to run
        // This seems to cause code to fail later

        // Replace stdout with the pipe.
        int status = dup2(fno, STDOUT_FILENO);
        if (status != -1) {
          success = true;
        }
      }
    }
    if (!success) {
      // Failed . Probably not important enough to stop the simulation
      std::cerr << _("Could not run bout-log-color. Make sure it is in your PATH\n");
    }
    return success;
  }
#endif // BOUT_USE_COLOR
  return false;
}

void setupOutput(const std::string& data_dir, const std::string& log_file, int verbosity,
                 int MYPE) {
  {
    Output& output = *Output::getInstance();
    if (MYPE == 0) {
      output.enable(); // Enable writing to stdout
    } else {
      output.disable(); // No writing to stdout
    }
    /// Open an output file to echo everything to
    /// On processor 0 anything written to output will go to stdout and the file
    if (output.open("{:s}/{:s}.{:d}", data_dir, log_file, MYPE)) {
      throw BoutException(_("Could not open {:s}/{:s}.{:d} for writing"), data_dir,
                          log_file, MYPE);
    }
  }

  output_error.enable(verbosity > 0);
  output_warn.enable(verbosity > 1);
  output_progress.enable(verbosity > 2);
  output_info.enable(verbosity > 3);
  output_verbose.enable(verbosity > 4);
  // Only actually enabled if also compiled with ENABLE_OUTPUT_DEBUG
  output_debug.enable(verbosity > 5);

  // The backward-compatible output object same as output_progress
  output.enable(verbosity > 2);
}

void setRunStartInfo(Options& options) {
  auto& runinfo = options["run"];

  // Note: have to force value, since may already be set if a previously
  // output BOUT.settings file was used as input
  runinfo["version"].force(bout::version::full, "Output");
  runinfo["revision"].force(bout::version::revision, "Output");

  time_t start_time = time(nullptr);
  runinfo["started"].force(ctime(&start_time), "Output");
}

void setRunFinishInfo(Options& options) {
  time_t end_time = time(nullptr);
  options["run"]["finished"].force(ctime(&end_time), "Output");
}

void addBuildFlagsToOptions(Options& options) {
  output_progress << "Setting up output (experimental output) file\n";

  options["BOUT_VERSION"].force(bout::version::as_double);
  options["use_check_level"].force(bout::build::check_level);
  options["use_openmp_schedule"].force(bout::build::openmp_schedule);
  options["has_fftw"].force(bout::build::has_fftw);
  options["has_gettext"].force(bout::build::has_gettext);
  options["has_lapack"].force(bout::build::has_lapack);
  options["has_netcdf"].force(bout::build::has_netcdf);
  options["has_adios"].force(bout::build::has_adios);
  options["has_petsc"].force(bout::build::has_petsc);
  options["has_hypre"].force(bout::build::has_hypre);
  options["has_umpire"].force(bout::build::has_umpire);
  options["has_caliper"].force(bout::build::has_caliper);
  options["has_raja"].force(bout::build::has_raja);
  options["has_pretty_function"].force(bout::build::has_pretty_function);
  options["has_pvode"].force(bout::build::has_pvode);
  options["has_scorep"].force(bout::build::has_scorep);
  options["has_slepc"].force(bout::build::has_slepc);
  options["has_sundials"].force(bout::build::has_sundials);
  options["use_backtrace"].force(bout::build::use_backtrace);
  options["use_color"].force(bout::build::use_color);
  options["use_openmp"].force(bout::build::use_openmp);
  options["openmp_threads"].force(omp_get_max_threads());
  options["use_output_debug"].force(bout::build::use_output_debug);
  options["use_sigfpe"].force(bout::build::use_sigfpe);
  options["use_signal"].force(bout::build::use_signal);
  options["use_track"].force(bout::build::use_track);
  options["has_cuda"].force(bout::build::has_cuda);
  options["use_metric_3d"].force(bout::build::use_metric_3d);
  options["use_msgstack"].force(bout::build::use_msgstack);
}

void writeSettingsFile(Options& options, const std::string& data_dir,
                       const std::string& settings_file) {
  OptionsReader::getInstance()->write(&options, "{}/{}", data_dir, settings_file);
}

} // namespace experimental
} // namespace bout

int BoutFinalise(bool write_settings) {

  // Output the settings, showing which options were used
  // This overwrites the file written during initialisation
  if (write_settings) {
    try {
      using namespace bout::experimental;
      auto& options = Options::root();

      setRunFinishInfo(options);

      const auto data_dir = options["datadir"].withDefault(std::string{DEFAULT_DIR});
      const auto set_file = options["settingsfile"].withDefault("BOUT.settings");

      if (BoutComm::rank() == 0) {
        writeSettingsFile(options, data_dir, set_file);
      }
    } catch (const BoutException& e) {
      output_error << _("Error whilst writing settings") << e.what() << endl;
    }
  }

  if (Options::root()["time_report:show"].withDefault(false)) {
    output.write("\nTimer report \n\n");
    Timer::printTimeReport();
    output.write("\n");
  }

  // Delete the mesh
  delete bout::globals::mesh;

  // Make sure all processes have finished writing before exit
  bout::globals::mpi->MPI_Barrier(BoutComm::get());

  // Laplacian inversion
  Laplacian::cleanup();

  // Delete field memory
  Array<BoutReal>::cleanup();
  Array<dcomplex>::cleanup();
  Array<fcmplx>::cleanup();
  Array<int>::cleanup();
  Array<unsigned long>::cleanup();

  // Cleanup boundary factory
  BoundaryFactory::cleanup();

  CoordinatesAccessor::clear();

  // Cleanup timer
  Timer::cleanup();

  // Options tree
  Options::cleanup();
  OptionsReader::cleanup();

  // Call SlepcFinalize if not already called
  SlepcLib::cleanup();

  // Call PetscFinalize if not already called
  PetscLib::cleanup();

  // Call HYPER_Finalize if not already called
  bout::HypreLib::cleanup();

#if BOUT_HAS_ADIOS
  bout::ADIOSFinalize();
#endif

  // MPI communicator, including MPI_Finalize()
  BoutComm::cleanup();

  // Debugging message stack
  msg_stack.clear();

  // Delete the MPI wrapper
  delete bout::globals::mpi;

  return 0;
}

/*!*************************************************************************
 * SOLUTION MONITOR FUNCTION
 *
 * Called each timestep by the solver
 **************************************************************************/

BoutMonitor::BoutMonitor(BoutReal timestep) : BoutMonitor(timestep, Options::root()) {}

BoutMonitor::BoutMonitor(BoutReal timestep, Options& options)
    : Monitor(timestep),
      wall_limit(options["wall_limit"]
                     .doc(_("Wall time limit in hours. By default (< 0), no limit"))
                     .withDefault(-1.0)
                 * 60. * 60.),
      stop_check(options["stopCheck"]
                     .doc(_("Check if a file exists, and exit if it does."))
                     .withDefault(false)),
      stop_check_name(
          fmt::format("{}/{}", Options::root()["datadir"].withDefault(DEFAULT_DIR),
                      options["stopCheckName"]
                          .doc(_("Name of file whose existence triggers a stop"))
                          .withDefault("BOUT.stop"))) {}

int BoutMonitor::call(Solver* solver, BoutReal t, [[maybe_unused]] int iter, int NOUT) {
  TRACE("BoutMonitor::call({:e}, {:d}, {:d})", t, iter, NOUT);

  // Increment Solver's iteration counter, and set the global `iteration`
  iteration = solver->incrementIterationCounter();

  // Set the global variables. This is done because they need to be
  // written to the output file before the first step (initial condition)
  simtime = t;

  /// Collect timing information
  run_data.wtime = Timer::resetTime("run");
  run_data.ncalls = solver->resetRHSCounter();
  run_data.ncalls_e = solver->resetRHSCounter_e();
  run_data.ncalls_i = solver->resetRHSCounter_i();

  const bool output_split = solver->splitOperator();
  run_data.wtime_rhs = Timer::resetTime("rhs");
  run_data.wtime_invert = Timer::resetTime("invert");
  // Time spent communicating (part of RHS)
  run_data.wtime_comms = Timer::resetTime("comms");
  // Time spend on I/O
  run_data.wtime_io = Timer::resetTime("io");

  run_data.calculateDerivedMetrics();

  output_progress.print("\r"); // Only goes to screen

  // First time the monitor has been called
  static bool first_time = true;
  if (first_time) {

    // Record the starting time
    mpi_start_time = bout::globals::mpi->MPI_Wtime() - run_data.wtime;

    first_time = false;

    // Print the column header for timing info
    if (!output_split) {
      output_progress.write(_("Sim Time  |  RHS evals  | Wall Time |  Calc    Inv   Comm "
                              "   I/O   SOLVER\n\n"));
    } else {
      output_progress.write(_("Sim Time  |  RHS_e evals  | RHS_I evals  | Wall Time |  "
                              "Calc    Inv   Comm    I/O   SOLVER\n\n"));
    }
  }

  run_data.writeProgress(simtime, output_split);

  // This bit only to screen, not log file

  run_data.t_elapsed = bout::globals::mpi->MPI_Wtime() - mpi_start_time;

  output_progress.print("{:c}  Step {:d} of {:d}. Elapsed {:s}", get_spin(), iteration,
                        NOUT, time_to_hms(run_data.t_elapsed));
  output_progress.print(
      " ETA {:s}",
      time_to_hms(run_data.wtime * static_cast<BoutReal>(NOUT - iteration - 2)));

  // Write dump file
  Options run_data_output;
  run_data.outputVars(run_data_output);
  solver->writeToModelOutputFile(run_data_output);

  if (wall_limit > 0.0) {
    // Check if enough time left

    BoutReal t_remain = mpi_start_time + wall_limit - bout::globals::mpi->MPI_Wtime();
    if (t_remain < run_data.wtime * 2) {
      // Less than 2 time-steps left
      output_warn.write(_("Only {:e} seconds ({:.2f} steps) left. Quitting\n"), t_remain,
                        t_remain / run_data.wtime);
      user_requested_exit = true;
    } else {
      output_progress.print(" Wall {:s}", time_to_hms(t_remain));
    }
  }

  // Check if the user has created the stop file and if so trigger an exit
  if (stop_check) {
    std::ifstream f(stop_check_name);
    if (f.good()) {
      output.write("\nStop file {} exists -- triggering exit\n", stop_check_name);
      user_requested_exit = true;
    }
  }

  return 0;
}

/**************************************************************************
 * Global error handling
 **************************************************************************/

/// Signal handler - handles all signals
void bout_signal_handler(int sig) {
  // Set signal handler back to default to prevent possible infinite loop
  signal(SIGSEGV, SIG_DFL);
  // print number of process to stderr, so the user knows which log to check
  fmt::print(stderr, FMT_STRING("\nSighandler called on process {:d} with sig {:d}\n"),
             BoutComm::rank(), sig);

  switch (sig) {
  case SIGSEGV:
    throw BoutException("\n****** SEGMENTATION FAULT CAUGHT ******\n\n");
    break;
  case SIGFPE:
    throw BoutException("\n****** Floating Point Exception "
                        "(FPE) caught ******\n\n");
    break;
  case SIGINT:
    throw BoutException("\n****** SigInt caught ******\n\n");
    break;
#ifndef _MSC_VER
  case SIGKILL:
    throw BoutException("\n****** SigKill caught ******\n\n");
    break;
  case SIGUSR1:
    user_requested_exit = true;
    break;
#endif
  default:
    throw BoutException("\n****** Signal {:d}  caught ******\n\n", sig);
    break;
  }
}

/**************************************************************************
 * Utilities
 **************************************************************************/

/// Write a time in h:mm:ss.s format
std::string time_to_hms(BoutReal t) {
  int h, m;

  h = static_cast<int>(t / 3600);
  t -= 3600. * static_cast<BoutReal>(h);
  m = static_cast<int>(t / 60);
  t -= 60 * static_cast<BoutReal>(m);

  return fmt::format(FMT_STRING("{:d}:{:02d}:{:04.1f}"), h, m, t);
}

/// Produce a spinning bar character
char get_spin() {
  static int i = 0;
  char c = '|'; // Doesn't need to be assigned; squash warning

  switch (i) {
  case 0:
    c = '|';
    break;
  case 1:
    c = '/';
    break;
  case 2:
    c = '-';
    break;
  case 3:
    c = '\\';
    break;
  }
  i = (i + 1) % 4;
  return c;
}

/**************************************************************************
 * Functions for writing run information
 **************************************************************************/

void RunMetrics::outputVars(Options& output_options) const {
  Timer time("io");
  output_options["wall_time"].assignRepeat(t_elapsed, "t", true, "Output");
  output_options["wtime"].assignRepeat(wtime, "t", true, "Output");
  output_options["ncalls"].assignRepeat(ncalls, "t", true, "Output");
  output_options["ncalls_e"].assignRepeat(ncalls_e, "t", true, "Output");
  output_options["ncalls_i"].assignRepeat(ncalls_i, "t", true, "Output");
  output_options["wtime_rhs"].assignRepeat(wtime_rhs, "t", true, "Output");
  output_options["wtime_invert"].assignRepeat(wtime_invert, "t", true, "Output");
  output_options["wtime_comms"].assignRepeat(wtime_comms, "t", true, "Output");
  output_options["wtime_io"].assignRepeat(wtime_io, "t", true, "Output");
  output_options["wtime_per_rhs"].assignRepeat(wtime_per_rhs, "t", true, "Output");
  output_options["wtime_per_rhs_e"].assignRepeat(wtime_per_rhs_e, "t", true, "Output");
  output_options["wtime_per_rhs_i"].assignRepeat(wtime_per_rhs_i, "t", true, "Output");
}

void RunMetrics::calculateDerivedMetrics() {
  // Terrible hack avoid divide-by-zero, needed because SLEPc solver
  // doesn't call `run_rhs` which increments `ncalls`. Better fix is
  // change `Solver::addMonitor` API to take a name so that we can
  // replace `BoutMonitor` with a different implementation. Currently
  // not possible because `Solver::removeMonitor` needs the pointer to
  // the specific instance
  if (ncalls == 0) {
    return;
  }
  wtime_per_rhs = wtime / ncalls;
  wtime_per_rhs_e = wtime / ncalls_e;
  wtime_per_rhs_i = wtime / ncalls_i;
}

void RunMetrics::writeProgress(BoutReal simtime, bool output_split) {
  if (!output_split) {
    output_progress.write(
        "{:.3e}      {:5d}       {:.2e}   {:5.1f}  {:5.1f}  {:5.1f}  {:5.1f}  {:5.1f}\n",
        simtime, ncalls, wtime, 100. * (wtime_rhs - wtime_comms - wtime_invert) / wtime,
        100. * wtime_invert / wtime,                    // Inversions
        100. * wtime_comms / wtime,                     // Communications
        100. * wtime_io / wtime,                        // I/O
        100. * (wtime - wtime_io - wtime_rhs) / wtime); // Everything else

  } else {
    output_progress.write("{:.3e}      {:5d}            {:5d}       {:.2e}   {:5.1f}  "
                          "{:5.1f}  {:5.1f}  {:5.1f}  {:5.1f}\n",
                          simtime, ncalls_e, ncalls_i, wtime,
                          100. * (wtime_rhs - wtime_comms - wtime_invert) / wtime,
                          100. * wtime_invert / wtime, // Inversions
                          100. * wtime_comms / wtime,  // Communications
                          100. * wtime_io / wtime,     // I/O
                          100. * (wtime - wtime_io - wtime_rhs)
                              / wtime); // Everything else
  }
}
