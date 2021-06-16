/*
 */
#include <bout/physicsmodel.hxx>
#include <bout.hxx>

#include <string>

std::string getCustomOutputName(Options& options) {
  const std::string default_filename = fmt::format("{}.{}", options.name(), "dmp");
  return fmt::format("{}/{}.{}.nc", options["datadir"].withDefault("data"),
                     options["file"].withDefault(default_filename), BoutComm::rank());
}

// Monitor to write out 1d Data
class Monitor1dDump : public Monitor {
public:
  Monitor1dDump(BoutReal timestep, std::string section_name)
      : Monitor(timestep),
        output_file(getCustomOutputName(Options::root()[section_name]),
                    Options::root()[section_name]["append"].withDefault(false)
                    ? bout::OptionsNetCDF::FileMode::append
                    : bout::OptionsNetCDF::FileMode::replace) {}
  int call(Solver*, BoutReal _time, int, int) override {
    Options output;
    output["t_array"].assignRepeat(_time);
    for (const auto& item: dump.getData()) {
      bout::utils::visit(bout::OptionsConversionVisitor{output, item.name}, item.value);
      if (item.repeat) {
        output[item.name].attributes["time_dimension"] = "t";
      }
    }
    output_file.write(output);

    return 0;
  }
  void add(BoutReal& data, const std::string& name) {
    dump.add(data, name, true);
  }

private:
  bout::DataFileFacade dump;
  bout::OptionsNetCDF output_file;
};

/// An example of using multiple monitors on different timescales
class MonitorExample : public PhysicsModel {
protected:
  int init(bool) {
    SOLVE_FOR2(n, T);

    // Our monitor writes out data every 100th of a time unit
    // note that this is independent of time_step.
    // Especially if the monitor influences the physics, and isn't
    // just used to output data, this is probably what you want.
    probes = bout::utils::make_unique<Monitor1dDump>(.01, "probes");

    // In case the monitor should be relative to the timestep, the
    // timestep needs to be read first:
    const BoutReal timestep = Options::root()["timestep"].withDefault(-1);

    // There is no 'slow' section in BOUT.inp, therefore it will write
    // to data/slow.dmp.0.nc
    slow = bout::utils::make_unique<Monitor1dDump>(timestep * 2, "slow");

    // Now we can add the monitors
    solver->addMonitor(probes.get());
    solver->addMonitor(slow.get());

    // The derived monitor Monitor1dData can dump 1d data, which we
    // can now add:
    probes->add(n(mesh->xstart, mesh->ystart, 0), "n_up");
    probes->add(T(mesh->xstart, mesh->ystart, 0), "T_up");

    // Add the corner value. T is already present in BOUT.dmp - but
    // as it is a different file, this doesn't cause issues
    slow->add(T(mesh->xstart, mesh->ystart, 0), "T");
    return 0;
  }

  int rhs(BoutReal) {
    ddt(n) = -T;
    ddt(T) = -n;
    return 0;
  }

private:
  Field3D n, T;

  std::unique_ptr<Monitor1dDump> probes{nullptr}, slow{nullptr};
};

BOUTMAIN(MonitorExample);
