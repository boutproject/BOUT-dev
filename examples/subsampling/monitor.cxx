/*
 */
#include <bout/physicsmodel.hxx>
#include <bout.hxx>

// Simplify the datafile with sane defaults
// Note that you should never ever copy the datafile.
// The constructor takes a pointer to an options object, or the name
// of a section.
// Unlike the default Datafile, SimpleDatafile always ends up in the
// correct folder, and should throw if it possible to detect the error.
class SimpleDatafile : public Datafile {
public:
  SimpleDatafile(std::string section)
      : SimpleDatafile(Options::root()[section], section){};
  SimpleDatafile(Options& ops, std::string section = "Default") : Datafile(&ops) {

    // Open a file for the output
    std::string datadir = "data";
    if (ops.isSet("path")) {
      datadir = ops["path"].as<std::string>();
    } else {
      datadir = Options::root()["datadir"].withDefault(datadir);
    }

    const std::string default_filename = section + ".dmp";
    const std::string file = ops["file"].withDefault(default_filename);

    bool append;
    if (ops.isSet("append")) {
      append = ops["append"].withDefault(false);
    } else {
      append = Options::root()["append"].withDefault(false);
    }

    const std::string dump_ext = "nc";

    if (append) {
      if (!this->opena("{:s}/{:s}.{:s}", datadir, file, dump_ext))
        throw BoutException("Failed to open file for appending!");
    } else {
      if (!this->openw("{:s}/{:s}.{:s}", datadir, file, dump_ext))
        throw BoutException("Failed to open file for writing!");
    }
  }
};

// Monitor to write out 1d Data
class Monitor1dDump : public Monitor {
public:
  Monitor1dDump(BoutReal timestep, std::string section_name)
      : Monitor(timestep), dump(bout::utils::make_unique<SimpleDatafile>(section_name)) {
    dump->add(time, "t_array", true);
  };
  int call(Solver*, BoutReal _time, int, int) override {
    time = _time;
    dump->write();
    return 0;
  }
  void add(BoutReal& data, const std::string& name) {
    dump->add(data, name.c_str(), true);
  }

private:
  BoutReal time;
  std::unique_ptr<Datafile> dump{nullptr};
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
