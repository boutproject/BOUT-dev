/*
 */
#include <bout/physicsmodel.hxx>
#include <bout.hxx>

// simplify the datafile with sane defaults
// Note that you should never ever copy the datafile ...
// The constructor takes a pointer to an options object, or the name
// of a section.
// Unlike the default Datafile, SimpleDatafile always ends up in the
// correct folder, and should throw if it possible to detect the error.
class SimpleDatafile: public Datafile{
public:
  SimpleDatafile(std::string section)
      : SimpleDatafile(Options::root()[section], section){};
  SimpleDatafile(Options& ops, std::string section = "Default") : Datafile(&ops) {
    // Open a file for the output
    std::string datadir = "data";
    if (ops.isSet("path")) {
      datadir = ops["path"].withDefault(datadir); // default never needed
    } else {
      datadir = Options::root()["datadir"].withDefault(
          datadir); // I need to know data is default :(
    }
    std::string file;
    file = section+".dmp";
    file = ops["file"].withDefault(file);
    bool append;
    if (ops.isSet("append")) {
      append = ops["append"].withDefault(false);
    } else {
      append = Options::root()["append"].withDefault(
          false); // I hope that is the correct default
    }
    std::string dump_ext="nc"; // bad style, but I only use nc

    if(append) {
      if (!this->opena("%s/%s.%s", datadir.c_str(), file.c_str(), dump_ext.c_str()))
        throw BoutException("Failed to open file for appending!");
      //output.write("opend succesfully for appending\n");
    }else {
      if (!this->openw("%s/%s.%s", datadir.c_str(), file.c_str(), dump_ext.c_str()))
        throw BoutException("Failed to open file for writing!");
      //output.write("opend succesfully for writing\n");
    }
  }
};

// Monitor to write out 1d Data
class Monitor1dDump:public Monitor{
public:
  Monitor1dDump(BoutReal timestep,std::string section_name)
    :Monitor(timestep),
    dump(new SimpleDatafile(section_name))
  {
    dump->add(time,"t_array",true);
  };
  int call(Solver * solver,double _time,int i1,int i2) override{
    time=_time;
    dump->write(); // this should throw if it doesn't work
    return 0;
  }
  void add(BoutReal &data,std::string name){
    dump->add(data,name.c_str(),true);
  }
private:
  BoutReal time;
  Datafile * dump;
};



class MonitorExample : public PhysicsModel {
protected:
  int init(bool restarting) {
    SOLVE_FOR2(n,T);
    // our monitor writes out data every 100th of a time unit
    // note that this is independent of time_step.
    // Especially if the monitor influences the physics, and isn't
    // just used to output data, this is probably what you want.
    probes=new Monitor1dDump(.01, "probes");
    // In case the monitor should be relative to the timestep, the
    // timestep needs to be read first:
    BoutReal timestep;
    timestep = Options::root()["timestep"].withDefault(-1);
    // There is no 'slow' section in BOUT.inp, therfore it will write
    // to data/slow.dmp.0.nc
    Monitor1dDump * slow=new Monitor1dDump(timestep*2,"slow");
    // now we can add the monitors
    solver->addMonitor(probes);
    solver->addMonitor(slow);
    // the derived monitor Monitor1dData can dump 1d data, which we can
    // now add:
    probes->add(n(mesh->xstart,mesh->ystart,0),"n_up");
    probes->add(T(mesh->xstart,mesh->ystart,0),"T_up");
    // add the corner value
    slow->add(T(0,0,0),"T"); // T is already present in BOUT.dmp - but
                             // as it is a differnt file, this doesn't
                             // cause issues.
    return 0;
  }

  int rhs(BoutReal t) {
    ddt(n) = -T;
    ddt(T) = -n;
    return 0;
  }
private:
  Field3D n,T;

  Monitor1dDump * probes, * slow;
};

BOUTMAIN(MonitorExample);

