#include <bout/physicsmodel.hxx>
#include <bout/sys/timer.hxx>

int main(int argc, char **argv) {
  BoutInitialise(argc, argv);

  Field3D f{0.};
  f.setBoundary("f");
  Options* opt = Options::getRoot();
  int ntests;
  OPTION(opt, ntests, 10000);


  {
    Timer timer("boundary");
    for (int i=0; i<ntests; i++) {
      f.applyBoundary();
    }
  }

  std::string name;
  opt->getSection("f")->get("bndry_all",name,"error");

  double time = Timer::getTime("boundary");
  printf("%-20s | %12.8f s", name.c_str(), time);
  
  return 0;
}

