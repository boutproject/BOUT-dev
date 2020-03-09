/*
 * Laplacian inversion
 *
 */

#include <bout.hxx>
#include <invert_laplace.hxx>
#include <field_factory.hxx>
#include <cstring>

int main(int argc, char **argv) {

  // Initialise BOUT++, setting up mesh
  BoutInitialise(argc, argv);

  FieldFactory f(mesh);

// String containing "not diagonally dominant" error
  std::string ndd_message("not diagonally dominant");
  int local_error=0;

  Field3D input = f.create3D("(1-gauss(x-0.5,0.2))*gauss(y-pi)*gauss(z-pi)");
  Field2D a = f.create2D("gauss(x) * sin(y)");
  Field2D c = f.create2D("sin(x) * gauss(x-0.5) * gauss(y-pi)");
  Field2D d = f.create2D("y - pi/2");
  SAVE_ONCE4(input, a, c, d);

#define RUN_TEST(test,error,set_command,test_command) {				\
  local_error = 0;								\
  error = 0;									\
  Laplacian *lap = Laplacian::create();						\
  lap->resetSolver();								\
  set_command;									\
  try {										\
    test_command;								\
  } catch (BoutException &e) {							\
    std::string message(e.what());						\
    if (message.find(ndd_message) != std::string::npos){			\
      local_error = 1;								\
    }										\
    test = 0.0; 								\
  }										\
  SAVE_ONCE(test);								\
  MPI_Reduce(&local_error,&error,1,MPI_INTEGER,MPI_MAX, 0, MPI_COMM_WORLD);	\
  SAVE_ONCE(error);								\
  delete lap;									\
}

#define RUN_TEST_OLD_IFACE(test,error,test_command) {				\
  local_error = 0;								\
  error = 0;									\
  try {										\
    test_command;								\
  } catch (BoutException &e) {							\
    std::string message(e.what());						\
    if (message.find(ndd_message) != std::string::npos){			\
      local_error = 1;								\
    }										\
    test = 0.0; 								\
  }										\
  SAVE_ONCE(test);								\
  MPI_Reduce(&local_error,&error,1,MPI_INTEGER,MPI_MAX, 0, MPI_COMM_WORLD);	\
  SAVE_ONCE(error);								\
}

  Field3D flag0;
  int flag0_error;
  RUN_TEST(flag0,flag0_error,
    lap->setFlags(0),
    flag0 = lap->solve(input)
  );

  Field3D flag3;
  int flag3_error;
  RUN_TEST(flag3,flag3_error,
    lap->setFlags(3),
    flag3 = lap->solve(input)
  );
  
  Field3D flag0a;
  int flag0a_error;
  RUN_TEST(flag0a,flag0a_error,
    lap->setCoefA(a); lap->setFlags(0);,
    flag0a = lap->solve(input)
  );

  Field3D flag3a;
  int flag3a_error;
  RUN_TEST(flag3a,flag3a_error,
    lap->setCoefA(a); lap->setFlags(3);,
    flag3a = lap->solve(input)
  );

  Field3D res0a  = Delp2(flag0a);
  SAVE_ONCE(res0a);

  Field3D flag0ac;
  int flag0ac_error;
  RUN_TEST_OLD_IFACE(flag0ac,flag0ac_error,
    flag0ac = invert_laplace(input, 0, &a, &c);
  );

  Field3D flag3ac;
  int flag3ac_error;
  RUN_TEST_OLD_IFACE(flag3ac,flag3ac_error,
    flag3ac = invert_laplace(input, 3, &a, &c);
  );

  Field3D flag0ad;
  int flag0ad_error;
  RUN_TEST_OLD_IFACE(flag0ad,flag0ad_error,
    flag0ad = invert_laplace(input, 0, &a, nullptr, &d);
  );

  Field3D flag3ad;
  int flag3ad_error;
  RUN_TEST_OLD_IFACE(flag3ad,flag3ad_error,
    flag3ad = invert_laplace(input, 3, &a, nullptr, &d);
  );

  /// Test new interface and INVERT_IN/OUT_SET flags

  Field2D set_to = f.create2D("cos(2*y)*(x - 0.5)");
  SAVE_ONCE(set_to);

  Field3D flagis;
  int flagis_error;
  RUN_TEST(flagis,flagis_error,
    lap->setFlags(4096);,
    flagis = lap->solve(input, set_to);
  );

  Field3D flagos;
  int flagos_error;
  RUN_TEST(flagos,flagos_error,
    lap->setFlags(8192);,
    flagos = lap->solve(input, set_to);
  );

  Field3D flagisa;
  int flagisa_error;
  RUN_TEST(flagisa,flagisa_error,
    lap->setCoefA(a); lap->setFlags(4096);,
    flagisa = lap->solve(input, set_to);
  );

  Field3D flagosa;
  int flagosa_error;
  RUN_TEST(flagosa,flagosa_error,
    lap->setCoefA(a); lap->setFlags(8192);,
    flagosa = lap->solve(input, set_to);
  );

  Field3D flagisac;
  int flagisac_error;
  RUN_TEST(flagisac,flagisac_error,
    lap->setCoefA(a); lap->setFlags(4096); lap->setCoefC(c),
    flagisac = lap->solve(input, set_to);
  );

  Field3D flagosac;
  int flagosac_error;
  RUN_TEST(flagosac,flagosac_error,
    lap->setCoefA(a); lap->setFlags(8192); lap->setCoefC(c),
    flagosac = lap->solve(input, set_to);
  );

  Field3D flagisad;
  int flagisad_error;
  RUN_TEST(flagisad,flagisad_error,
    lap->setCoefA(a); lap->setFlags(4096); lap->setCoefC(1.0); lap->setCoefD(d),
    flagisad = lap->solve(input, set_to);
  );

  Field3D flagosad;
  int flagosad_error;
  RUN_TEST(flagosad,flagosad_error,
    lap->setCoefA(a); lap->setFlags(8192); lap->setCoefC(1.0); lap->setCoefD(d),
    flagosad = lap->solve(input, set_to);
  );

  // Write and close the output file

  dump.write();
  dump.close();

  output << "\nFinished running test. Triggering error to quit\n\n";

  MPI_Barrier(BoutComm::get()); // Wait for all processors to write data

  BoutFinalise();
  return 0;
}
