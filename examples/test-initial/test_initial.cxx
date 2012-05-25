/*
 * Initial profiles regression test
 * 
 * Check that initial profiles are valid, and do
 * not depend on number of processors
 *
 */

#include <bout.hxx>
#include <boutmain.hxx>
#include <initialprofiles.hxx>

int physics_init(bool restarting)
{
  int status;
  
  // Test option 0
  Field2D var0_2D;
  Field3D var0_3D;
  status = initial_profile("var0", var0_2D);
  if(finite(var0_2D) && (!status)) {
    output << "Option 0 2D initialisation: SUCCESS\n";
  }else
    output << "Option 0 2D initialisation: FAILED\n";
  status = initial_profile("var0", var0_3D);
  if(finite(var0_3D) && (!status)) {
    output << "Option 0 3D initialisation: SUCCESS\n";
  }else
    output << "Option 0 3D initialisation: FAILED\n";
  dump.add(var0_2D, "var0_2D");
  dump.add(var0_3D, "var0_3D");
  
  // Test option 1
  Field2D var1_2D;
  Field3D var1_3D;
  status = initial_profile("var1", var1_2D);
  if(finite(var1_2D) && (!status)) {
    output << "Option 1 2D initialisation: SUCCESS\n";
  }else
    output << "Option 1 2D initialisation: FAILED\n";
  status = initial_profile("var1", var1_3D);
  if(finite(var1_3D) && (!status)) {
    output << "Option 1 3D initialisation: SUCCESS\n";
  }else
    output << "Option 1 3D initialisation: FAILED\n";
  dump.add(var1_2D, "var1_2D");
  dump.add(var1_3D, "var1_3D");

  // Test option 2
  Field2D var2_2D;
  Field3D var2_3D;
  status = initial_profile("var2", var2_2D);
  if(finite(var2_2D) && (!status)) {
    output << "Option 2 2D initialisation: SUCCESS\n";
  }else
    output << "Option 2 2D initialisation: FAILED\n";
  status = initial_profile("var1", var2_3D);
  if(finite(var2_3D) && (!status)) {
    output << "Option 2 3D initialisation: SUCCESS\n";
  }else
    output << "Option 2 3D initialisation: FAILED\n";
  dump.add(var2_2D, "var2_2D");
  dump.add(var2_3D, "var2_3D");

  // Test option 3
  Field2D var3_2D;
  Field3D var3_3D;
  status = initial_profile("var3", var3_2D);
  if(finite(var3_2D) && (!status)) {
    output << "Option 3 2D initialisation: SUCCESS\n";
  }else
    output << "Option 3 2D initialisation: FAILED\n";
  status = initial_profile("var3", var3_3D);
  if(finite(var3_3D) && (!status)) {
    output << "Option 3 3D initialisation: SUCCESS\n";
  }else
    output << "Option 3 3D initialisation: FAILED\n";
  dump.add(var3_2D, "var3_2D");
  dump.add(var3_3D, "var3_3D");

  // Test option 4
  Field2D var4_2D;
  Field3D var4_3D;
  status = initial_profile("var4", var4_2D);
  if(finite(var4_2D) && (!status)) {
    output << "Option 4 2D initialisation: SUCCESS\n";
  }else
    output << "Option 4 2D initialisation: FAILED\n";
  status = initial_profile("var4", var4_3D);
  if(finite(var4_3D) && (!status)) {
    output << "Option 4 3D initialisation: SUCCESS\n";
  }else
    output << "Option 4 3D initialisation: FAILED\n";
  dump.add(var4_2D, "var4_2D");
  dump.add(var4_3D, "var4_3D");
  
  // Write data to file
  dump.write();
  dump.close();
  
  // Need to wait for all processes to finish writing
  MPI_Barrier(BoutComm::get());

  // Send an error code so quits
  return 1;
}

int physics_run(BoutReal t)
{
  // Doesn't do anything
  return 1;
}
