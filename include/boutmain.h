/*!************************************************************************
 * Main function
 **************************************************************************/

/*
 * This is a sloppy way to include main but for now it is a stop gap for converting all the examples
 */

int main(int argc, char **argv)
{
  if(bout_init(argc, argv)) {
    fprintf(stderr, "ERROR INITIALISING BOUT++. ABORTING\n");
    return 1;
  }

  bout_run();

  bout_finish();

  return(0);
}
