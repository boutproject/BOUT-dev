
#include "bout.hxx"
#include "bout/sys/gettext.hxx"

int main(int argc, char**argv) {
  // Setting the i18n environment
  if (setlocale (LC_ALL, "") == NULL) {
    throw BoutException("Failed to set locale");
  }
  bindtextdomain ("test", "/home/bd512/BOUT-dev/locale/");
  textdomain ("test");

  BoutInitialise(argc, argv);
  
  output.write(_("Hello World\n"));

  BoutFinalise();
  return 0;
}
