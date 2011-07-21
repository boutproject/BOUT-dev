#include <globals.hxx>
#include <optionsreader.hxx>
#include <boutexception.hxx>

// Interface for option file parsers
#include "options/optionparser.hxx"

// Individual parsers
#include "options/options_ini.hxx"

OptionsReader* OptionsReader::instance = NULL;

OptionsReader* OptionsReader::getInstance() {
  if (instance == NULL) instance = new OptionsReader(); // Create the singleton object

  return instance;
}

void OptionsReader::read(Options *options, const char *file, ...) {
  if(file == (const char*) NULL) throw new BoutException("OptionsReader::read passed NULL filename\n");

  va_list ap;  // List of arguments
  char filename[512];

  va_start(ap, file);
  vsprintf(filename, file, ap);
  va_end(ap);

  output.write("Reading options file %s\n", filename);

  // Need to decide what file format to use
  OptionParser *parser = new OptionINI();

  parser->read(options, filename);

  delete parser;
}

void OptionsReader::parseCommandLine(Options *options, int argc, char **argv) {
  // A key/value pair, separated by a '=' or a switch
  // and sections separated with an '_' but don't start with a '-'

  string buffer,buffer_peek,buffer_peek_peek;

  // Go through command-line arguments
  for (size_t i=1;i<argc;i++) {

    // Reset the section
    options = options->getRoot();

    buffer = argv[i];
    size_t dashpos = buffer.find_first_of("-");

    // Test to see if the user put spaces around the '=' sign
    if (i < argc-2) {
      buffer_peek = argv[i+1];
      buffer_peek_peek = argv[i+2];
      size_t test_next_for_equal_sign = buffer_peek.find_first_of("=");
      size_t test_next_next_for_equal_sign = buffer_peek_peek.find_first_of("=");

      // If the user did, then concatenate all the strings
      if (test_next_for_equal_sign != string::npos) {
        buffer.append(buffer_peek);
        buffer.append(buffer_peek_peek);
        i = i+2;
      } else if (dashpos != string::npos && test_next_next_for_equal_sign == string::npos) {
        buffer.append("=");
        buffer.append(buffer_peek);
        i = i+1;
      }
    } else if (i < argc-1 && dashpos != string::npos) {
      buffer_peek = argv[i+1];
      size_t test_next_for_dash = buffer_peek.find_first_of("-");

      if (test_next_for_dash == string::npos) {
        buffer.append("=");
        buffer.append(buffer_peek);
        i = i+1;
      }
    }

    if (dashpos != string::npos) buffer = buffer.substr(dashpos+1);

    size_t startpos = buffer.find_first_of("=");

    if (startpos == string::npos) {
      // Just set a flag to true
      // e.g. "restart" or "append" on command line

      options->set(buffer, true, "Command line");
    } else {
      size_t endpos = buffer.find_last_of("=");

      if(startpos != endpos) throw BoutException("\tMultiple '=' in command-line argument '%s'\n", buffer.c_str());

      string key = buffer.substr(0, startpos);
      string value = buffer.substr(startpos+1);

      size_t scorepos = key.find_first_of("_");

      if (scorepos != string::npos) {
        size_t endscore = key.find_last_of("_");

        // if (scorepos != endscore) throw BoutException("\nMultiple '_' in command-line argument '%s'\n", key.c_str());

        string section = key.substr(0,scorepos);
        key = key.substr(scorepos+1);

        options = options->getSection(section);
      }

      if(key.empty() || value.empty()) throw BoutException("\tEmpty key or value in command line '%s'\n", buffer.c_str());

      options->set(key, value, "Command line");
    }
  }
}
