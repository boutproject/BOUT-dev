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
  
  string buffer;
  
  // Go through command-line arguments
  for(size_t i=1;i<argc;i++) {
    buffer = argv[i];
    size_t startpos = buffer.find_first_of("=");
    
    if(startpos == string::npos) {
      // Just set a flag to true
      // e.g. "restart" or "append" on command line
      
      options->set(buffer, true, "Command line");
    }else {
      size_t endpos = buffer.find_last_of("=");
      
      if(startpos != endpos) {
        throw BoutException("\tMultiple '=' in command-line argument '%s'\n", buffer.c_str());
      }
      
      string key = buffer.substr(0, startpos);
      string value = buffer.substr(startpos+1);
      if(key.empty() || value.empty()) {
        throw BoutException("\tEmpty key or value in command line '%s'\n", buffer.c_str());
      }
      options->set(key, value, "Command line");
    }
  }
}
