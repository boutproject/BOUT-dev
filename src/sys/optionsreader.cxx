#include <optionsreader.hxx>
#include <boutexception.hxx>
#include <utils.hxx>

// Interface for option file parsers
#include "options/optionparser.hxx"

// Individual parsers
#include "options/options_ini.hxx"

#include <output.hxx>

OptionsReader* OptionsReader::instance = NULL;

OptionsReader* OptionsReader::getInstance() {
  if (instance == NULL) instance = new OptionsReader(); // Create the singleton object

  return instance;
}

void OptionsReader::read(Options *options, const char *file, ...) {
  if(file == (const char*) NULL) throw new BoutException("OptionsReader::read passed NULL filename\n");

  int buf_len=512;
  char * filename=new char[buf_len];

  myvsnprintf(filename,buf_len, file, ap);

  output.write("Reading options file %s\n", filename);

  // Need to decide what file format to use
  OptionParser *parser = new OptionINI();

  parser->read(options, filename);

  delete[] filename;
  delete parser;
}

void OptionsReader::parseCommandLine(Options *options, int argc, char **argv) {
  // A key/value pair, separated by a '=' or a switch
  // and sections separated with an '_' but don't start with a '-'

  string buffer;

  // Go through command-line arguments
  for (size_t i=1;i<argc;i++) {

    // Reset the section
    options = options->getRoot();

    buffer = argv[i];
    // Test if name starts with a '-', and remove if found
    if (buffer[0] == '-')
      buffer = buffer.substr(1);  // Remove the first character (-)
    
    // Test to see if the user put spaces around the '=' sign
    if (i < argc-1) {
      if(buffer[buffer.length()-1] == '=') {
        // Space after '=' sign
        
        i++;
        buffer.append(argv[i]);
        
      }else if(argv[i+1][0] == '=') {
        // Space before '=' sign
        
        i++;
        buffer.append(argv[i]);
        
        if((argv[i][1] == 0) && (i < argc-2)) {
          // End of string, so space after '=' sign too
          
          i++;
          buffer.append(argv[i]);
        }
      }
    }
    
    size_t startpos = buffer.find_first_of("=");

    if (startpos == string::npos) {
      // Just set a flag to true
      // e.g. "restart" or "append" on command line

      options->set(buffer, true, "Command line");
    } else {
      size_t endpos = buffer.find_last_of("=");

      if(startpos != endpos) throw BoutException("\tMultiple '=' in command-line argument '%s'\n", buffer.c_str());

      string key = trim(buffer.substr(0, startpos));
      string value = trim(buffer.substr(startpos+1));
      
      size_t scorepos;
      while((scorepos = key.find_first_of(":")) != string::npos) {
	// sub-section
	string section = key.substr(0,scorepos);
	key = trim(key.substr(scorepos+1));
	options = options->getSection(section);
      }
      
      if(key.empty() || value.empty()) throw BoutException("\tEmpty key or value in command line '%s'\n", buffer.c_str());

      options->set(key, value, "Command line");
    }
  }
}
