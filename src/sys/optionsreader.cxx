#include <optionsreader.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>
#include <bout/assert.hxx>
#include <utils.hxx>

// Interface for option file parsers
#include "options/optionparser.hxx"

// Individual parsers
#include "options/options_ini.hxx"

#include <output.hxx>

OptionsReader *OptionsReader::instance = nullptr;

OptionsReader* OptionsReader::getInstance() {
  if (instance == nullptr)
    instance = new OptionsReader(); // Create the singleton object

  return instance;
}

void OptionsReader::read(Options *options, const std::string& filename) {
  TRACE("OptionsReader::read");
  if (filename.empty()) {
    throw BoutException("OptionsReader::read passed empty filename\n");
  }

  output_info << "Reading options file " << filename << "\n";

  OptionINI{}.read(options, filename);
}

void OptionsReader::write(Options *options, const std::string& filename) {
  TRACE("OptionsReader::write");
  if (filename.empty()) {
    throw BoutException("OptionsReader::write passed empty filename\n");
  }

  output_info.write(_("Writing options to file {:s}\n"), filename);

  OptionINI{}.write(options, filename);
}

void OptionsReader::parseCommandLine(Options* options, int argc, char** argv) {
  return parseCommandLine(options, std::vector<std::string>(argv, argv + argc));
}

void OptionsReader::parseCommandLine(Options *options, const std::vector<std::string>& argv) {

  // A key/value pair, separated by a '=' or a switch
  // and sections separated with an '_' but don't start with a '-'

  std::string buffer;

  const auto argc = argv.size();

  // Go through command-line arguments
  for (auto i = std::vector<std::string>::size_type{1}; i < argc; i++) {

    // Reset the section
    options = Options::getRoot();

    buffer = argv[i];
    if (buffer.length() == 0) {
      continue;
    }
    // Test if name starts with a '-', and remove if found
    if (buffer[0] == '-') {
      buffer = buffer.substr(1); // Remove the first character (-)
      if (buffer.length() == 0) {
        throw BoutException(_("Invalid command line option '-' found - maybe check whitespace?"));
      }
    }
    // Test to see if the user put spaces around the '=' sign
    if (i < argc - 1) {
      if (buffer[buffer.length() - 1] == '=') {
        // Space after '=' sign

        i++;
        buffer.append(argv[i]);
      } else if (argv[i + 1][0] == '=') {
        // Space before '=' sign

        i++;
        buffer.append(argv[i]);

        if ((argv[i][1] == 0) && (i < argc - 1)) {
          // End of string, so space after '=' sign too
          i++;
          buffer.append(argv[i]);
        }
      }
    }

    size_t startpos = buffer.find_first_of('=');

    if (startpos == std::string::npos) {
      // Just set a flag to true
      // e.g. "restart" or "append" on command line

      options->set(buffer, true, "Command line");
    } else {
      size_t endpos = buffer.find_last_of('=');

      if (startpos != endpos) {
        throw BoutException(_("\tMultiple '=' in command-line argument '{:s}'\n"),
                            buffer);
      }

      std::string key = trim(buffer.substr(0, startpos));
      std::string value = trim(buffer.substr(startpos+1));
      
      size_t scorepos;
      while((scorepos = key.find_first_of(':')) != std::string::npos) {
	// sub-section
	std::string section = key.substr(0,scorepos);
	key = trim(key.substr(scorepos+1));
	options = options->getSection(section);
      }

      if (key.empty() || value.empty()) {
        throw BoutException(_("\tEmpty key or value in command line '{:s}'\n"), buffer);
      }

      options->set(key, value, _("Command line"));
    }
  }
}
