#include <options.hxx>
#include <boutexception.hxx>
#include <utils.hxx>
#include <output.hxx>
#include <field_factory.hxx> // Used for parsing expressions

#include <iomanip>
#include <sstream>

/// The source label given to default values
const string DEFAULT_SOURCE{"default"};

Options::Options() : parent_instance(nullptr) {}

Options *Options::root_instance = nullptr;

Options &Options::root() {
  if (root_instance == nullptr) {
    // Create the singleton
    root_instance = new Options();
  }
  return *root_instance;
}

void Options::cleanup() {
  if (root_instance == nullptr)
    return;
  delete root_instance;
  root_instance = nullptr;
}

Options& Options::operator[](const string &name) {
  // Mark this object as being a section
  is_section = true;

  if (name.empty()) {
    return *this;
  }

  // Find and return if already exists
  auto it = children.find(lowercase(name));
  if (it != children.end()) {
    return it->second;
  }

  // Doesn't exist yet, so add
  string secname = name;
  if (!full_name.empty()) { // prepend the section name
    secname = full_name + ":" + secname;
  }

  // emplace returns a pair with iterator first, boolean (insert yes/no) second
  auto pair_it = children.emplace(lowercase(name), Options{this, secname});
  
  return pair_it.first->second;
}

const Options &Options::operator[](const string &name) const {
  if (!is_section) {
    throw BoutException("Option %s is not a section", full_name.c_str());
  }

  if (name.empty()) {
    return *this;
  }

  // Find and return if already exists
  auto it = children.find(lowercase(name));
  if (it == children.end()) {
    // Doesn't exist
    throw BoutException("Option %s:%s does not exist", full_name.c_str(), name.c_str());
  }

  return it->second;
}

template<>
void Options::assign<bool>(bool val, const string &source) {
  if (val) {
    _set("true", source, false);
  } else {
    _set("false", source, false);
  }
}

template<>
void Options::assign<BoutReal>(BoutReal val, const string &source) {
  std::stringstream ss;
  // Make sure the precision is large enough to hold a BoutReal
  ss << std::scientific << std::setprecision(17) << val;
  _set(ss.str(), source, false);
}

void Options::_set(string val, string source, bool force) {
  if (isSet()) {
    // Check if current value the same as new value
    if (value.value != val) {
      if (force or value.source != source) {
        output_warn << "\tOption " << full_name << " = " << value.value
                    << " (" << value.source << ") overwritten with:"
                    << "\n"
                    << "\t\t" << full_name << " = " << val
                    << " (" << source << ")\n";
      } else {
        throw BoutException("Options: Setting a value from same source (%s) to new value "
                            "'%s' - old value was '%s'.",
                            source.c_str(), val.c_str(), value.value.c_str());
      }
    }
  }
  
  value.value = std::move(val);
  value.source = std::move(source);
  value.used = false;
  is_value = true;
}

bool Options::isSet() const {
  // Check if no value
  if (!is_value) {
    return false;
  }

  // Ignore if set from default
  if (value.source == DEFAULT_SOURCE) {
    return false;
  }
  
  return true;
}

template<>
int Options::withDefault<int>(int def) {  
  if (!is_value) {
    // Option not found
    // Set the option, with source "default". This is to ensure that:
    //   a) the same option has a consistent default value
    //   b) the value used can be recorded in the output settings file
    assign(def, DEFAULT_SOURCE);
    value.used = true; // Mark the option as used
    
    output_info << "\tOption " << full_name << " = " << def
                << " (default)" << endl;
    return def;
  }
  
  // Use FieldFactory to evaluate expression
  // Parse the string, giving this Option pointer for the context
  // then generate a value at t,x,y,z = 0,0,0,0
  auto gen = FieldFactory::get()->parse( value.value, this );
  if (!gen) {
    throw BoutException("Couldn't get integer from %s = '%s'", full_name.c_str(), value.value.c_str());
  }
  BoutReal rval = gen->generate(0, 0, 0, 0);

  // Convert to int by rounding
  int val = ROUND(rval);

  // Check that the value is close to an integer
  if (fabs(rval - static_cast<BoutReal>(val)) > 1e-3) {
    throw BoutException("Value for %s = %e is not an integer", full_name.c_str(), rval);
  }
  
  // Check if this was previously set as a default option
  if (value.source == DEFAULT_SOURCE) {
    // Check that the default values are the same
    if (def != val) {
      throw BoutException("Inconsistent default values for '%s': '%d' then '%d'",
                          full_name.c_str(), val, def);
    }
  }

  value.used = true;

  output_info << "\tOption " << full_name << " = " << val;
  if (!value.source.empty()) {
    // Specify the source of the setting
    output_info << " (" << value.source << ")";
  }
  output_info << endl;

  return val;
}

template<>
BoutReal Options::withDefault<BoutReal>(BoutReal def) {
  if (!is_value) {
    assign(def, DEFAULT_SOURCE);
    value.used = true; // Mark the option as used
    
    output_info << "\tOption " << full_name << " = " << def
                << " (" << DEFAULT_SOURCE << ")" << endl;
    return def;
  }
  
  // Use FieldFactory to evaluate expression
  // Parse the string, giving this Option pointer for the context
  // then generate a value at t,x,y,z = 0,0,0,0
  std::shared_ptr<FieldGenerator>  gen = FieldFactory::get()->parse( value.value, this );
  if (!gen) {
    throw BoutException("Couldn't get BoutReal from %s = '%s'", full_name.c_str(), value.value.c_str());
  }
  BoutReal val = gen->generate(0, 0, 0, 0);

  // Check if this was previously set as a default option
  if (value.source == DEFAULT_SOURCE) {
    // Check that the default values are the same
    if (fabs(def - val) > 1e-10) {
      throw BoutException("Inconsistent default values for '%s': '%e' then '%e'",
                          full_name.c_str(), val, def);
    }
  }

  // Mark this option as used
  value.used = true;

  output_info << "\tOption " << full_name << " = " << val;
  if (!value.source.empty()) {
    // Specify the source of the setting
    output_info << " (" << value.source << ")";
  }
  output_info << endl;

  return val;
}

template<> 
bool Options::withDefault<bool>(bool def) {
  if (!is_value) {
    assign(def, DEFAULT_SOURCE);
    value.used = true; // Mark the option as used
    
    if (def) {
      output_info << "\tOption " << full_name << " = true";
    } else {
      output_info << "\tOption " << full_name << " = false";
    }
    output_info << "   (" << DEFAULT_SOURCE << ")" << endl;
    return def;
  }
  
  value.used = true;

  bool val;
  char c = static_cast<char>(toupper((value.value)[0]));
  if ((c == 'Y') || (c == 'T') || (c == '1')) {
    val = true;
    output_info << "\tOption " << full_name << " = true";
  } else if ((c == 'N') || (c == 'F') || (c == '0')) {
    val = false;
    output_info << "\tOption " << full_name << " = false";
  } else {
    throw BoutException("\tOption '%s': Boolean expected. Got '%s'\n", full_name.c_str(),
                        value.value.c_str());
  }
  if (!value.source.empty()) {
    // Specify the source of the setting
    output_info << " (" << value.source << ")";
  }
  output_info << endl;

  return val;
}

template<> 
std::string Options::withDefault<std::string>(std::string def) {
  if (!is_value) {
    _set(def, DEFAULT_SOURCE, false);
    value.used = true; // Mark the option as used
    
    output_info << "\tOption " << full_name << " = " << def
                << " (" << DEFAULT_SOURCE << ")" << endl;
    return def;
  }
  
  // Check if this was previously set as a default option
  if (value.source == DEFAULT_SOURCE) {
    // Check that the default values are the same
    if (value.value != def) {
      throw BoutException("Inconsistent default values for '%s': '%s' then '%s'",
                          full_name.c_str(), value.value.c_str(), def.c_str());
    }
  }

  value.used = true;

  output_info << "\tOption " << full_name << " = " << value.value;
  if (!value.source.empty()) {
    // Specify the source of the setting
    output_info << " (" << value.source << ")";
  }
  output_info << endl;

  return value.value;
}

void Options::printUnused() const {
  bool allused = true;
  // Check if any options are unused
  for (const auto &it : children) {
    if ( it.second.is_value &&
         !it.second.value.used) {
      allused = false;
      break;
    }
  }
  if (allused) {
    output_info << "All options used\n";
  } else {
    output_info << "Unused options:\n";
    for (const auto &it : children) {
      if (it.second.is_value &&
          !it.second.value.used) {
        output_info << "\t" << full_name << ":" << it.first << " = " << it.second.value.value;
        if (!it.second.value.source.empty())
          output_info << " (" << it.second.value.source << ")";
        output_info << endl;
      }
    }
  }
  for (const auto &it : children) {
    if (it.second.is_section) {
      it.second.printUnused();
    }
  }
}

void Options::cleanCache() {
  FieldFactory::get()->cleanCache();
}

std::map<string, Options::OptionValue> Options::values() const {
  std::map<string, OptionValue> options;
  for (const auto &it : children) {
    if (it.second.is_value) {
      options[it.first] = it.second.value;
    }
  }
  return options;
}

std::map<string, const Options*> Options::subsections() const {
  std::map<string, const Options*> sections;
  for (const auto &it : children) {
    if (it.second.is_section) {
      sections[it.first] = &it.second;
    }
  }
  return sections;
}
