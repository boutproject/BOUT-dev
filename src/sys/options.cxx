#include <options.hxx>
#include <boutexception.hxx>
#include <utils.hxx>
#include <output.hxx>
#include <field_factory.hxx> // Used for parsing expressions

#include <iomanip>
#include <sstream>

/// The source label given to default values
const string DEFAULT_SOURCE{"default"};

Options::Options() : parent(nullptr) {}

Options::~Options() {
  // Delete sub-sections
  for (const auto &it : sections) {
    delete it.second;
  }
}

Options *Options::root = nullptr;

Options *Options::getRoot() {
  if (root == nullptr) {
    // Create the singleton
    root = new Options();
  }
  return root;
}

void Options::cleanup() {
  if (root == nullptr)
    return;
  delete root;
  root = nullptr;
}

void Options::set(const string &key, const int &val, const string &source, bool force) {
  std::stringstream ss;
  ss << val;
  _set(key, ss.str(), source, force);
}

void Options::set(const string &key, const bool &val, const string &source, bool force) {
  if (val) {
    _set(key, "true", source, force);
  } else {
    _set(key, "false", source, force);
  }
}

void Options::set(const string &key, BoutReal val, const string &source, bool force) {
  std::stringstream ss;
  // Make sure the precision is large enough to hold a BoutReal
  ss << std::scientific << std::setprecision(17) << val;
  _set(key, ss.str(), source, force);
}

void Options::_set(const string &key, const string &val, const string &source,
                   bool force) {
  OptionValue opt;
  opt.value = val;
  opt.source = source;
  opt.used = false;
  string key_l = lowercase(key);
  if (this->isSet(key_l)) {
    const OptionValue &oldopt = options[key_l];
    if (oldopt.value != opt.value) {
      if (force or oldopt.source != opt.source) {
        output_warn << "\tOption " << sectionName << ":" << key << " = " << oldopt.value
                    << " (" << oldopt.source << ") overwritten with:"
                    << "\n"
                    << "\t\t" << sectionName << ":" << key << " = " << opt.value
                    << " (" << opt.source << ")\n";
      } else {
        throw BoutException("Options: Setting a value from same source (%s) to new value "
                            "'%s' - old value was '%s'.",
                            source.c_str(), opt.value.c_str(), oldopt.value.c_str());
      }
    }
  }
  options[key_l] = opt;
}

bool Options::isSet(const string &key) {
  std::map<string, OptionValue>::iterator it(options.find(lowercase(key)));
  return it != options.end();
}

void Options::get(const string &key, int &val, int def) {
  std::map<string, OptionValue>::iterator it(options.find(lowercase(key)));
  if (it == options.end()) {
    // Option not found
    // Set the option, with source "default". This is to ensure that:
    //   a) the same option has a consistent default value
    //   b) the value used can be recorded in the output settings file
    set(key, def, DEFAULT_SOURCE);
    options[lowercase(key)].used = true; // Mark the option as used

    // Set the return value
    val = def;
    output_info << "\tOption " << sectionName << ":" << key << " = " << def
                << " (default)" << endl;
    return;
  }

  // Use FieldFactory to evaluate expression
  // Parse the string, giving this Option pointer for the context
  // then generate a value at t,x,y,z = 0,0,0,0
  std::shared_ptr<FieldGenerator>  gen = FieldFactory::get()->parse( it->second.value, this );
  if (!gen) {
    throw BoutException("Couldn't get integer from %s:%s = '%s'", sectionName.c_str(),
                        key.c_str(), it->second.value.c_str());
  }
  BoutReal rval = gen->generate(0, 0, 0, 0);

  // Convert to int by rounding
  val = ROUND(rval);

  // Check that the value is close to an integer
  if (fabs(rval - static_cast<BoutReal>(val)) > 1e-3) {
    throw BoutException("Value for %s:%s = %e is not an integer", sectionName.c_str(),
                        key.c_str(), rval);
  }

  // Check if this was previously set as a default option
  if (it->second.source == DEFAULT_SOURCE) {
    // Check that the default values are the same
    if (def != val) {
      throw BoutException("Inconsistent default values for '%s': '%d' then '%d'",
                          key.c_str(), val, def);
    }
  }

  it->second.used = true;

  output_info << "\tOption " << sectionName << ":" << it->first << " = " << val;
  if (!it->second.source.empty()) {
    // Specify the source of the setting
    output_info << " (" << it->second.source << ")";
  }
  output_info << endl;
}

void Options::get(const string &key, BoutReal &val, BoutReal def) {
  std::map<string, OptionValue>::iterator it(options.find(lowercase(key)));
  if (it == options.end()) {
    set(key, def, "default");
    options[lowercase(key)].used = true; // Mark the option as used

    val = def;

    output_info << "\tOption " << sectionName << ":" << key << " = " << def
                << " (default)" << endl;
    return;
  }

  // Use FieldFactory to evaluate expression
  // Parse the string, giving this Option pointer for the context
  // then generate a value at t,x,y,z = 0,0,0,0
  std::shared_ptr<FieldGenerator>  gen = FieldFactory::get()->parse( it->second.value, this );
  if (!gen) {
    throw BoutException("Couldn't get BoutReal from %s:%s = '%s'", sectionName.c_str(),
                        key.c_str(), it->second.value.c_str());
  }
  val = gen->generate(0, 0, 0, 0);

  // Check if this was previously set as a default option
  if (it->second.source == DEFAULT_SOURCE) {
    // Check that the default values are the same
    if (fabs(def - val) > 1e-10) {
      throw BoutException("Inconsistent default values for '%s': '%e' then '%e'",
                          key.c_str(), val, def);
    }
  }

  // Mark this option as used
  it->second.used = true;

  output_info << "\tOption " << sectionName << ":" << it->first << " = " << val;
  if (!it->second.source.empty()) {
    // Specify the source of the setting
    output_info << " (" << it->second.source << ")";
  }
  output_info << endl;
}

void Options::get(const string &key, bool &val, bool def) {
  std::map<string, OptionValue>::iterator it(options.find(lowercase(key)));
  if (it == options.end()) {
    set(key, def, "default");
    options[lowercase(key)].used = true; // Mark the option as used

    val = def;

    if (def) {
      output_info << "\tOption " << sectionName << ":" << key << " = true   (default)"
                  << endl;
    } else {
      output_info << "\tOption " << sectionName << ":" << key << " = false  (default)"
                  << endl;
    }
    return;
  }

  it->second.used = true;
  
  char c = static_cast<char>(toupper((it->second.value)[0]));
  if ((c == 'Y') || (c == 'T') || (c == '1')) {
    val = true;
    output_info << "\tOption " << sectionName << ":" << it->first << " = true";
  } else if ((c == 'N') || (c == 'F') || (c == '0')) {
    val = false;
    output_info << "\tOption " << sectionName << ":" << it->first << " = false";
  } else {
    throw BoutException("\tOption '%s': Boolean expected. Got '%s'\n", it->first.c_str(),
                        it->second.value.c_str());
  }
  if (!it->second.source.empty()) {
    // Specify the source of the setting
    output_info << " (" << it->second.source << ")";
  }
  output_info << endl;
}

void Options::get(const string &key, string &val, const string &def) {
  std::map<string, OptionValue>::iterator it(options.find(lowercase(key)));
  if (it == options.end()) {
    set(key, def, "default");
    options[lowercase(key)].used = true; // Mark the option as used

    val = def;

    output_info << "\tOption " << sectionName << ":" << key << " = " << def
                << " (default)" << endl;
    return;
  }

  val = it->second.value;

  // Check if this was previously set as a default option
  if (it->second.source == DEFAULT_SOURCE) {
    // Check that the default values are the same
    if (val != def) {
      throw BoutException("Inconsistent default values for '%s': '%s' then '%s'",
                          key.c_str(), val.c_str(), def.c_str());
    }
  }

  it->second.used = true;

  output_info << "\tOption " << sectionName << ":" << it->first << " = " << val;
  if (!it->second.source.empty()) {
    // Specify the source of the setting
    output_info << " (" << it->second.source << ")";
  }
  output_info << endl;
}

Options *Options::getSection(const string &name) {
  if (name.empty()) {
    return this;
  }
  std::map<string, Options *>::iterator it(sections.find(lowercase(name)));
  if (it != sections.end()) {
    return it->second;
  }

  // Doesn't exist yet, so add
  string secname = name;
  if (!sectionName.empty()) { // prepend the section name
    secname = sectionName + ":" + secname;
  }
  Options *o = new Options(this, secname);
  sections[lowercase(name)] = o;
  return o;
}

string Options::str() {
  // Note: name of parent already prepended in getSection
  return sectionName;
}

void Options::printUnused() {
  bool allused = true;
  // Check if any options are unused
  for (const auto &it : options) {
    if (!it.second.used) {
      allused = false;
      break;
    }
  }
  if (allused) {
    output_info << "All options used\n";
  } else {
    output_info << "Unused options:\n";
    for (const auto &it : options) {
      if (!it.second.used) {
        output_info << "\t" << sectionName << ":" << it.first << " = " << it.second.value;
        if (!it.second.source.empty())
          output_info << " (" << it.second.source << ")";
        output_info << endl;
      }
    }
  }
  for (const auto &it : sections) {
    it.second->printUnused();
  }
}

void Options::cleanCache() {
  FieldFactory::get()->cleanCache();
}
