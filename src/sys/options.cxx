#include <options.hxx>
#include <boutexception.hxx>
#include <utils.hxx>
#include <sstream>
#include <output.hxx>

#include <field_factory.hxx> // Used for parsing expressions

Options::~Options() {
  // Delete sub-sections
  for(const auto& it : sections) {
    delete it.second;
  }
}

Options* Options::root = NULL;

Options* Options::getRoot() {
  if(root == NULL) {
    // Create the singleton
    root = new Options();
  }
  return root;
}

void Options::cleanup() {
  if(root == NULL)
    return;
  delete root;
  root = NULL;
}

void Options::set(const string &key, const int &val, const string &source) {
  stringstream ss;
  ss << val;
  set(key, ss.str(), source);
}

void Options::set(const string &key, BoutReal val, const string &source) {
  if(val) {
    set(key, "true", source);
  }else
    set(key, "false", source);
}

void Options::set(const string &key, const bool &val, const string &source) {
  stringstream ss;
  ss << val;
  set(key, ss.str(), source);
}

void Options::set(const string &key, const string &val, const string &source) {
  OptionValue opt;
  opt.value = val;
  opt.source = source;
  opt.used = false;
  
  options[lowercase(key)] = opt;
}

bool Options::isSet(const string &key) {
  map<string, OptionValue>::iterator it(options.find(lowercase(key)));
  return it != options.end();
}

void Options::get(const string &key, int &val, const int &def, bool log) {
  map<string, OptionValue>::iterator it(options.find(lowercase(key)));
  if(it == options.end()) {
    val = def;
    if(log) {
      output << "\tOption " << sectionName << ":" << key << " = " << def << " (default)" << endl;
    }
    return;
  }
  
  // Use FieldFactory to evaluate expression
  // Parse the string, giving this Option pointer for the context
  // then generate a value at t,x,y,z = 0,0,0,0
  FieldGenerator* gen = FieldFactory::get()->parse( it->second.value, this );
  if(!gen) {
    throw BoutException("Couldn't get integer from %s:%s = '%s'", 
                        sectionName.c_str(), key.c_str(), it->second.value.c_str());
  }
  BoutReal rval = gen->generate(0,0,0,0);
  
  // Convert to int by rounding
  val = ROUND(rval);
  
  // Check that the value is close to an integer
  if(fabs(rval - static_cast<BoutReal>(val)) > 1e-3) {
    throw BoutException("Value for %s:%s = %e is not an integer",
                        sectionName.c_str(), key.c_str(), rval);
  }
  
  it->second.used = true;
  if(log) {
    output << "\tOption " << sectionName << ":" << it->first << " = " << val;
    if(!it->second.source.empty()) {
      // Specify the source of the setting
      output << " (" << it->second.source << ")";
    }
    output << endl;
  }
}

void Options::get(const string &key, BoutReal &val, BoutReal def, bool log) {
  map<string, OptionValue>::iterator it(options.find(lowercase(key)));
  if(it == options.end()) {
    val = def;
    if(log)
      output << "\tOption " << sectionName << ":" << key << " = " << def << " (default)" << endl;
    return;
  }
  
  // Use FieldFactory to evaluate expression
  // Parse the string, giving this Option pointer for the context
  // then generate a value at t,x,y,z = 0,0,0,0
  FieldGenerator* gen = FieldFactory::get()->parse( it->second.value, this );
  if(!gen) {
    throw BoutException("Couldn't get BoutReal from %s:%s = '%s'", 
                        sectionName.c_str(), key.c_str(), it->second.value.c_str());
  }
  val = gen->generate(0,0,0,0);
  
  // Mark this option as used
  it->second.used = true;
  
  if(log) {
    output << "\tOption " << sectionName << ":" << it->first << " = " << val;
    if(!it->second.source.empty()) {
      // Specify the source of the setting
      output << " (" << it->second.source << ")";
    }
    output << endl;
  }
}

void Options::get(const string &key, bool &val, const bool &def, bool log) {
  map<string, OptionValue>::iterator it(options.find(lowercase(key)));
  if(it == options.end()) {
    val = def;
    if(log) {
      if(def) {
        output << "\tOption " << sectionName << ":" << key << " = true   (default)" << endl;
      }else
        output << "\tOption " << sectionName << ":" << key << " = false  (default)" << endl;
    }
    return;
  }
  
  it->second.used = true;
  
  char c = toupper((it->second.value)[0]);
  if((c == 'Y') || (c == 'T') || (c == '1')) {
    val = true;
    if(log)
      output << "\tOption " << sectionName << ":" << it->first << " = true";
  } else if((c == 'N') || (c == 'F') || (c == '0')) {
    val = false;
    if(log)
      output << "\tOption " << sectionName << ":" << it->first << " = false";
  } else
    throw BoutException("\tOption '%s': Boolean expected. Got '%s'\n", 
                        it->first.c_str(), it->second.value.c_str());
  if(!it->second.source.empty() && log) {
    // Specify the source of the setting
    output << " (" << it->second.source << ")";
  }
  if(log)
    output << endl;
}

void Options::get(const string &key, string &val, const string &def, bool log) {
  map<string, OptionValue>::iterator it(options.find(lowercase(key)));
  if(it == options.end()) {
    val = def;
    if(log)
      output << "\tOption " << sectionName << ":" << key << " = " << def << " (default)" << endl;
    return;
  }
  
  val = it->second.value;
  it->second.used = true;
  
  if(log) {
    output << "\tOption " << sectionName << ":" << it->first << " = " << val;
    if(!it->second.source.empty()) {
      // Specify the source of the setting
      output << " (" << it->second.source << ")";
    }
    output << endl;
  }
}

Options* Options::getSection(const string &name) {
  if(name.empty()) {
    return this;
  }
  map<string, Options*>::iterator it(sections.find(lowercase(name)));
  if(it != sections.end())
    return it->second;
  
  // Doesn't exist yet, so add
  string secname = name;
  if(!sectionName.empty()) // prepend the section name
    secname = sectionName + ":" + secname;
  Options *o = new Options(this, secname);
  sections[lowercase(name)] = o;
  return o;
}

string Options::str() {
  if(parent == NULL) {
    return sectionName;
  }
  string name = parent->str();
  if(name.length() > 0)
    name += string(":");
  return  name + sectionName;
}

void Options::printUnused() {
  bool allused = true;
  // Check if any options are unused
  for(const auto& it : options) {
    if(!it.second.used) {
      allused = false;
      break;
    }
  }
  if(allused) {
    output << "All options used\n";
  }else {
    output << "Unused options:\n";
    for(const auto& it : options) {
      if(!it.second.used) {
        output << "\t" << sectionName << ":" << it.first << " = " << it.second.value;
        if(!it.second.source.empty())
          output << " (" << it.second.source << ")";
        output << endl;
      }
    }
  }
  for(const auto& it : sections) {
    it.second->printUnused();
  }
}
