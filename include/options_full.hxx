#pragma once

#include "options.hxx"
template <typename T>
void Options::_set(T val, std::string source, bool force){
  // If already set, and not time evolving then check for changing values
  // If a variable has a "time_dimension" attribute then it is assumed
  // that updates to the value is ok and don't need to be forced.
  if (isSet() && (attributes.find("time_dimension") == attributes.end())) {
    // Check if current value the same as new value
    if (!bout::utils::variantEqualTo(value, val)) {
      if (force or !bout::utils::variantEqualTo(attributes["source"], source)) {
	output_warn << _("\tOption ") << full_name << " = "
		    << bout::utils::variantToString(value) << " ("
		    << bout::utils::variantToString(attributes["source"])
		    << _(") overwritten with:") << "\n"
		    << "\t\t" << full_name << " = " << toString(val) << " (" << source
		    << ")\n";
      } else {
	throw BoutException(
	      _("Options: Setting a value from same source ({:s}) to new value "
		"'{:s}' - old value was '{:s}'."),
	      source, toString(val), bout::utils::variantToString(value));
      }
    }
  }
  value = std::move(val);
  attributes["source"] = std::move(source);
  value_used = false;
  is_value = true;
}

#ifndef OPT_FULL_EXTERN
#define OPT_FULL_EXTERN extern
#endif

OPT_FULL_EXTERN template void Options::_set<bool>(bool, std::string, bool);
OPT_FULL_EXTERN template void Options::_set<int>(int, std::string, bool);
OPT_FULL_EXTERN template void Options::_set<BoutReal>(BoutReal, std::string, bool);
OPT_FULL_EXTERN template void Options::_set<std::string>(std::string, std::string, bool);

template <typename T> T Options::withDefault(T def) {

  // Set the type
  attributes["type"] = bout::utils::typeName<T>();

  if (!is_value) {
    // Option not found
    assign(def, DEFAULT_SOURCE);
    value_used = true; // Mark the option as used

    output_info << _("\tOption ") << full_name << " = " << def << " (" << DEFAULT_SOURCE
		<< ")" << std::endl;
    return def;
  }
  T val = as<T>(def);
  // Check if this was previously set as a default option
  if (bout::utils::variantEqualTo(attributes.at("source"), DEFAULT_SOURCE)) {
    // Check that the default values are the same
    if (!similar(val, def)) {
      throw BoutException("Inconsistent default values for '{:s}': '{:s}' then '{:s}'",
			  full_name, bout::utils::variantToString(value),
			  toString(def));
    }
  }
  return val;
}

#define OWD(type) \
  OPT_FULL_EXTERN template type Options::withDefault<type>(type);	\
  OPT_FULL_EXTERN template type Options::withDefault<type>(type) const;

OWD(bool)
OWD(int)
OWD(BoutReal)
OWD(long)
OWD(long long)
OWD(uint)
OWD(std::string)
OWD(Field2D)
OWD(Field3D)

Options& Options::withDefault(const Options& def) {
  // if def is a section, then it does not make sense to try to use it as a default for
  // a value
  ASSERT0(def.is_value);

  if (!is_value) {
    // Option not found
    *this = def;

    output_info << _("\tOption ") << full_name << " = " << def.full_name << " ("
		<< DEFAULT_SOURCE << ")" << std::endl;
  } else {
    // Check if this was previously set as a default option
    if (bout::utils::variantEqualTo(attributes.at("source"), DEFAULT_SOURCE)) {
      // Check that the default values are the same
      if (!similar(bout::utils::variantToString(value),
		   bout::utils::variantToString(def.value))) {
	throw BoutException(
	      "Inconsistent default values for '{:s}': '{:s}' then '{:s}'", full_name,
	      bout::utils::variantToString(value),
	      bout::utils::variantToString(def.value));
      }
    }
  }
  return *this;
}


template <typename T> T Options::withDefault(T def) const {
  if (!is_value) {
    // Option not found
    output_info << _("\tOption ") << full_name << " = " << def << " (" << DEFAULT_SOURCE
		<< ")" << std::endl;
    return def;
  }
  T val = as<T>(def);
  // Check if this was previously set as a default option
  if (bout::utils::variantEqualTo(attributes.at("source"), DEFAULT_SOURCE)) {
    // Check that the default values are the same
    if (!similar(val, def)) {
      throw BoutException("Inconsistent default values for '{:s}': '{:s}' then '{:s}'",
			  full_name, bout::utils::variantToString(value),
			  toString(def));
    }
  }
  return val;
}

template <typename T>
T Options::as(const T& UNUSED(similar_to)) const {
  if (!is_value) {
    throw BoutException("Option {:s} has no value", full_name);
  }

  T val;

  // Try casting. This will throw std::bad_cast if it can't be done
  try {
    val = bout::utils::variantStaticCastOrThrow<ValueType, T>(value);
  } catch (const std::bad_cast &e) {
    // If the variant is a string then we may be able to parse it

    if (bout::utils::holds_alternative<std::string>(value)) {
      std::stringstream ss(bout::utils::get<std::string>(value));
      ss >> val;

      // Check if the parse failed
      if (ss.fail()) {
	throw BoutException("Option {:s} could not be parsed ('{:s}')", full_name,
			    bout::utils::variantToString(value));
      }

      // Check if there are characters remaining
      std::string remainder;
      std::getline(ss, remainder);
      for (const char &ch : remainder) {
	if (!std::isspace(static_cast<unsigned char>(ch))) {
	  // Meaningful character not parsed
	  throw BoutException("Option {:s} could not be parsed", full_name);
	}
      }
    } else {
      // Another type which can't be casted
      throw BoutException("Option {:s} could not be converted to type {:s}", full_name,
			  typeid(T).name());
    }
  }

  // Mark this option as used
  value_used = true; // Note this is mutable

  output_info << "\tOption " << full_name  << " = " << val;
  if (attributes.count("source")) {
    // Specify the source of the setting
    output_info << " (" << bout::utils::variantToString(attributes.at("source")) << ")";
  }
  output_info << endl;

  return val;
}
