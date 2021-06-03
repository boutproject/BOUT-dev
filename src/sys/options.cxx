#include <boutexception.hxx>
#include <field_factory.hxx> // Used for parsing expressions
#include <options.hxx>
#include <output.hxx>
#include <utils.hxx>

#include <fmt/format.h>

#include <algorithm>
#include <iomanip>
#include <sstream>

/// The source label given to default values
const std::string Options::DEFAULT_SOURCE{_("default")};

/// Name of the attribute to indicate an Option should always count as
/// having been used
constexpr auto conditionally_used_attribute = "conditionally used";

Options *Options::root_instance{nullptr};

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

Options::Options(const Options& other)
    : value(other.value), attributes(other.attributes),
      parent_instance(other.parent_instance), full_name(other.full_name),
      is_section(other.is_section), children(other.children), value_used(other.value_used) {

  // Ensure that this is the parent of all children,
  // otherwise will point to the original Options instance
  for (auto& child : children) {
    child.second.parent_instance = this;
  }
}

template <>
Options::Options(const char* value) {
  assign<std::string>(value);
}

Options::Options(std::initializer_list<std::pair<std::string, Options>> values) {
  // Yes, this looks bad, but bear with me...  The _only_ way to
  // construct a nested initializer_list is inside-out, from the
  // bottom of the tree structure. Unfortunately, this is very much
  // not how you would want to construct the tree of options, as we
  // don't have the parent section's name as we construct each
  // child. Therefore, when we _do_ construct the parent, we'll need
  // to recursively step down the tree, prepending the parent's name
  // to each child. Rather than have a private member to do that, we
  // use a lambda. And to make that lambda recursive, we need to have
  // a nested lambda.
  auto append_section_name = [](auto& children, const std::string& section_name) {
    auto append_impl = [](auto& children, const std::string& section_name, auto& append_ref) mutable -> void {
      for (auto& child : children) {
        child.second.full_name = fmt::format("{}:{}", section_name, child.second.full_name);
        if (child.second.is_section) {
          append_ref(child.second.children, section_name, append_ref);
        }
      }
    };
    append_impl(children, section_name, append_impl);
  };

  for (auto& value : values) {
    (*this)[value.first] = value.second;
    // value.second was constructed from the "bare" `Options<T>(T)` so
    // doesn't have `full_name` set. This clobbers
    // `(*this)[value.first].full_name` in the copy constructor, so we
    // need to explicitly set it again
    (*this)[value.first].full_name = value.first;
    append_section_name((*this)[value.first].children, value.first);
  }
}

Options& Options::operator[](const std::string& name) {
  TRACE("Options::operator[]");

  if (isValue()) {
    throw BoutException(
        _("Trying to index Option '{0}' with '{1}', but '{0}' is a value, not a section.\n"
          "This is likely the result of clashing input options, and you may have to rename one of them.\n"),
        full_name, name);
  }

  if (name.empty()) {
    return *this;
  }

  // If name is compound, e.g. "section:subsection", then split the name
  auto subsection_split = name.find(":");
  if (subsection_split != std::string::npos) {
    return (*this)[name.substr(0, subsection_split)][name.substr(subsection_split + 1)];
  }

  // Find and return if already exists
  auto it = children.find(name);
  if (it != children.end()) {
    return it->second;
  }

  // Doesn't exist yet, so add
  std::string secname = name;
  if (!full_name.empty()) { // prepend the section name
    secname = full_name + ":" + secname;
  }

  // emplace returns a pair with iterator first, boolean (insert yes/no) second
  auto pair_it = children.emplace(name, Options{this, secname});

  return pair_it.first->second;
}

const Options& Options::operator[](const std::string& name) const {
  TRACE("Options::operator[] const");

  if (isValue()) {
    throw BoutException(
        _("Trying to index Option '{0}' with '{1}', but '{0}' is a value, not a section.\n"
          "This is likely the result of clashing input options, and you may have to rename one of them.\n"),
        full_name, name);
  }

  if (name.empty()) {
    return *this;
  }

  // If name is compound, e.g. "section:subsection", then split the name
  auto subsection_split = name.find(":");
  if (subsection_split != std::string::npos) {
    return (*this)[name.substr(0, subsection_split)][name.substr(subsection_split + 1)];
  }

  // Find and return if already exists
  auto it = children.find(name);
  if (it == children.end()) {
    // Doesn't exist
    throw BoutException(_("Option {:s}:{:s} does not exist"), full_name, name);
  }

  return it->second;
}

std::multiset<Options::FuzzyMatch>
Options::fuzzyFind(const std::string& name, std::string::size_type distance) const {
  std::multiset<Options::FuzzyMatch> matches;

  // Add option.full_name to matches if possible_match is within
  // distance of name, including a possible extra_cost. Returns true
  // if it was a fuzzy match
  auto insert_if_match = [&](const Options& option, const std::string& possible_match,
                             std::string::size_type extra_cost = 0) -> bool {
    if (not option.isValue()) {
      // Don't match section names
      return false;
    }
    if ((name != possible_match) and (lowercase(name) == lowercase(possible_match))) {
      // Differs only in case: pretty good match
      matches.insert({option, 1 + extra_cost});
      return true;
    }
    const auto fuzzy_distance = editDistance(name, possible_match) + extra_cost;
    if (fuzzy_distance <= distance) {
      // Insert the full_name with parent sections, not the possible_match
      matches.insert({option, fuzzy_distance});
      return true;
    }
    return false;
  };

  insert_if_match(*this, full_name);

  for (const auto& child : children) {
    // Check against fully-qualified name
    if (not insert_if_match(child.second, child.second.full_name)) {
      // Try again without parent sections, but make it cost a little bit more.
      // Could make it cost more per wrong section by counting ":" in full_name
      insert_if_match(child.second, child.second.name(), 1);
    }

    if (child.second.is_section) {
      // Recurse down the tree
      auto child_matches = child.second.fuzzyFind(name, distance);
      matches.insert(child_matches.begin(), child_matches.end());
    }
  }

  return matches;
}

Options& Options::operator=(const Options& other) {
  // Note: Here can't do copy-and-swap because pointers to parents are stored

  value = other.value;
  attributes = other.attributes;
  full_name = other.full_name;
  is_section = other.is_section;
  children = other.children;
  value_used = other.value_used;

  // Ensure that this is the parent of all children,
  // otherwise will point to the original Options instance
  for (auto& child : children) {
    child.second.parent_instance = this;
  }
  return *this;
}

bool Options::isSet() const {
  // Only values can be set/unset
  if (is_section) {
    return false;
  }

  // Ignore if set from default
  if (bout::utils::variantEqualTo(attributes.at("source"), DEFAULT_SOURCE)) {
    return false;
  }

  return true;
}

bool Options::isSection(const std::string& name) const {
  if (name == "") {
    // Test this object
    return is_section;
  }

  // Is there a child section?
  auto it = children.find(name);
  if (it == children.end()) {
    return false;
  } else {
    return it->second.isSection();
  }
}

template <> std::string Options::as<std::string>(const std::string& UNUSED(similar_to)) const {
  if (is_section) {
    throw BoutException(_("Option {:s} has no value"), full_name);
  }

  // Mark this option as used
  value_used = true;

  std::string result = bout::utils::variantToString(value);
  
  output_info << _("\tOption ") << full_name << " = " << result;
  if (attributes.count("source")) {
    // Specify the source of the setting
    output_info << " (" << bout::utils::variantToString(attributes.at("source")) << ")";
  }
  output_info << endl;

  return result;
}

template <> int Options::as<int>(const int& UNUSED(similar_to)) const {
  if (is_section) {
    throw BoutException(_("Option {:s} has no value"), full_name);
  }

  int result;

  if (bout::utils::holds_alternative<int>(value)) {
    result = bout::utils::get<int>(value);
    
  } else {
    // Cases which get a BoutReal then check if close to an integer
    BoutReal rval;
    
    if (bout::utils::holds_alternative<BoutReal>(value)) {
      rval = bout::utils::get<BoutReal>(value);
    
    } else if (bout::utils::holds_alternative<std::string>(value)) {
      // Use FieldFactory to evaluate expression
      // Parse the string, giving this Option pointer for the context
      // then generate a value at t,x,y,z = 0,0,0,0
      auto gen = FieldFactory::get()->parse(bout::utils::get<std::string>(value), this);
      if (!gen) {
        throw BoutException(_("Couldn't get integer from option {:s} = '{:s}'"),
                            full_name, bout::utils::variantToString(value));
      }
      rval = gen->generate({});
    } else {
      // Another type which can't be converted
      throw BoutException(_("Value for option {:s} is not an integer"), full_name);
    }
    
    // Convert to int by rounding
    result = ROUND(rval);
    
    // Check that the value is close to an integer
    if (fabs(rval - static_cast<BoutReal>(result)) > 1e-3) {
      throw BoutException(_("Value for option {:s} = {:e} is not an integer"), full_name,
                          rval);
    }
  }

  value_used = true;

  output_info << _("\tOption ") << full_name << " = " << result;
  if (attributes.count("source")) {
    // Specify the source of the setting
    output_info << " (" << bout::utils::variantToString(attributes.at("source")) << ")";
  }
  output_info << endl;

  return result;
}

template <> BoutReal Options::as<BoutReal>(const BoutReal& UNUSED(similar_to)) const {
  if (is_section) {
    throw BoutException(_("Option {:s} has no value"), full_name);
  }

  BoutReal result;
  
  if (bout::utils::holds_alternative<int>(value)) {
    result = static_cast<BoutReal>(bout::utils::get<int>(value));
    
  } else if (bout::utils::holds_alternative<BoutReal>(value)) {
    result = bout::utils::get<BoutReal>(value);
      
  } else if (bout::utils::holds_alternative<std::string>(value)) {
    
    // Use FieldFactory to evaluate expression
    // Parse the string, giving this Option pointer for the context
    // then generate a value at t,x,y,z = 0,0,0,0
    auto gen = FieldFactory::get()->parse(bout::utils::get<std::string>(value), this);
    if (!gen) {
      throw BoutException(_("Couldn't get BoutReal from option {:s} = '{:s}'"), full_name,
                          bout::utils::get<std::string>(value));
    }
    result = gen->generate({});
  } else {
    throw BoutException(_("Value for option {:s} cannot be converted to a BoutReal"),
                        full_name);
  }
  
  // Mark this option as used
  value_used = true;
  
  output_info << _("\tOption ") << full_name << " = " << result;
  if (attributes.count("source")) {
    // Specify the source of the setting
    output_info << " (" << bout::utils::variantToString(attributes.at("source")) << ")";
  }
  output_info << endl;
  
  return result;
}

template <> bool Options::as<bool>(const bool& UNUSED(similar_to)) const {
  if (is_section) {
    throw BoutException(_("Option {:s} has no value"), full_name);
  }
  
  bool result;
  
  if (bout::utils::holds_alternative<bool>(value)) {
    result = bout::utils::get<bool>(value);
  
  } else if(bout::utils::holds_alternative<std::string>(value)) {
    // case-insensitve check, so convert string to lower case
    const auto strvalue = lowercase(bout::utils::get<std::string>(value));
  
    if ((strvalue == "y") or (strvalue == "yes") or (strvalue == "t")
        or (strvalue == "true") or (strvalue == "1")) {
      result = true;
    } else if ((strvalue == "n") or (strvalue == "no") or (strvalue == "f")
        or (strvalue == "false") or (strvalue == "0")) {
      result = false;
    } else {
      throw BoutException(_("\tOption '{:s}': Boolean expected. Got '{:s}'\n"), full_name,
                          strvalue);
    }
  } else {
    throw BoutException(_("Value for option {:s} cannot be converted to a bool"),
                        full_name);
  }
  
  value_used = true;
  
  output_info << _("\tOption ") << full_name << " = " << toString(result);
  
  if (attributes.count("source")) {
    // Specify the source of the setting
    output_info << " (" << bout::utils::variantToString(attributes.at("source")) << ")";
  }
  output_info << endl;

  return result;
}

template <> Field3D Options::as<Field3D>(const Field3D& similar_to) const {
  if (is_section) {
    throw BoutException("Option {:s} has no value", full_name);
  }

  // Mark value as used
  value_used = true;

  if (bout::utils::holds_alternative<Field3D>(value)) {
    Field3D stored_value = bout::utils::get<Field3D>(value);
    
    // Check that meta-data is consistent
    ASSERT1_FIELDS_COMPATIBLE(stored_value, similar_to);
    
    return stored_value;
  }

  if (bout::utils::holds_alternative<Field2D>(value)) {
    const auto& stored_value = bout::utils::get<Field2D>(value);

    // Check that meta-data is consistent
    ASSERT1_FIELDS_COMPATIBLE(stored_value, similar_to);

    return Field3D(stored_value);
  }
  
  try {
    BoutReal scalar_value = bout::utils::variantStaticCastOrThrow<ValueType, BoutReal>(value);
    
    // Get metadata from similar_to, fill field with scalar_value
    return filledFrom(similar_to, scalar_value);
  } catch (const std::bad_cast&) {
    
    // Convert from a string using FieldFactory
    if (bout::utils::holds_alternative<std::string>(value)) {
      return FieldFactory::get()->create3D(bout::utils::get<std::string>(value), this,
                                           similar_to.getMesh(),
                                           similar_to.getLocation());
    } else if (bout::utils::holds_alternative<Tensor<BoutReal>>(value)) {
      auto localmesh = similar_to.getMesh();
      if (!localmesh) {
        throw BoutException("mesh must be supplied when converting Tensor to Field3D");
      }

      // Get a reference, to try and avoid copying
      const auto& tensor = bout::utils::get<Tensor<BoutReal>>(value);
      
      // Check if the dimension sizes are the same as a Field3D
      if (tensor.shape() == std::make_tuple(localmesh->LocalNx,
                                            localmesh->LocalNy,
                                            localmesh->LocalNz)) {
        return Field3D(tensor.getData(), localmesh, similar_to.getLocation(),
                       {similar_to.getDirectionY(), similar_to.getDirectionZ()});
      }
      // If dimension sizes not the same, may be able
      // to select a region from it using Mesh e.g. if this
      // is from the input grid file.

    }
  }
  throw BoutException(_("Value for option {:s} cannot be converted to a Field3D"),
                      full_name);
}

template <> Field2D Options::as<Field2D>(const Field2D& similar_to) const {
  if (is_section) {
    throw BoutException("Option {:s} has no value", full_name);
  }
  
  // Mark value as used
  value_used = true;

  if (bout::utils::holds_alternative<Field2D>(value)) {
    Field2D stored_value = bout::utils::get<Field2D>(value);
    
    // Check that meta-data is consistent
    ASSERT1_FIELDS_COMPATIBLE(stored_value, similar_to);

    return stored_value;
  }
  
  try {
    BoutReal scalar_value = bout::utils::variantStaticCastOrThrow<ValueType, BoutReal>(value);

    // Get metadata from similar_to, fill field with scalar_value
    return filledFrom(similar_to, scalar_value);
  } catch (const std::bad_cast&) {
    
    // Convert from a string using FieldFactory
    if (bout::utils::holds_alternative<std::string>(value)) {
      return FieldFactory::get()->create2D(bout::utils::get<std::string>(value), this,
                                           similar_to.getMesh(),
                                           similar_to.getLocation());
    } else if (bout::utils::holds_alternative<Matrix<BoutReal>>(value)) {
      auto localmesh = similar_to.getMesh();
      if (!localmesh) {
        throw BoutException("mesh must be supplied when converting Matrix to Field2D");
      }

      // Get a reference, to try and avoid copying
      const auto& matrix = bout::utils::get<Matrix<BoutReal>>(value);

      // Check if the dimension sizes are the same as a Field3D
      if (matrix.shape() == std::make_tuple(localmesh->LocalNx,
                                            localmesh->LocalNy)) {
        return Field2D(matrix.getData(), localmesh, similar_to.getLocation(),
                       {similar_to.getDirectionY(), similar_to.getDirectionZ()});
      }
    }
  }
  throw BoutException(_("Value for option {:s} cannot be converted to a Field2D"),
                      full_name);
}

template <>
FieldPerp Options::as<FieldPerp>(const FieldPerp& similar_to) const {
  if (is_section) {
    throw BoutException("Option {:s} has no value", full_name);
  }

  // Mark value as used
  value_used = true;

  if (bout::utils::holds_alternative<FieldPerp>(value)) {
    FieldPerp stored_value = bout::utils::get<FieldPerp>(value);

    // Check that meta-data is consistent
    ASSERT1_FIELDS_COMPATIBLE(stored_value, similar_to);

    return stored_value;
  }

  try {
    BoutReal scalar_value =
        bout::utils::variantStaticCastOrThrow<ValueType, BoutReal>(value);

    // Get metadata from similar_to, fill field with scalar_value
    return filledFrom(similar_to, scalar_value);
  } catch (const std::bad_cast&) {

    const CELL_LOC location = hasAttribute("cell_location")
                                  ? CELL_LOCFromString(attributes.at("cell_location"))
                                  : similar_to.getLocation();

    // Convert from a string using FieldFactory
    if (bout::utils::holds_alternative<std::string>(value)) {
      return FieldFactory::get()->createPerp(bout::utils::get<std::string>(value), this,
                                             similar_to.getMesh(), location);
    } else if (bout::utils::holds_alternative<Matrix<BoutReal>>(value)) {
      auto localmesh = similar_to.getMesh();
      if (!localmesh) {
        throw BoutException("mesh must be supplied when converting Matrix to FieldPerp");
      }

      // Get a reference, to try and avoid copying
      const auto& matrix = bout::utils::get<Matrix<BoutReal>>(value);

      // Check if the dimension sizes are the same as a FieldPerp
      if (matrix.shape() == std::make_tuple(localmesh->LocalNx, localmesh->LocalNz)) {
        const auto y_direction =
            hasAttribute("direction_y")
                ? YDirectionTypeFromString(attributes.at("direction_y"))
                : similar_to.getDirectionY();
        const auto z_direction =
            hasAttribute("direction_z")
                ? ZDirectionTypeFromString(attributes.at("direction_z"))
                : similar_to.getDirectionZ();

        auto result = FieldPerp(matrix.getData(), localmesh, location, -1,
                                {y_direction, z_direction});

        // Set the index after creating the field so as to not
        // duplicate the code in `FieldPerp::setIndexFromGlobal`
        if (hasAttribute("yindex_global")) {
          result.setIndexFromGlobal(attributes.at("yindex_global"));
        } else if (similar_to.getIndex() == -1) {
          // If `yindex_global` attribute wasn't present (might be an
          // older file), and `similar_to` doesn't have its index set
          // (might not have been passed, so be default constructed),
          // use the no-boundary form so that we get a default value
          // on a grid cell
          result.setIndex(localmesh->getLocalYIndexNoBoundaries(0));
        } else {
          result.setIndex(similar_to.getIndex());
        }
        return result;
      }
      // If dimension sizes not the same, may be able
      // to select a region from it using Mesh e.g. if this
      // is from the input grid file.
    }
  }
  throw BoutException(_("Value for option {:s} cannot be converted to a Field3D"),
                      full_name);
}

// Note: This is defined here rather than in the header
// to avoid using as<string> before specialising it.
bool Options::operator==(const char* other) const {
  return as<std::string>() == std::string(other);
}

bool Options::operator<(const char* other) const {
  return as<std::string>() < std::string(other);
}

Options Options::getUnused(const std::vector<std::string>& exclude_sources) const {
  // Check if the option should count as having been used due to its source
  const auto has_excluded_source = [&exclude_sources](const Options& option) -> bool {
    if (not option.hasAttribute("source")) {
      return false;
    }
    const auto source = option.attributes.at("source").as<std::string>();
    return std::find(exclude_sources.begin(), exclude_sources.end(), source)
           != exclude_sources.end();
  };

  const auto conditionally_used = [](const Options& option) -> bool {
    if (not option.hasAttribute(conditionally_used_attribute)) {
      return false;
    }
    return option.attributes.at(conditionally_used_attribute).as<bool>();
  };

  // Copy this object, and then we're going to chuck out everything
  // that has been used. This turns out to be easier than copying just
  // the unused options into an empty instance
  Options unused = *this;

  if (unused.isValue()) {
    // If this is from an excluded source, count it as being used
    if (has_excluded_source(unused) or conditionally_used(unused)) {
      unused.value_used = true;
    }
    // We don't have a nice way to "clear" the value, so if it was
    // used, mark it as no longer a value: if it has been used, this
    // does nothing
    unused.is_section = unused.value_used;
    return unused;
  }

  // This loop modifies the map in the loop, so we need to manually
  // manage the iterator
  for (auto child = unused.children.begin(); child != unused.children.end();) {
    // Remove the child if it's been used or if it's from a source we
    // should count as having been used
    if (child->second.isValue()
        and (child->second.value_used or has_excluded_source(child->second)
             or conditionally_used(child->second))) {
      child = unused.children.erase(child);
      continue;
    }

    if (child->second.is_section) {
      // Recurse down and replace this section by its "unused" version
      child->second = child->second.getUnused();
      // If all of its children have been used, then we can remove it
      // as well
      if (child->second.children.empty()) {
        child = unused.children.erase(child);
        continue;
      }
    }

    ++child;
  }

  if (unused.children.empty()) {
    // If all the children have been used, we don't want to print a
    // section name any more
    unused.full_name.clear();
  }

  return unused;
}

void Options::printUnused() const {
  Options unused = getUnused();

  // Two cases: single value, or a section.  If it's a single value,
  // we can check it directly. If it's a section, we can see if it has
  // any children
  if ((unused.isValue() and unused.value_used) or unused.children.empty()) {
    output_info << _("All options used\n");
    return;
  }

  output_info << _("Unused options:\n") << unused;
}

void Options::setConditionallyUsed() {
  attributes[conditionally_used_attribute] = true;
  for (auto& child : children) {
    child.second.setConditionallyUsed();
  }
}

void Options::cleanCache() { FieldFactory::get()->cleanCache(); }

std::map<std::string, Options::OptionValue> Options::values() const {
  std::map<std::string, OptionValue> options;
  for (const auto& it : children) {
    if (it.second.isValue()) {
      options.emplace(it.first, OptionValue { bout::utils::variantToString(it.second.value),
                                               bout::utils::variantToString(it.second.attributes.at("source")),
                                               it.second.value_used});
    }
  }
  return options;
}

std::map<std::string, const Options *> Options::subsections() const {
  std::map<std::string, const Options *> sections;
  for (const auto &it : children) {
    if (it.second.is_section) {
      sections[it.first] = &it.second;
    }
  }
  return sections;
}

std::vector<std::string> Options::getFlattenedKeys() const {
  std::vector<std::string> flattened_names;

  if (isValue() and not full_name.empty()) {
    flattened_names.push_back(full_name);
  }

  for (const auto& child : children) {
    if (child.second.isValue()) {
      flattened_names.push_back(child.second.full_name);
    }
    if (child.second.is_section) {
      const auto child_names = child.second.getFlattenedKeys();
      flattened_names.insert(flattened_names.end(), child_names.begin(),
                             child_names.end());
    }
  }

  return flattened_names;
}

fmt::format_parse_context::iterator
bout::details::OptionsFormatterBase::parse(fmt::format_parse_context& ctx) {

  auto it = ctx.begin();
  const auto end = ctx.end();
  while (it != end and *it != '}') {
    switch (*it) {
    case 'd':
      docstrings = true;
      break;
    case 'i':
      inline_section_names = true;
      break;
    case 'k':
      key_only = true;
      break;
    case 's':
      source = true;
      break;
    default:
      throw fmt::format_error("invalid format");
    }
    ++it;
  }

  // Keep a copy of the format string (without the last '}') so we can
  // pass it down to the subsections.
  const auto size = std::distance(ctx.begin(), it);
  format_string.reserve(size + 3);
  format_string.assign("{:");
  format_string.append(ctx.begin(), it);
  format_string.push_back('}');

  return it;
}

fmt::format_context::iterator
bout::details::OptionsFormatterBase::format(const Options& options,
                                            fmt::format_context& ctx) {
  const auto children = options.getChildren();
  const bool has_child_values =
      std::any_of(children.begin(), children.end(),
                  [](const auto& child) { return child.second.isValue(); });

  // Only print section headers if the section has a name and it has
  // non-section children
  const std::string section_name = options.str();
  if (not inline_section_names and not section_name.empty() and has_child_values) {
    fmt::format_to(ctx.out(), "\n[{}]\n", section_name);
  }

  // Get all the child values first
  for (const auto& child : children) {
    if (child.second.isValue()) {
      if (inline_section_names and not section_name.empty()) {
        fmt::format_to(ctx.out(), "{}:", section_name);
      }

      fmt::format_to(ctx.out(), "{}", child.first);

      if (not key_only) {
        const auto value = bout::utils::variantToString(child.second.value);
        // Convert empty strings to ""
        const std::string as_str = value.empty() ? "\"\"" : value;
        fmt::format_to(ctx.out(), " = {}", as_str);
      }

      const bool has_doc = child.second.attributes.count("doc");
      const bool has_source = child.second.attributes.count("source");
      const bool has_type = child.second.attributes.count("type");

      std::vector<std::string> comments;

      if (docstrings) {
        if (has_type) {
          comments.emplace_back(fmt::format(
              "type: {}", child.second.attributes.at("type").as<std::string>()));
        }

        if (has_doc) {
          comments.emplace_back(fmt::format(
              "doc: {}", child.second.attributes.at("doc").as<std::string>()));
        }
      }

      if (source and has_source) {
        const auto source = child.second.attributes.at("source").as<std::string>();
        if (not source.empty()) {
          comments.emplace_back(fmt::format("source: {}", source));
        }
      }

      if (not comments.empty()) {
        fmt::format_to(ctx.out(), "\t\t# {}", fmt::join(comments, ", "));
      }

      fmt::format_to(ctx.out(), "\n");
    }
  }

  // Now descend the tree, accumulating subsections
  for (const auto& subsection : options.subsections()) {
    fmt::format_to(ctx.out(), format_string, *subsection.second);
  }

  return ctx.out();
}

std::string toString(const Options& value) { return fmt::format("{}", value); }

namespace bout {
void checkForUnusedOptions() {
  auto& options = Options::root();
  const bool error_on_unused_options =
      options["input"]["error_on_unused_options"]
          .doc(
              "Error if there are any unused options before starting the main simulation")
          .withDefault(true);

  if (not error_on_unused_options) {
    return;
  }
  checkForUnusedOptions(options, options["datadir"].withDefault("data"),
                        options["optionfile"].withDefault("BOUT.inp"));
}

void checkForUnusedOptions(const Options& options, const std::string& data_dir,
                           const std::string& option_file) {
  Options unused = options.getUnused();
  if (not unused.getChildren().empty()) {

    // Construct a string with all the fuzzy matches for each unused option
    const auto keys = unused.getFlattenedKeys();
    std::string possible_misspellings;
    for (const auto& key : keys) {
      auto fuzzy_matches = options.fuzzyFind(key);
      // Remove unacceptable matches, including:
      // - exact matches
      // - other unused options
      // - options set internally by the library and not meant as user inputs
      bout::utils::erase_if(fuzzy_matches, [](const Options::FuzzyMatch& match) -> bool {
        const auto source = match.match.hasAttribute("source")
                                ? match.match.attributes.at("source").as<std::string>()
                                : "";
        const bool internal_source = (source == "Solver") or (source == "Output");

        return match.distance == 0 or (not match.match.valueUsed()) or internal_source;
      });

      if (fuzzy_matches.empty()) {
        continue;
      }
      possible_misspellings += fmt::format("\nUnused option '{}', did you mean:\n", key);
      for (const auto& match : fuzzy_matches) {
        possible_misspellings += fmt::format("\t{}\n", match.match.str());
      }
    }

    // Only display the possible matches if we actually have some to show
    const std::string additional_info =
        possible_misspellings.empty()
            ? ""
            : fmt::format("Suggested alternatives:\n{}", possible_misspellings);

    // Raw string to help with the formatting of the message, and a
    // separate variable so clang-format doesn't barf on the
    // exception
    const std::string unused_message = _(R"""(
There were unused input options:
-----
{:i}
-----
It's possible you've mistyped some options. BOUT++ input arguments are
now case-sensitive, and some have changed name. You can try running

    <BOUT++ directory>/bin/bout-v5-input-file-upgrader.py {}/{}

to automatically fix the most common issues. If these options above
are sometimes used depending on other options, you can call
`Options::setConditionallyUsed()`, for example:

    Options::root()["{}"].setConditionallyUsed();

to mark a section or value as depending on other values, and so ignore
it in this check. Alternatively, if you're sure the above inputs are
not a mistake, you can set 'input:error_on_unused_options=false' to
turn off this check for unused options. You can always set
'input:validate=true' to check inputs without running the full
simulation.

{})""");

    throw BoutException(unused_message, unused, data_dir, option_file,
                        unused.getChildren().begin()->first, additional_info);
  }
}
} // namespace bout
