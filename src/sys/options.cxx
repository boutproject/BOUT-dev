#include "bout/options.hxx"

#include "bout/array.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/field_factory.hxx" // Used for parsing expressions
#include "bout/fieldperp.hxx"
#include "bout/output.hxx"
#include "bout/sys/expressionparser.hxx"
#include "bout/sys/gettext.hxx"
#include "bout/sys/type_name.hxx"
#include "bout/sys/variant.hxx"
#include "bout/traits.hxx"
#include "bout/unused.hxx"
#include "bout/utils.hxx"

#include <fmt/core.h>
#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

/// The source label given to default values
const std::string Options::DEFAULT_SOURCE{_("default")};

std::string Options::getDefaultSource() { return DEFAULT_SOURCE; }

/// Name of the attribute to indicate an Option should always count as
/// having been used
constexpr auto conditionally_used_attribute = "conditionally used";

Options& Options::root() {
  static Options root_instance;
  return root_instance;
}

void Options::cleanup() { root() = Options{}; }

Options Options::copy() const {
  Options result;

  result.value = value;
  result.attributes = attributes;
  result.parent_instance = parent_instance;
  result.full_name = full_name;
  result.is_section = is_section;
  result.value_used = value_used;

  // Recursively copy children
  for (const auto& child_it : children) {
    auto pair_it = result.children.emplace(child_it.first, child_it.second.copy());
    Options& child = pair_it.first->second;
    child.parent_instance = &result;
  }
  return result;
}

Options::Options(Options&& other) noexcept
    : value(std::move(other.value)), attributes(std::move(other.attributes)),
      parent_instance(other.parent_instance), full_name(std::move(other.full_name)),
      is_section(other.is_section), children(std::move(other.children)),
      value_used(other.value_used) {

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

Options::Options(InitializerList values, Options* parent_instance,
                 std::string section_name)
    : parent_instance(parent_instance), full_name(std::move(section_name)),
      is_section(true) {
  for (const auto& value_it : values) {
    std::string child_name = full_name.empty()
                                 ? value_it.first
                                 : fmt::format("{}:{}", full_name, value_it.first);
    if (value_it.second.children.size() != 0) {
      // A section, so construct with an initializer_list
      children.emplace(value_it.first,
                       Options(value_it.second.children, this, std::move(child_name)));
    } else {
      // A value
      auto pair_it =
          children.emplace(value_it.first, Options(this, std::move(child_name)));
      Options& child = pair_it.first->second;
      child._set_no_check(value_it.second.value, "");
    }
  }
}

Options& Options::operator[](const std::string& name) {
  TRACE("Options::operator[]");

  if (isValue()) {
    throw BoutException(_("Trying to index Option '{0}' with '{1}', but '{0}' is a "
                          "value, not a section.\n"
                          "This is likely the result of clashing input options, and you "
                          "may have to rename one of them.\n"),
                        full_name, name);
  }

  if (name.empty()) {
    return *this;
  }

  // If name is compound, e.g. "section:subsection", then split the name
  auto subsection_split = name.find(':');
  if (subsection_split != std::string::npos) {
    return (*this)[name.substr(0, subsection_split)][name.substr(subsection_split + 1)];
  }

  // Find and return if already exists
  auto child = children.find(name);
  if (child != children.end()) {
    return child->second;
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
    throw BoutException(_("Trying to index Option '{0}' with '{1}', but '{0}' is a "
                          "value, not a section.\n"
                          "This is likely the result of clashing input options, and you "
                          "may have to rename one of them.\n"),
                        full_name, name);
  }

  if (name.empty()) {
    return *this;
  }

  // If name is compound, e.g. "section:subsection", then split the name
  auto subsection_split = name.find(':');
  if (subsection_split != std::string::npos) {
    return (*this)[name.substr(0, subsection_split)][name.substr(subsection_split + 1)];
  }

  // Find and return if already exists
  auto child = children.find(name);
  if (child == children.end()) {
    // Doesn't exist
    throw BoutException(_("Option {:s}:{:s} does not exist"), full_name, name);
  }

  return child->second;
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

Options& Options::operator=(Options&& other) noexcept {
  if (this == &other) {
    return *this;
  }

  // Note: Here can't do copy-and-swap because pointers to parents are stored

  value = std::move(other.value);
  attributes = std::move(other.attributes);
  full_name = std::move(other.full_name);
  is_section = other.is_section;
  children = std::move(other.children);
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
  if (name.empty()) {
    // Test this object
    return is_section;
  }

  // Is there a child section?
  const auto child = children.find(name);
  if (child == children.end()) {
    return false;
  }
  return child->second.isSection();
}

template <>
void Options::assign<>(Field2D val, std::string source) {
  attributes["cell_location"] = toString(val.getLocation());
  attributes["direction_y"] = toString(val.getDirectionY());
  attributes["direction_z"] = toString(val.getDirectionZ());
  _set_no_check(std::move(val), std::move(source));
}
template <>
void Options::assign<>(Field3D val, std::string source) {
  attributes["cell_location"] = toString(val.getLocation());
  attributes["direction_y"] = toString(val.getDirectionY());
  attributes["direction_z"] = toString(val.getDirectionZ());
  _set_no_check(std::move(val), std::move(source));
}
template <>
void Options::assign<>(FieldPerp val, std::string source) {
  attributes["cell_location"] = toString(val.getLocation());
  attributes["direction_y"] = toString(val.getDirectionY());
  attributes["direction_z"] = toString(val.getDirectionZ());
  attributes["yindex_global"] = val.getGlobalIndex();
  _set_no_check(std::move(val), std::move(source));
}
template <>
void Options::assign<>(Array<BoutReal> val, std::string source) {
  _set_no_check(std::move(val), std::move(source));
}
template <>
void Options::assign<>(Matrix<BoutReal> val, std::string source) {
  _set_no_check(std::move(val), std::move(source));
}
template <>
void Options::assign<>(Tensor<BoutReal> val, std::string source) {
  _set_no_check(std::move(val), std::move(source));
}

namespace {
/// Use FieldFactory to evaluate expression
double parseExpression(const Options::ValueType& value, const Options* options,
                       const std::string& type, const std::string& full_name) {
  try {
    // Parse the string, giving this Option pointer for the context
    // then generate a value at t,x,y,z = 0,0,0,0
    auto gen = FieldFactory::get()->parse(bout::utils::get<std::string>(value), options);
    if (!gen) {
      throw ParseException("FieldFactory did not return a generator for '{}'",
                           bout::utils::variantToString(value));
    }
    return gen->generate({});
  } catch (ParseException& error) {
    // Convert any exceptions to something a bit more useful
    throw BoutException(_("Couldn't get {} from option {:s} = '{:s}': {}"), type,
                        full_name, bout::utils::variantToString(value), error.what());
  }
}

/// Helper function to print `key = value` with optional source
template <class T>
void printNameValueSourceLine(const Options& option, const T& value) {
  output_info.write(_("\tOption {} = {}"), option.str(), value);
  if (option.hasAttribute("source")) {
    // Specify the source of the setting
    output_info.write(" ({})",
                      bout::utils::variantToString(option.attributes.at("source")));
  }
  output_info.write("\n");
}
} // namespace

template <>
std::string Options::as<std::string>(const std::string& UNUSED(similar_to)) const {
  if (is_section) {
    throw BoutException(_("Option {:s} has no value"), full_name);
  }

  // Mark this option as used
  value_used = true;

  std::string result = bout::utils::variantToString(value);

  printNameValueSourceLine(*this, result);

  return result;
}

template <>
int Options::as<int>(const int& UNUSED(similar_to)) const {
  if (is_section) {
    throw BoutException(_("Option {:s} has no value"), full_name);
  }

  int result = 0;

  if (bout::utils::holds_alternative<int>(value)) {
    result = bout::utils::get<int>(value);

  } else {
    // Cases which get a BoutReal then check if close to an integer
    BoutReal rval = BoutNaN;

    if (bout::utils::holds_alternative<BoutReal>(value)) {
      rval = bout::utils::get<BoutReal>(value);

    } else if (bout::utils::holds_alternative<std::string>(value)) {
      rval = parseExpression(value, this, "integer", full_name);

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

  printNameValueSourceLine(*this, result);

  return result;
}

template <>
BoutReal Options::as<BoutReal>(const BoutReal& UNUSED(similar_to)) const {
  if (is_section) {
    throw BoutException(_("Option {:s} has no value"), full_name);
  }

  BoutReal result = BoutNaN;

  if (bout::utils::holds_alternative<int>(value)) {
    result = static_cast<BoutReal>(bout::utils::get<int>(value));

  } else if (bout::utils::holds_alternative<BoutReal>(value)) {
    result = bout::utils::get<BoutReal>(value);

  } else if (bout::utils::holds_alternative<std::string>(value)) {
    result = parseExpression(value, this, "BoutReal", full_name);

  } else {
    throw BoutException(_("Value for option {:s} cannot be converted to a BoutReal"),
                        full_name);
  }

  // Mark this option as used
  value_used = true;

  printNameValueSourceLine(*this, result);

  return result;
}

template <>
bool Options::as<bool>(const bool& UNUSED(similar_to)) const {
  if (is_section) {
    throw BoutException(_("Option {:s} has no value"), full_name);
  }

  bool result = false;

  if (bout::utils::holds_alternative<bool>(value)) {
    result = bout::utils::get<bool>(value);

  } else if (bout::utils::holds_alternative<std::string>(value)) {
    // Parse as floating point because that's the only type the parser understands
    const BoutReal rval = parseExpression(value, this, "bool", full_name);

    // Check that the result is either close to 1 (true) or close to 0 (false)
    const int ival = ROUND(rval);
    if ((fabs(rval - static_cast<BoutReal>(ival)) > 1e-3) or (ival < 0) or (ival > 1)) {
      throw BoutException(_("Value for option {:s} = {:e} is not a bool"), full_name,
                          rval);
    }
    result = ival == 1;
  } else {
    throw BoutException(_("Value for option {:s} cannot be converted to a bool"),
                        full_name);
  }

  value_used = true;

  printNameValueSourceLine(*this, toString(result));

  return result;
}

template <>
Field3D Options::as<Field3D>(const Field3D& similar_to) const {
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

  if (bout::utils::holds_alternative<BoutReal>(value)
      or bout::utils::holds_alternative<int>(value)) {
    const BoutReal scalar_value =
        bout::utils::variantStaticCastOrThrow<ValueType, BoutReal>(value);

    // Get metadata from similar_to, fill field with scalar_value
    return filledFrom(similar_to, scalar_value);
  }

  // Convert from a string using FieldFactory
  if (bout::utils::holds_alternative<std::string>(value)) {
    return FieldFactory::get()->create3D(bout::utils::get<std::string>(value), this,
                                         similar_to.getMesh(), similar_to.getLocation());
  }
  if (bout::utils::holds_alternative<Tensor<BoutReal>>(value)) {
    Mesh* localmesh = similar_to.getMesh();
    if (localmesh == nullptr) {
      throw BoutException("mesh must be supplied when converting Tensor to Field3D");
    }

    // Get a reference, to try and avoid copying
    const auto& tensor = bout::utils::get<Tensor<BoutReal>>(value);

    // Check if the dimension sizes are the same as a Field3D
    if (tensor.shape()
        == std::make_tuple(localmesh->LocalNx, localmesh->LocalNy, localmesh->LocalNz)) {
      return Field3D(tensor.getData(), localmesh, similar_to.getLocation(),
                     {similar_to.getDirectionY(), similar_to.getDirectionZ()});
    }
    // If dimension sizes not the same, may be able
    // to select a region from it using Mesh e.g. if this
    // is from the input grid file.
    const auto [tx, ty, tz] = tensor.shape();
    throw BoutException("Size mismatch for option {:s}: Tensor ({}, {}, {}) cannot be "
                        "converted to Field3D ({}, {}, {})",
                        full_name, tx, ty, tz, localmesh->LocalNx, localmesh->LocalNy,
                        localmesh->LocalNz);
  }

  throw BoutException(_("Value for option {:s} cannot be converted to a Field3D"),
                      full_name);
}

template <>
Field2D Options::as<Field2D>(const Field2D& similar_to) const {
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

  if (bout::utils::holds_alternative<BoutReal>(value)
      or bout::utils::holds_alternative<int>(value)) {
    const BoutReal scalar_value =
        bout::utils::variantStaticCastOrThrow<ValueType, BoutReal>(value);

    // Get metadata from similar_to, fill field with scalar_value
    return filledFrom(similar_to, scalar_value);
  }

  // Convert from a string using FieldFactory
  if (bout::utils::holds_alternative<std::string>(value)) {
    return FieldFactory::get()->create2D(bout::utils::get<std::string>(value), this,
                                         similar_to.getMesh(), similar_to.getLocation());
  }
  if (bout::utils::holds_alternative<Matrix<BoutReal>>(value)) {
    Mesh* localmesh = similar_to.getMesh();
    if (localmesh == nullptr) {
      throw BoutException("mesh must be supplied when converting Matrix to Field2D");
    }

    // Get a reference, to try and avoid copying
    const auto& matrix = bout::utils::get<Matrix<BoutReal>>(value);

    // Check if the dimension sizes are the same as a Field3D
    if (matrix.shape() == std::make_tuple(localmesh->LocalNx, localmesh->LocalNy)) {
      return Field2D(matrix.getData(), localmesh, similar_to.getLocation(),
                     {similar_to.getDirectionY(), similar_to.getDirectionZ()});
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

  if (bout::utils::holds_alternative<BoutReal>(value)
      or bout::utils::holds_alternative<int>(value)) {
    const BoutReal scalar_value =
        bout::utils::variantStaticCastOrThrow<ValueType, BoutReal>(value);

    // Get metadata from similar_to, fill field with scalar_value
    return filledFrom(similar_to, scalar_value);
  }
  const CELL_LOC location = hasAttribute("cell_location")
                                ? CELL_LOCFromString(attributes.at("cell_location"))
                                : similar_to.getLocation();

  // Convert from a string using FieldFactory
  if (bout::utils::holds_alternative<std::string>(value)) {
    return FieldFactory::get()->createPerp(bout::utils::get<std::string>(value), this,
                                           similar_to.getMesh(), location);
  }
  if (bout::utils::holds_alternative<Matrix<BoutReal>>(value)) {
    Mesh* localmesh = similar_to.getMesh();
    if (localmesh == nullptr) {
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
  throw BoutException(_("Value for option {:s} cannot be converted to a FieldPerp"),
                      full_name);
}

namespace {
/// Visitor to convert an int, BoutReal or Array/Matrix/Tensor to the
/// appropriate container
template <class Container>
struct ConvertContainer {
  ConvertContainer(std::string error, Container similar_to_)
      : error_message(std::move(error)), similar_to(std::move(similar_to_)) {}

  Container operator()(int value) {
    Container result(similar_to);
    std::fill(std::begin(result), std::end(result), value);
    return result;
  }

  Container operator()(BoutReal value) {
    Container result(similar_to);
    std::fill(std::begin(result), std::end(result), value);
    return result;
  }

  Container operator()(const Container& value) { return value; }

  template <class Other>
  Container operator()([[maybe_unused]] const Other& value) {
    throw BoutException(error_message);
  }

private:
  std::string error_message;
  Container similar_to;
};
} // namespace

template <>
Array<BoutReal> Options::as<Array<BoutReal>>(const Array<BoutReal>& similar_to) const {
  if (is_section) {
    throw BoutException(_("Option {:s} has no value"), full_name);
  }

  Array<BoutReal> result = bout::utils::visit(
      ConvertContainer<Array<BoutReal>>{
          fmt::format(
              _("Value for option {:s} cannot be converted to an Array<BoutReal>"),
              full_name),
          similar_to},
      value);

  // Mark this option as used
  value_used = true;

  printNameValueSourceLine(*this, "Array<BoutReal>");

  return result;
}

template <>
Matrix<BoutReal> Options::as<Matrix<BoutReal>>(const Matrix<BoutReal>& similar_to) const {
  if (is_section) {
    throw BoutException(_("Option {:s} has no value"), full_name);
  }

  auto result = bout::utils::visit(
      ConvertContainer<Matrix<BoutReal>>{
          fmt::format(
              _("Value for option {:s} cannot be converted to an Matrix<BoutReal>"),
              full_name),
          similar_to},
      value);

  // Mark this option as used
  value_used = true;

  printNameValueSourceLine(*this, "Matrix<BoutReal>");

  return result;
}

template <>
Tensor<BoutReal> Options::as<Tensor<BoutReal>>(const Tensor<BoutReal>& similar_to) const {
  if (is_section) {
    throw BoutException(_("Option {:s} has no value"), full_name);
  }

  auto result = bout::utils::visit(
      ConvertContainer<Tensor<BoutReal>>{
          fmt::format(
              _("Value for option {:s} cannot be converted to an Tensor<BoutReal>"),
              full_name),
          similar_to},
      value);

  // Mark this option as used
  value_used = true;

  printNameValueSourceLine(*this, "Tensor<BoutReal>");

  return result;
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
  Options unused = this->copy();

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
      child->second = child->second.getUnused(exclude_sources);
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
  const Options unused = getUnused();

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

std::map<std::string, const Options*> Options::subsections() const {
  std::map<std::string, const Options*> sections;
  for (const auto& child : children) {
    if (child.second.is_section) {
      sections[child.first] = &child.second;
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

  const auto* closing_brace = std::find(ctx.begin(), ctx.end(), '}');
  std::for_each(ctx.begin(), closing_brace, [&](auto ctx_opt) {
    switch (ctx_opt) {
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
    case 'u':
      unused = true;
      break;
    default:
      throw fmt::format_error("invalid format for 'Options'");
    }
  });

  // Keep a copy of the format string (without the last '}') so we can
  // pass it down to the subsections.
  const auto size = std::distance(ctx.begin(), closing_brace);
  format_string.reserve(size + 3);
  format_string.assign("{:");
  format_string.append(ctx.begin(), closing_brace);
  format_string.push_back('}');

  return closing_brace;
}

fmt::format_context::iterator
bout::details::OptionsFormatterBase::format(const Options& options,
                                            fmt::format_context& ctx) {

  const auto conditionally_used = [](const Options& option) -> bool {
    if (not option.hasAttribute(conditionally_used_attribute)) {
      return false;
    }
    return option.attributes.at(conditionally_used_attribute).as<bool>();
  };

  if (options.isValue()) {
    const std::string section_name = options.str();
    const std::string name = (inline_section_names and not section_name.empty())
                                 ? section_name
                                 : options.name();
    fmt::format_to(ctx.out(), "{}", name);

    if (not key_only) {
      const auto value = bout::utils::variantToString(options.value);
      // Convert empty strings to ""
      const std::string as_str = value.empty() ? "\"\"" : value;
      fmt::format_to(ctx.out(), " = {}", as_str);
    }

    const bool has_doc = options.attributes.count("doc") != 0U;
    const bool has_source = options.attributes.count("source") != 0U;
    const bool has_type = options.attributes.count("type") != 0U;

    std::vector<std::string> comments;

    if (unused and not options.valueUsed()) {
      if (conditionally_used(options)) {
        comments.emplace_back("unused value (marked conditionally used)");
      } else {
        comments.emplace_back("unused value (NOT marked conditionally used)");
      }
    }

    if (docstrings) {
      if (has_type) {
        comments.emplace_back(
            fmt::format("type: {}", options.attributes.at("type").as<std::string>()));
      }

      if (has_doc) {
        comments.emplace_back(
            fmt::format("doc: {}", options.attributes.at("doc").as<std::string>()));
      }
    }

    if (source and has_source) {
      const auto source = options.attributes.at("source").as<std::string>();
      if (not source.empty()) {
        comments.emplace_back(fmt::format("source: {}", source));
      }
    }

    if (not comments.empty()) {
      fmt::format_to(ctx.out(), "\t\t# {}", fmt::join(comments, ", "));
    }
    return ctx.out();
  }

  // Only print section headers if the section has a name and it has
  // non-section children
  const auto& children = options.getChildren();
  const bool has_child_values =
      std::any_of(children.begin(), children.end(),
                  [](const auto& child) { return child.second.isValue(); });
  const std::string section_name = options.str();
  if (not inline_section_names and not section_name.empty() and has_child_values) {
    fmt::format_to(ctx.out(), "\n[{}]\n", section_name);
  }

  // Get all the child values first
  for (const auto& child : children) {
    if (child.second.isValue()) {
      fmt::format_to(ctx.out(), format_string, child.second);
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
  const Options unused = options.getUnused();
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
        possible_misspellings += fmt::format("\t{:idk}\n", match.match);
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
    const std::string unused_message = _(R"(
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

{})");

    throw BoutException(unused_message, unused, data_dir, option_file,
                        unused.getChildren().begin()->first, additional_info);
  }
}
} // namespace bout
