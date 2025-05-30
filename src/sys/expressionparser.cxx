/**************************************************************************
 * Parses strings containing expressions, returning a tree of generators
 *
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#include <bout/sys/expressionparser.hxx>

#include <utility>

#include "bout/sys/gettext.hxx"
#include "bout/utils.hxx"

using std::list;
using std::string;
using std::stringstream;

using namespace std::string_literals;

using bout::generator::Context;

// Note: Here rather than in header to avoid many deprecated warnings
// Remove in future and make this function pure virtual
double FieldGenerator::generate(const Context& ctx) {
  return generate(ctx.x(), ctx.y(), ctx.z(), ctx.t());
}

/////////////////////////////////////////////
namespace { // These classes only visible in this file

// Basic generators: Numerical value, 'x', 'y' and 'z'

class FieldX : public FieldGenerator {
public:
  FieldGeneratorPtr clone(const list<FieldGeneratorPtr> UNUSED(args)) override {
    return std::make_shared<FieldX>();
  }
  double generate(const Context& ctx) override { return ctx.x(); }
  std::string str() const override { return "x"s; }
};

class FieldY : public FieldGenerator {
public:
  FieldGeneratorPtr clone(const list<FieldGeneratorPtr> UNUSED(args)) override {
    return std::make_shared<FieldY>();
  }
  double generate(const Context& ctx) override { return ctx.y(); }
  std::string str() const override { return "y"s; }
};

class FieldZ : public FieldGenerator {
public:
  FieldGeneratorPtr clone(const list<FieldGeneratorPtr> UNUSED(args)) override {
    return std::make_shared<FieldZ>();
  }
  double generate(const Context& ctx) override { return ctx.z(); }
  std::string str() const override { return "z"; }
};

class FieldT : public FieldGenerator {
public:
  FieldGeneratorPtr clone(const list<FieldGeneratorPtr> UNUSED(args)) override {
    return std::make_shared<FieldT>();
  }
  double generate(const Context& ctx) override { return ctx.t(); }
  std::string str() const override { return "t"s; }
};

class FieldParam : public FieldGenerator {
public:
  FieldParam(const std::string name) : name(name) {}
  double generate(const Context& ctx) override {
    return ctx.get(name); // Get a parameter
  }
  std::string str() const override { return "{"s + name + "}"s; }

private:
  std::string name; // The name of the parameter to look up
};

/// Define a new context to evaluate an expression in
/// This is essentially a dynamic scope mechanism
class FieldContext : public FieldGenerator {
public:
  using variable_list = std::vector<std::pair<string, FieldGeneratorPtr>>;

  /// Create with a list of context variables to modify
  /// and an expression to evaluate in that new context
  FieldContext(variable_list variables, FieldGeneratorPtr expr)
      : variables(std::move(variables)), expr(std::move(expr)) {}

  double generate(const Context& ctx) override {
    // Create a new context
    Context new_context{ctx};

    // Set values in the context by evaluating the generators
    for (auto const& var : variables) {
      new_context.set(var.first, var.second->generate(ctx));
    }

    // Evaluate the expression in the new context
    return expr->generate(new_context);
  }
  std::string str() const override {
    auto result = "["s;
    for (auto const& var : variables) {
      result += var.first + "="s + var.second->str() + ","s;
    }
    result += "]("s + expr->str() + ")"s;
    return result;
  }

private:
  variable_list variables; ///< A list of context variables to modify
  FieldGeneratorPtr expr;  ///< The expression to evaluate in the new context
};

/// Sum expressions in a loop, given symbol, count and expression
class FieldSum : public FieldGenerator {
public:
  /// Loop a symbol counter SYM from 0 to count-1
  /// The count is calculated by evaluating COUNTEXPR, which must be a non-negative integer
  /// Each iteration the expression EXPR is evaluated and the results summed.
  FieldSum(const std::string& sym, FieldGeneratorPtr countexpr, FieldGeneratorPtr expr)
      : sym(sym), countexpr(countexpr), expr(expr) {}

  double generate(const Context& ctx) override {
    // Get the count by evaluating the count expression
    BoutReal countval = countexpr->generate(ctx);
    int count = ROUND(countval);

    // Check that the count is a non-negaitve integer
    if (fabs(countval - static_cast<BoutReal>(count)) > 1e-4) {
      throw BoutException("Count {:e} is not an integer in sum expression", countval);
    }
    if (count < 0) {
      throw BoutException("Negative count {:d} in sum expression", count);
    }

    BoutReal result{0.0};
    Context new_context{ctx}; // Make a copy, so the counter value can be set
    for (int i = 0; i < count; i++) {
      // Evaluate the expression, setting the given symbol to the loop counter
      new_context.set(sym, i);
      result += expr->generate(new_context);
    }
    return result;
  }

  std::string str() const override {
    return "sum("s + sym + ","s + countexpr->str() + ","s + expr->str() + ")"s;
  }

private:
  std::string sym;
  FieldGeneratorPtr countexpr, expr;
};

} // namespace

FieldGeneratorPtr FieldBinary::clone(const list<FieldGeneratorPtr> args) {
  if (args.size() != 2) {
    throw ParseException("Binary operator expecting 2 arguments. Got {{}}", args.size());
  }

  return std::make_shared<FieldBinary>(args.front(), args.back(), op);
}

/// Convert a real value to a Boolean
/// Throw exception if `rval` isn't close to 0 or 1
bool toBool(BoutReal rval) {
  int ival = ROUND(rval);
  if ((fabs(rval - static_cast<BoutReal>(ival)) > 1e-3) or (ival < 0) or (ival > 1)) {
    throw BoutException(_("Boolean operator argument {:e} is not a bool"), rval);
  }
  return ival == 1;
}

BoutReal FieldBinary::generate(const Context& ctx) {
  BoutReal lval = lhs->generate(ctx);
  BoutReal rval = rhs->generate(ctx);

  switch (op) {
  case '|': // Logical OR
    return (toBool(lval) or toBool(rval)) ? 1.0 : 0.0;
  case '&': // Logical AND
    return (toBool(lval) and toBool(rval)) ? 1.0 : 0.0;
  case '>': // Comparison
    return (lval > rval) ? 1.0 : 0.0;
  case '<':
    return (lval < rval) ? 1.0 : 0.0;
  case '+':
    return lval + rval;
  case '-':
    return lval - rval;
  case '*':
    return lval * rval;
  case '/':
    return lval / rval;
  case '^':
    return pow(lval, rval);
  }
  // Unknown operator.
  throw ParseException("Unknown binary operator '{:c}'", op);
}

class LogicalNot : public FieldGenerator {
public:
  /// Logically negate a boolean expression
  LogicalNot(FieldGeneratorPtr expr) : expr(expr) {}

  /// Evaluate expression, check it's a bool, and return 1 or 0
  double generate(const Context& ctx) override {
    return toBool(expr->generate(ctx)) ? 0.0 : 1.0;
  }

  std::string str() const override { return "!"s + expr->str(); }

private:
  FieldGeneratorPtr expr;
};

/////////////////////////////////////////////

ExpressionParser::ExpressionParser() {
  // Add standard binary operations
  addBinaryOp('|', std::make_shared<FieldBinary>(nullptr, nullptr, '|'), 3);
  addBinaryOp('&', std::make_shared<FieldBinary>(nullptr, nullptr, '&'), 5);
  addBinaryOp('<', std::make_shared<FieldBinary>(nullptr, nullptr, '<'), 7);
  addBinaryOp('>', std::make_shared<FieldBinary>(nullptr, nullptr, '>'), 7);
  addBinaryOp('+', std::make_shared<FieldBinary>(nullptr, nullptr, '+'), 10);
  addBinaryOp('-', std::make_shared<FieldBinary>(nullptr, nullptr, '-'), 10);
  addBinaryOp('*', std::make_shared<FieldBinary>(nullptr, nullptr, '*'), 20);
  addBinaryOp('/', std::make_shared<FieldBinary>(nullptr, nullptr, '/'), 20);
  addBinaryOp('^', std::make_shared<FieldBinary>(nullptr, nullptr, '^'), 30);

  // Add standard generators
  addGenerator("x", std::make_shared<FieldX>());
  addGenerator("y", std::make_shared<FieldY>());
  addGenerator("z", std::make_shared<FieldZ>());
  addGenerator("t", std::make_shared<FieldT>());
}

void ExpressionParser::addGenerator(const string& name, FieldGeneratorPtr g) {
  gen[name] = std::move(g);
}

void ExpressionParser::addBinaryOp(char sym, FieldGeneratorPtr b, int precedence) {
  bin_op[sym] = std::make_pair(b, precedence);
  // Add to string of reserved characters
  reserved_chars += sym;
}

FieldGeneratorPtr ExpressionParser::parseString(const string& input) const {
  // Allocate a new lexer
  LexInfo lex(input, reserved_chars);

  // Parse
  auto expr = parseExpression(lex);

  // Check for remaining characters
  if (lex.curtok != 0) {
    throw ParseException("Tokens remaining unparsed in '{:s}'", input);
  }

  return expr;
}

//////////////////////////////////////////////////////////
// Private functions

std::multiset<ExpressionParser::FuzzyMatch>
ExpressionParser::fuzzyFind(const std::string& name,
                            std::string::size_type max_distance) const {
  std::multiset<ExpressionParser::FuzzyMatch> matches;
  for (const auto& key : gen) {
    if ((key.first != name) and (lowercase(key.first) == lowercase(name))) {
      // Differs only in case: pretty good match
      matches.insert({key.first, 1});
      continue;
    }
    const auto fuzzy_distance = editDistance(key.first, name);
    if (fuzzy_distance <= max_distance) {
      matches.insert({key.first, fuzzy_distance});
    }
  }
  return matches;
}

FieldGeneratorPtr ExpressionParser::parseIdentifierExpr(LexInfo& lex) const {
  // Make a nice error message if we couldn't find the identifier
  const auto generatorNotFoundErrorMessage = [&](const std::string& name) -> std::string {
    const std::string message_template = _(
        R"(Couldn't find generator '{}'. BOUT++ expressions are now case-sensitive, so you
may need to change your input file.
{})");

    // Start position of the current identifier: by this point, we've either
    // moved one character past the token, or we're still at the start
    const auto start =
        std::max(std::stringstream::off_type{0},
                 lex.ss.tellg() - std::stringstream::pos_type(lex.curident.length() + 1));

    // Try to put a nice little '^~~~~' underlining the problem
    // identifier. NOTE: this will be longer than required for multibyte UTF-8
    // characters, but will still point to the start of the identifier
    const std::string problem_bit =
        fmt::format("  {0}\n  {1: >{2}}{3:~^{4}}", lex.ss.str(), "^", start, "",
                    lex.curident.length() - 1);

    auto possible_matches = fuzzyFind(name);
    // Remove exact matches -- not helpful
    bout::utils::erase_if(possible_matches,
                          [](const auto& match) -> bool { return match.distance == 0; });

    // No matches, just point out the error
    if (possible_matches.empty()) {
      return fmt::format(message_template, name, problem_bit);
    }

    // Give the first suggestion as a possible alternative
    std::string error_message = fmt::format(message_template, name, problem_bit);
    error_message += fmt::format(_("\n  {1: ^{2}}{0}\n  Did you mean '{0}'?"),
                                 possible_matches.begin()->name, "", start);
    return error_message;
  };

  string name = lex.curident;
  lex.nextToken();

  // sum(symbol, count, expr)
  // e.g. sum(i, 10, {i}^2) -> 0 + 1^2 + 2^2 + 3^2 + ... + 9^2
  if (name == "sum") {
    if (lex.curtok != '(') {
      throw ParseException("Expecting '(' after 'sum' in 'sum(symbol, count, expr)'");
    }
    lex.nextToken();

    if ((lex.curtok != -2) && (lex.curtok != -3)) {
      throw ParseException("Expecting symbol in 'sum(symbol, count, expr)'");
    }
    string sym = lex.curident;
    lex.nextToken();

    if (lex.curtok != ',') {
      throw ParseException("Expecting , after symbol {:s} in 'sum(symbol, count, expr)'",
                           sym);
    }
    lex.nextToken();

    auto countexpr = parseExpression(lex);

    if (lex.curtok != ',') {
      throw ParseException(
          "Expecting , after count expression in 'sum(symbol, count, expr)'");
    }
    lex.nextToken();

    auto expr = parseExpression(lex);

    if (lex.curtok != ')') {
      throw ParseException("Expecting ) after expr in 'sum(symbol, count, expr)'");
    }
    lex.nextToken();
    return std::make_shared<FieldSum>(sym, countexpr, expr);
  }

  if (lex.curtok == '(') {
    // Argument list. Find if a generator or function

    auto it = gen.find(name);
    if (it == gen.end()) {
      throw ParseException(generatorNotFoundErrorMessage(name));
    }

    // Parse arguments (if any)
    list<FieldGeneratorPtr> args;

    lex.nextToken();
    if (lex.curtok == ')') {
      // Empty list
      lex.nextToken();
      return it->second->clone(args);
    }
    do {
      // Should be an expression
      args.push_back(parseExpression(lex));

      // Now either a comma or ')'

      if (lex.curtok == ')') {
        // Finished list
        lex.nextToken();
        return it->second->clone(args);
      }
      if (lex.curtok != ',') {
        throw ParseException("Expecting ',' or ')' in function argument list ({:s})\n",
                             name);
      }
      lex.nextToken();
    } while (true);

  } else {
    // No arguments. Search in generator list
    auto it = gen.find(name);
    if (it == gen.end()) {
      // Not in internal map. Try to resolve
      FieldGeneratorPtr g = resolve(name);
      if (g == nullptr) {
        throw ParseException(generatorNotFoundErrorMessage(name));
      }
      return g;
    }
    list<FieldGeneratorPtr> args;
    return it->second->clone(args);
  }
}

FieldGeneratorPtr ExpressionParser::parseParenExpr(LexInfo& lex) const {
  lex.nextToken(); // eat '('

  FieldGeneratorPtr g = parseExpression(lex);

  if ((lex.curtok != ')') && (lex.curtok != ']')) {
    throw ParseException("Expecting ')' or ']' but got curtok={:d} ({:c})",
                         static_cast<int>(lex.curtok), static_cast<char>(lex.curtok));
  }

  lex.nextToken(); // eat ')'
  return g;
}

// This will return a pointer to a FieldContext.
FieldGeneratorPtr ExpressionParser::parseContextExpr(LexInfo& lex) const {
  lex.nextToken(); // eat '['

  FieldContext::variable_list variables;

  while (lex.curtok != ']') {
    if (lex.curtok == 0) {
      throw ParseException("Expecting ']' in context expression");
    }

    // Definition, ident = expression
    // First comes the identifier symbol
    if (lex.curtok != -2) {
      throw ParseException(
          "Expecting an identifier in context expression, but got curtok={:d} ({:c})",
          static_cast<int>(lex.curtok), lex.curtok);
    }
    string symbol = lex.curident;
    lex.nextToken();

    // Now should be '='
    if (lex.curtok != '=') {
      throw ParseException(
          "Expecting '=' after '{:s}' in context expression, but got curtok={:d} ({:c})",
          symbol, static_cast<int>(lex.curtok), lex.curtok);
    }
    lex.nextToken();

    // Should be expression
    FieldGeneratorPtr value = parseExpression(lex);

    variables.push_back(std::make_pair(symbol, value));

    if (lex.curtok == ',') {
      lex.nextToken(); // Skip comma
    }
  }
  lex.nextToken(); // eat ']'

  // Should now be '('
  if (lex.curtok != '(') {
    throw ParseException(
        "Expecting '(' after ] context expression,  but got curtok={:d} ({:c})",
        static_cast<int>(lex.curtok), static_cast<char>(lex.curtok));
  }

  // Get the next expression to evaluate, put into FieldContext
  // Note: Ensure that only the first expression in parentheses is parsed
  //       by calling parseParenExpr rather than parseExpression
  return std::make_shared<FieldContext>(variables, parseParenExpr(lex));
}

FieldGeneratorPtr ExpressionParser::parsePrimary(LexInfo& lex) const {
  switch (lex.curtok) {
  case -1: {         // a number
    lex.nextToken(); // Eat number
    return std::make_shared<FieldValue>(lex.curval);
  }
  case -2: {
    return parseIdentifierExpr(lex);
  }
  case -3: {
    // A parameter, passed as an argument to generate
    auto gen = std::make_shared<FieldParam>(lex.curident);
    lex.nextToken();
    return gen;
  }
  case '-': {
    // Unary minus
    // Don't eat the minus, and return an implicit zero
    return std::make_shared<FieldValue>(0.0);
  }
  case '!': {
    // Logical not
    lex.nextToken(); // Eat '!'
    return std::make_shared<LogicalNot>(parsePrimary(lex));
  }
  case '(': {
    return parseParenExpr(lex);
  }
  case '[': {
    // Define a new context (scope).
    return parseContextExpr(lex);
  }
  }
  throw ParseException("Unexpected token {:d} ({:c})", static_cast<int>(lex.curtok),
                       static_cast<char>(lex.curtok));
}

FieldGeneratorPtr ExpressionParser::parseBinOpRHS(LexInfo& lex, int ExprPrec,
                                                  FieldGeneratorPtr lhs) const {

  while (true) {
    // Check for end of input
    if ((lex.curtok == 0) || (lex.curtok == ')') || (lex.curtok == ',')
        || (lex.curtok == ']')) {
      return lhs;
    }

    // Next token should be a binary operator
    auto it = bin_op.find(lex.curtok);

    if (it == bin_op.end()) {
      throw ParseException("Unexpected binary operator '{:c}'",
                           static_cast<char>(lex.curtok));
    }

    FieldGeneratorPtr op = it->second.first;
    int TokPrec = it->second.second;

    if (TokPrec < ExprPrec) {
      return lhs;
    }

    lex.nextToken(); // Eat binop

    FieldGeneratorPtr rhs = parsePrimary(lex);

    if ((lex.curtok == 0) || (lex.curtok == ')') || (lex.curtok == ',')
        || (lex.curtok == ']')) {
      // Done

      list<FieldGeneratorPtr> args;
      args.push_front(lhs);
      args.push_back(rhs);
      return op->clone(args);
    }

    // Find next binop
    it = bin_op.find(lex.curtok);

    if (it == bin_op.end()) {
      throw ParseException("Unexpected character '{:c}' ({:d})",
                           static_cast<char>(lex.curtok), static_cast<int>(lex.curtok));
    }

    int NextPrec = it->second.second;
    if (TokPrec < NextPrec) {
      rhs = parseBinOpRHS(lex, TokPrec + 1, rhs);
    }

    // Merge lhs and rhs into new lhs
    list<FieldGeneratorPtr> args;
    args.push_front(lhs);
    args.push_back(rhs);
    lhs = op->clone(args);
  }
}

FieldGeneratorPtr ExpressionParser::parseExpression(LexInfo& lex) const {
  FieldGeneratorPtr lhs = parsePrimary(lex);
  return parseBinOpRHS(lex, 0, lhs);
}

//////////////////////////////////////////////////////////
// LexInfo

ExpressionParser::LexInfo::LexInfo(const std::string& input, std::string reserved_chars)
    : reserved_chars(std::move(reserved_chars)) {
  ss.clear();
  ss.str(input); // Set the input stream
  ss.seekg(0, std::ios_base::beg);

  LastChar = static_cast<signed char>(ss.get()); // First char from stream
  nextToken();                                   // Get first token
}

char ExpressionParser::LexInfo::nextToken() {
  while (isspace(static_cast<unsigned char>(LastChar))) {
    LastChar = static_cast<signed char>(ss.get());
  }

  if (!ss.good()) {
    curtok = 0;
    return 0;
  }

  // Handle numbers
  if (isdigit(static_cast<unsigned char>(LastChar))
      || (LastChar == '.')) { // Number: [0-9.]+
    bool gotdecimal = false, gotexponent = false;
    std::string NumStr;

    while (true) {
      if (LastChar == '.') {
        if (gotdecimal || gotexponent) {
          throw ParseException("Unexpected '.' in number expression");
        }
        gotdecimal = true;
      } else if ((LastChar == 'E') || (LastChar == 'e')) {
        if (gotexponent) {
          throw ParseException(
              "ExpressionParser error: Unexpected extra 'e' in number expression");
        }
        gotexponent = true;
        // Next character should be a '+' or '-' or digit
        NumStr += 'e';
        LastChar = static_cast<signed char>(ss.get());
        if ((LastChar != '+') && (LastChar != '-')
            && !isdigit(static_cast<unsigned char>(LastChar))) {
          throw ParseException(
              "ExpressionParser error: Expecting '+', '-' or number after 'e'");
        }
      } else if (!isdigit(static_cast<unsigned char>(LastChar))) {
        break;
      }

      NumStr += LastChar;
      LastChar = static_cast<signed char>(ss.get());
    }

    curval = std::stod(NumStr);
    curtok = -1;
    return curtok;
  }

  // Symbols can contain anything else which is not reserved
  if ((LastChar == '`') || (reserved_chars.find(LastChar) == std::string::npos)) {

    // Special case: If the last token returned was a number
    // then insert a multiplication ("*")
    if (curtok == -1) {
      curtok = '*';
      return curtok;
    }

    curident.clear();
    do {
      if (LastChar == '\\') {
        // Escape character.
        // Whatever the next character is, include it in the identifier
        // Note: Even though this only treats one character specially,
        //       it should still work for utf8 since all chars are
        //       allowed except reserved_chars and whitespace
        LastChar = static_cast<signed char>(ss.get());
        if (LastChar == EOF) {
          throw ParseException("Unexpected end of input after \\ character");
        }

        curident += LastChar;
      } else if (LastChar == '`') {
        // An escaped symbol
        // Include all characters until the next ` (backquote)
        LastChar = static_cast<signed char>(ss.get()); // Skip the `
        do {
          curident += LastChar;
          LastChar = static_cast<signed char>(ss.get());
          if (LastChar == EOF) {
            throw ParseException("Unexpected end of input; expecting ` (backquote)");
          }
        } while (LastChar != '`');
        // Final ` will not be added to the symbol
      } else {
        curident += LastChar;
      }
      LastChar = static_cast<signed char>(ss.get());
    } while ((LastChar != EOF && !isspace(static_cast<unsigned char>(LastChar))
              && (reserved_chars.find(LastChar) == std::string::npos))
             || (LastChar == '\\') || (LastChar == '`'));
    curtok = -2;
    return curtok;
  }

  if (LastChar == '(') {
    // Special case: If the last token returned was a number
    // then insert a multiplication ("*") before the opening bracket
    if (curtok == -1) {
      curtok = '*';
      return curtok;
    }
  }

  if (LastChar == '{') {
    // A special quoted name, which is turned into a FieldParam
    // and used to look up an input parameter

    // Special case: If the last token returned was a number
    // then insert a multiplication ("*") before the opening brace
    if (curtok == -1) {
      curtok = '*';
      return curtok;
    }

    curident.clear();

    LastChar = static_cast<signed char>(ss.get()); // Skip the {
    do {
      curident += LastChar;
      LastChar = static_cast<signed char>(ss.get());
      if (LastChar == EOF) {
        throw ParseException("Unexpected end of input; expecting }}");
      }
      if (LastChar == '{') {
        throw ParseException("Unexpected opening brace {{; expecting }}");
      }
    } while (LastChar != '}');
    LastChar = static_cast<signed char>(ss.get());
    curtok = -3;
    return curtok;
  }

  // LastChar is unsigned, explicitly cast
  curtok = LastChar;
  LastChar = static_cast<signed char>(ss.get());
  return curtok;
}
