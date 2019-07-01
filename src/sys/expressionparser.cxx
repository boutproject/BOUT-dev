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

#include <utils.hxx> // for lowercase

using std::list;
using std::string;
using std::stringstream;

using bout::generator::Context;

// Note: Here rather than in header to avoid many deprecated warnings
// Remove in future and make this function pure virtual
double FieldGenerator::generate(const Context& pos) {
  return generate(pos.x(), pos.y(), pos.z(), pos.t());
}


/////////////////////////////////////////////
namespace { // These classes only visible in this file

// Basic generators: Numerical value, 'x', 'y' and 'z'

class FieldX : public FieldGenerator {
public:
  FieldGeneratorPtr clone(const list<FieldGeneratorPtr> UNUSED(args)) override {
    return std::make_shared<FieldX>();
  }
  double generate(const Context& pos) override {
    return pos.x();
  }
  std::string str() const override { return std::string("x"); }
};

class FieldY : public FieldGenerator {
public:
  FieldGeneratorPtr clone(const list<FieldGeneratorPtr> UNUSED(args)) override {
    return std::make_shared<FieldY>();
  }
  double generate(const Context& pos) override {
    return pos.y();
  }
  std::string str() const override { return std::string("y"); }
};

class FieldZ : public FieldGenerator {
public:
  FieldGeneratorPtr clone(const list<FieldGeneratorPtr> UNUSED(args)) override {
    return std::make_shared<FieldZ>();
  }
  double generate(const Context& pos) override {
    return pos.z();
  }
  std::string str() const override { return std::string("z"); }
};

class FieldT : public FieldGenerator {
public:
  FieldGeneratorPtr clone(const list<FieldGeneratorPtr> UNUSED(args)) override {
    return std::make_shared<FieldT>();
  }
  double generate(const Context& pos) override {
    return pos.t();
  }
  std::string str() const override { return std::string("t"); }
};

class FieldParam : public FieldGenerator {
public:
  FieldParam(const std::string name) : name(name) {}
  double generate(const Context& pos) override {
    return pos.get(name); // Get a parameter
  }
  std::string str() const override { return std::string("{") + name + std::string("}"); }
private:
  std::string name; // The name of the parameter to look up
};
} // namespace

FieldGeneratorPtr FieldBinary::clone(const list<FieldGeneratorPtr> args) {
  if (args.size() != 2)
    throw ParseException("Binary operator expecting 2 arguments. Got '%lu'",
                         static_cast<unsigned long>(args.size()));

  return std::make_shared<FieldBinary>(args.front(), args.back(), op);
}

BoutReal FieldBinary::generate(const Context& pos) {
  BoutReal lval = lhs->generate(pos);
  BoutReal rval = rhs->generate(pos);
  switch(op) {
  case '+': return lval + rval;
  case '-': return lval - rval;
  case '*': return lval * rval;
  case '/': return lval / rval;
  case '^': return pow(lval, rval);
  }
  // Unknown operator.
  throw ParseException("Unknown binary operator '%c'", op);
}

/////////////////////////////////////////////

ExpressionParser::ExpressionParser() {
  // Add standard binary operations
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
  gen[name] = g;
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
  return parseExpression(lex);
}

//////////////////////////////////////////////////////////
// Private functions

FieldGeneratorPtr ExpressionParser::parseIdentifierExpr(LexInfo& lex) const {
  string name = lowercase(lex.curident);
  lex.nextToken();

  if (lex.curtok == '(') {
    // Argument list. Find if a generator or function

    auto it = gen.find(name);
    if (it == gen.end())
      throw ParseException("Couldn't find generator '%s'", name.c_str());

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
        throw ParseException("Expecting ',' or ')' in function argument list (%s)\n",
                             name.c_str());
      }
      lex.nextToken();
    } while (true);

  } else {
    // No arguments. Search in generator list
    auto it = gen.find(name);
    if (it == gen.end()) {
      // Not in internal map. Try to resolve
      FieldGeneratorPtr g = resolve(name);
      if (g == nullptr)
        throw ParseException("Couldn't find generator '%s'", name.c_str());
      return g;
    }
    list<FieldGeneratorPtr> args;
    return it->second->clone(args);
  }
}

FieldGeneratorPtr ExpressionParser::parseParenExpr(LexInfo& lex) const {
  lex.nextToken(); // eat '('

  FieldGeneratorPtr g = parseExpression(lex);

  if ((lex.curtok != ')') && (lex.curtok != ']'))
    throw ParseException("Expecting ')' or ']' but got curtok=%d (%c)",
                         static_cast<int>(lex.curtok), lex.curtok);

  lex.nextToken(); // eat ')'
  return g;
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
  case '(':
  case '[':
    return parseParenExpr(lex);
  }
  throw ParseException("Unexpected token %d (%c)", static_cast<int>(lex.curtok),
                       lex.curtok);
}

FieldGeneratorPtr ExpressionParser::parseBinOpRHS(LexInfo& lex, int ExprPrec,
                                                  FieldGeneratorPtr lhs) const {

  while (true) {
    // Check for end of input
    if ((lex.curtok == 0) || (lex.curtok == ')') || (lex.curtok == ','))
      return lhs;

    // Next token should be a binary operator
    auto it = bin_op.find(lex.curtok);

    if (it == bin_op.end())
      throw ParseException("Unexpected binary operator '%c'", lex.curtok);

    FieldGeneratorPtr op = it->second.first;
    int TokPrec = it->second.second;

    if (TokPrec < ExprPrec)
      return lhs;

    lex.nextToken(); // Eat binop

    FieldGeneratorPtr rhs = parsePrimary(lex);

    if ((lex.curtok == 0) || (lex.curtok == ')') || (lex.curtok == ',')) {
      // Done

      list<FieldGeneratorPtr> args;
      args.push_front(lhs);
      args.push_back(rhs);
      return op->clone(args);
    }

    // Find next binop
    it = bin_op.find(lex.curtok);

    if (it == bin_op.end())
      throw ParseException("Unexpected character '%c' (%d)", lex.curtok, static_cast<int>(lex.curtok));

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
  while (isspace(LastChar))
    LastChar = static_cast<signed char>(ss.get());

  if (!ss.good()) {
    curtok = 0;
    return 0;
  }

  // Handle numbers
  if (isdigit(LastChar) || (LastChar == '.')) { // Number: [0-9.]+
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
        if ((LastChar != '+') && (LastChar != '-') && !isdigit(LastChar)) {
          throw ParseException(
              "ExpressionParser error: Expecting '+', '-' or number after 'e'");
        }
      } else if (!isdigit(LastChar))
        break;

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
    } while ((LastChar != EOF && !isspace(LastChar)
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
        throw ParseException("Unexpected end of input; expecting }");
      }
      if (LastChar == '{') {
        throw ParseException("Unexpected opening brace {; expecting }");
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

//////////////////////////////////////////////////////////
// ParseException

ParseException::ParseException(const char* s, ...) {
  if (s == nullptr)
    return;

  int buf_len = 1024;
  char* buffer = new char[buf_len];
  bout_vsnprintf(buffer, buf_len, s);

  message.assign(buffer);
  delete[] buffer;
}

const char* ParseException::what() const noexcept { return message.c_str(); }
