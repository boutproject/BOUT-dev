/*!************************************************************************
 * Provides a message stack to print more useful error
 * messages.
 *
 **************************************************************************
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

class MsgStack;

#ifndef __MSG_STACK_H__
#define __MSG_STACK_H__

#include "unused.hxx"

#include <exception>
#include <stdarg.h>
#include <string>
#include <vector>

/// The maximum length (in chars) of messages, not including terminating '0'
#define MSG_MAX_SIZE 127

/*!
 * Message stack
 *
 * Implements a stack of messages which can be pushed onto the top
 * and popped off the top. This is used for debugging: messages are put
 * into this stack at the start of a section of code, and removed at the end.
 * If an error occurs in between push and pop, then the message can be printed.
 *
 * This code is only enabled if CHECK > 1. If CHECK is disabled then this
 * message stack code reverts to empty functions which should be removed by
 * the optimiser
 */
class MsgStack {
public:
  MsgStack() : position(0){};
  ~MsgStack() { clear(); }

#if CHECK > 1
  int push(const char *s, ...); ///< Add a message to the stack. Returns a message id

  int setPoint(); ///< get a message point

  void pop();       ///< Remove the last message
  void pop(int id); ///< Remove all messages back to msg \p id
  void clear();     ///< Clear all message

  void dump();           ///< Write out all messages (using output)
  std::string getDump(); ///< Write out all messages to a string
#else
  /// Dummy functions which should be optimised out
  int push(const char *UNUSED(s), ...) { return 0; }

  int setPoint() { return 0; }

  void pop() {}
  void pop(int UNUSED(id)) {}
  void clear() {}

  void dump() {}
  std::string getDump() { return ""; }
#endif

private:
  char buffer[256]; ///< Buffer for vsnprintf

  std::vector<std::string> stack;               ///< Message stack;
  std::vector<std::string>::size_type position; ///< Position in stack
};

/*!
 * This is a way to define a global object,
 * so that it is declared extern in all files except one
 * where GLOBALORIGIN is defined.
 */
#ifndef GLOBALORIGIN
#define GLOBAL extern
#else
#define GLOBAL
#endif

/// Global object. Will eventually replace with better system
GLOBAL MsgStack msg_stack;

#undef GLOBAL

/*!
 * MsgStackItem
 *
 * Simple class to manage pushing and popping messages
 * from the message stack. Pushes a message in the
 * constructor, and pops the message on destruction.
 */
class MsgStackItem {
public:
  // Not currently used anywhere
  MsgStackItem(const char *msg) { point = msg_stack.push(msg); }
  // Not currently used anywhere
  MsgStackItem(const char *msg, const char *file, int line) {
    point = msg_stack.push("%s on line %d of '%s'", msg, line, file);
  }
  MsgStackItem(const char *file, int line, const char *msg, ...) {
    va_list args;
    va_start(args, msg);
    vsnprintf(buffer, MSG_MAX_SIZE, msg, args);
    point = msg_stack.push("%s on line %d of '%s'", buffer, line, file);
    va_end(args);
  }
  ~MsgStackItem() {
    // If an exception has occurred, don't pop the message
    if (!std::uncaught_exception()) {
      msg_stack.pop(point);
    }
  }

private:
  int point;
  char buffer[256];
};

/// To concatenate strings for a variable name
#define CONCATENATE_DIRECT(s1, s2) s1##s2
/// Need to use two levels due to macro strangeness
#define CONCATENATE(s1, s2) CONCATENATE_DIRECT(s1, s2)

/*!
 * The TRACE macro provides a convenient way to put messages onto the msg_stack
 * It pushes a message onto the stack, and pops it when the scope ends
 *
 * Example
 * -------
 *
 * {
 *   TRACE("Starting calculation")
 *
 * } // Scope ends, message popped
 */
#if CHECK > 0
/* Would like to have something like TRACE(message, ...) so that we can directly refer
   to the (required) first argument, which is the main message string. However because
   we want to allow TRACE("Message with no args") we have to deal with the case where
   __VA_ARGS__ is empty. There's a GCC specific extension such that
    //#define TRACE(message, ...) MsgStackItem CONCATENATE(msgTrace_ , __LINE__) (message,
   __FILE__, __LINE__, ##__VA_ARGS__) //## is non-standard here would achieve this for us.
   However to be more portable have to instead just reorder the arguments from the
   original MsgStackItem constructor so that the message is the last of the required
   arguments and the optional arguments follow from there.
 */
#define TRACE(...)                                                                       \
  MsgStackItem CONCATENATE(msgTrace_, __LINE__)(__FILE__, __LINE__, __VA_ARGS__)
#else
#define TRACE(...)
#endif

#endif // __MSG_STACK_H__
