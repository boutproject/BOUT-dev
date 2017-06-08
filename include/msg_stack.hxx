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

#include "string_utils.hxx"

#include <exception>
#include <stdio.h>
#include <string>
#include <vector>

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
  MsgStack() {}
  ~MsgStack() {}

#if CHECK > 1
  /// Add a message to the stack. Returns a message id
  template <typename... Args> int push(const std::string &message, Args... args) {
    return push(string_format(message, args...));
  }
  int push(const std::string &message);

  int setPoint();        ///< get a message point

  void pop();            ///< Remove the last message
  void pop(int id);      ///< Remove all messages back to msg <id>
  void clear();          ///< Clear all message

  void dump();           ///< Write out all messages (using output)
  std::string getDump(); ///< Write out all messages to a string
#else
  /// Dummy functions which should be optimised out
  template <typename... Args> int push(const std::string &message, Args... args) {
    return 0;
  }
  int push(const std::string &message) { return 0; }

  int setPoint() { return 0; }

  void pop() {}
  void pop(int id) {}
  void clear() {}

  void dump() {}
#endif

private:
  std::vector<std::string> message_stack; ///< The stack of messages
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
  MsgStackItem(const char *msg) { point = msg_stack.push(msg); }
  MsgStackItem(const char *msg, const char *file, int line) {
    point = msg_stack.push("%s on line %d of '%s'", msg, line, file);
  }
  ~MsgStackItem() {
    // If an exception has occurred, don't pop the message
    if (!std::uncaught_exception())
      msg_stack.pop(point);
  }

private:
  int point;
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
#ifdef CHECK
#define TRACE(message) MsgStackItem CONCATENATE(msgTrace_ , __LINE__) (message, __FILE__, __LINE__)
#else
#define TRACE(message)
#endif

#endif // __MSG_STACK_H__
