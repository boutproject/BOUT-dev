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

#include <stdio.h>

#define MSG_MAX_SIZE 127

typedef struct {
  char str[MSG_MAX_SIZE+1];
}msg_item_t;

class MsgStack {
 public:
  MsgStack();
  ~MsgStack();
  
#if CHECK > 1
  int push(const char *s, ...); ///< Add a message to the stack. Returns a message id
  
  int setPoint();     ///< get a message point
  
  void pop();          ///< Remove the last message
  void pop(int id);    ///< Remove all messages back to msg <id>
  void clear();        ///< Clear all message
  
  void dump();         ///< Write out all messages (using output)
#else
  /// Dummy functions which should be optimised out
  int push(const char *s, ...) {return 0;}
  
  int setPoint() {return 0;}
  
  void pop() {}
  void pop(int id) {}
  void clear() {}
  
  void dump() {}
#endif
  
 private:
  char buffer[256];
  
  msg_item_t *msg;  ///< Message stack;
  int nmsg;    ///< Current number of messages
  int size;    ///< Size of the stack
};

#ifndef GLOBALORIGIN
#define GLOBAL extern
#else
#define GLOBAL
#endif

/// Global object. Will eventually replace with better system
GLOBAL MsgStack msg_stack;

#undef GLOBAL

#endif // __MSG_STACK_H__

