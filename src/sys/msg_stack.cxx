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

#include "globals.hxx" // For the output object
#include "msg_stack.hxx"
#include <string.h>
#include <stdarg.h>

MsgStack::MsgStack()
{
  nmsg = 0;
  size = 0;
}

MsgStack::~MsgStack()
{
  clear();
}

int MsgStack::push(const char *s, ...)
{
#if CHECK > 1
  va_list ap;  // List of arguments
  msg_item_t *m;

  if(size > nmsg) {
    m = &msg[nmsg];
  }else {
    // need to allocate more memory
    if(size == 0) {
      msg = (msg_item_t*) malloc(sizeof(msg_item_t)*10);
      size = 10;
      m = msg;
    }else {
      msg = (msg_item_t*) realloc(msg, sizeof(msg_item_t)*(size + 10));
      m = &msg[size];
      size += 10;
    }
  }

  if(s != NULL) {

    va_start(ap, s);
      vsprintf(buffer, s, ap);
    va_end(ap);
    
    strncpy(m->str, buffer, MSG_MAX_SIZE);
    m->str[MSG_MAX_SIZE] = '\0';

    //output.write("Pushing '%s' -> %d\n", buffer, nmsg);
  }else
    m->str[0] = '\0';

  nmsg++;
  return nmsg-1;
#else
  return 0;
#endif
}

int MsgStack::setPoint()
{
#if CHECK > 1
  // Create an empty message

  return push(NULL);
#else
  return 0;
#endif
}

void MsgStack::pop()
{
#if CHECK > 1
  
  if(nmsg <= 0)
    return;
  
  //output.write("Popping %d\n", nmsg);

  nmsg--;
#endif
}

void MsgStack::pop(int id)
{
#if CHECK > 1
  if(id < 0)
    id = 0;

  if(id > nmsg)
    return;

  nmsg = id;

#endif
}

void MsgStack::clear()
{
#if CHECK > 1
  nmsg = 0;
#endif
}

void MsgStack::dump()
{
#if CHECK > 1
  output.write("====== Back trace ======\n");

  for(int i=nmsg-1;i>=0;i--) {
    if(msg[i].str[0] != '\0') {
      output.write(" -> ");
      output.write(msg[i].str);
      output.write("\n");
    }
  }
#endif
}
