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

#include <msg_stack.hxx>
#include <output.hxx>
#include <string.h>
#include <string>
#include <stdarg.h>

#if CHECK > 1
int MsgStack::push(const std::string &message) {
  message_stack.push_back(message);
  return message_stack.size() - 1;
}

int MsgStack::setPoint() {
  // Create an empty message
  return push("");
}

void MsgStack::pop() {
  if (message_stack.size() <= 0) {
    return;
  }

  message_stack.pop_back();
}

void MsgStack::pop(int id) {
  if (id < 0) {
    id = 0;
  }

  if (id > message_stack.size()) {
    return;
  }

  // Erase from message id to end
  message_stack.erase(std::begin(message_stack) + id, std::end(message_stack));
}

void MsgStack::clear() { message_stack.clear(); }

void MsgStack::dump() { output << this->getDump(); }

std::string MsgStack::getDump() {
  std::string result = "====== Back trace ======\n";
  // Loop backwards over message stack
  for (auto index = std::rbegin(message_stack); index != std::rend(message_stack);
       ++index) {
    result += " -> " + *index + "\n";
  }
  return result;
}

#endif
