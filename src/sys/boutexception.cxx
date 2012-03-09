#include <boutexception.hxx>
#include <msg_stack.hxx>
#include <iostream>
#include <stdarg.h>
#include <output.hxx>

BoutException::~BoutException() throw()
{
}

BoutException::BoutException(const char* s, ...)
{
  va_list ap;  // List of arguments

  if(s == (const char*) NULL)
    return;
  
  char buffer[1024];
  va_start(ap, s);
    vsprintf(buffer, s, ap);
  va_end(ap);
  
  message.assign(buffer);
}

const char* BoutException::what() const throw()
{
  #ifdef CHECK
    /// Print out the message stack to help debugging
    msg_stack.dump();
  #else
    output.write("Enable checking (-DCHECK flag) to get a trace\n");
  #endif
  return message.c_str();
}
