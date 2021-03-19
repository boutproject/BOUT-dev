#include "boutexception.hxx"
#include "msg_stack.hxx"

void troublemaker() {
  AUTO_TRACE();
  throw BoutException("test");
}
void f() {
  AUTO_TRACE();
  troublemaker();
}
void e() {
  AUTO_TRACE();
  f();
}

int main() {
  e();
  return 0;
}
