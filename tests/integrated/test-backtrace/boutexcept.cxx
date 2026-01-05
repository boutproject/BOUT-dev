#include "bout/boutexception.hxx"
#include "bout/msg_stack.hxx"

void troublemaker() { throw BoutException("test"); }
void f() { troublemaker(); }
void e() { f(); }

int main() {
  e();
  return 0;
}
