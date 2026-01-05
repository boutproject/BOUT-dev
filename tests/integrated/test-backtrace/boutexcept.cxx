#include "bout/boutexception.hxx"

void troublemaker() { throw BoutException("test"); }
void f() { troublemaker(); }
void e() { f(); }

int main() {
  e();
  return 0;
}
