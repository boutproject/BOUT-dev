#include "bout/boutexception.hxx"

namespace {
void troublemaker() { throw BoutException("test"); }
void f() { troublemaker(); }
void e() { f(); }
} // namespace

int main() {
  e();
  return 0;
}
