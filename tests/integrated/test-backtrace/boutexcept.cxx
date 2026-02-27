#include "bout/boutexception.hxx"

// Can't use anonymous namespace, or the compiler will inline everything,
// defeating the point of this test.
// NOLINTBEGIN(misc-use-internal-linkage)
void troublemaker() { throw BoutException("test"); }
void f() { troublemaker(); }
void e() { f(); }
// NOLINTEND(misc-use-internal-linkage)

int main() {
  e();
  return 0;
}
