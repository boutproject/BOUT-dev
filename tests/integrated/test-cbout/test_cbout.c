#include <cbout/bout.h>

int main() {
  Field3D *f = Field3D_new();

  Field3D_delete(f);
  return 0;
}
