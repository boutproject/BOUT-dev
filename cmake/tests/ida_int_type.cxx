#include <ida/ida_bbdpre.h>
extern int ida_bbd_rhs(IDAINT, double, N_Vector, N_Vector, N_Vector, void*);
int main() { IDABBDPrecInit(nullptr, 0, 0, 0, 0, 0, 0, ida_bbd_rhs, nullptr); }
