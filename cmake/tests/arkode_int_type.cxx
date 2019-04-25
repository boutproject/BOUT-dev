#include <arkode/arkode_bbdpre.h>
extern int arkode_bbd_rhs(ARKODEINT, double, N_Vector, N_Vector, void*);
int main() { ARKBBDPrecInit(nullptr, 0, 0, 0, 0, 0, 0, arkode_bbd_rhs, nullptr); }
