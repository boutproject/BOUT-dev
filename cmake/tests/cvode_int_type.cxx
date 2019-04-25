#include <cvode/cvode_bbdpre.h>
extern int cvode_bbd_rhs(CVODEINT, double, N_Vector, N_Vector, void*);
int main() { CVBBDPrecInit(nullptr, 0, 0, 0, 0, 0, 0, cvode_bbd_rhs, nullptr); }
