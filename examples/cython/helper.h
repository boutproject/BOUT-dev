#include <field3d.hxx>

void c_set_f3d_all(Field3D * f3d, double * data);
void c_get_f3d_all(Field3D * f3d, double * data);
Field3D * fadd( Field3D*,Field3D*);
Field3D * fmul( Field3D*,Field3D*);
int getNx(Field3D * a);
int getNy(Field3D * a);
int getNz(Field3D * a);

Mesh * c_get_global_mesh();
