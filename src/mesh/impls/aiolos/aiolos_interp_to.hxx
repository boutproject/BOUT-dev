
// to be included from aiolosmesh.hxx

void interp_to_CtoL_Field3D_x(BoutReal *__restrict__ result_ptr,
                              const BoutReal *__restrict__ in_ptr) const;
void interp_to_CtoL_Field3D_y(BoutReal *__restrict__ result_ptr,
                              const BoutReal *__restrict__ in_ptr) const;
void interp_to_CtoL_Field3D_z(BoutReal *__restrict__ result_ptr,
                              const BoutReal *__restrict__ in_ptr) const;
void interp_to_LtoC_Field3D_x(BoutReal *__restrict__ result_ptr,
                              const BoutReal *__restrict__ in_ptr) const;
void interp_to_LtoC_Field3D_y(BoutReal *__restrict__ result_ptr,
                              const BoutReal *__restrict__ in_ptr) const;
void interp_to_LtoC_Field3D_z(BoutReal *__restrict__ result_ptr,
                              const BoutReal *__restrict__ in_ptr) const;
