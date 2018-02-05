typedef double BoutReal;// Do the actual multiplication of Field3D and Field3D
void autogen_Field3D_Field3D_Field3D_multiplication(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs[i] * rhs[i];
       }
}



// Provide the C function to update Field3D by multiplication with Field3D
void autogen_Field3D_Field3D_multiplication(
  BoutReal * restrict lhs, const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
        {
        lhs[i] *= rhs[i];
      }
}
// Do the actual division of Field3D and Field3D
void autogen_Field3D_Field3D_Field3D_division(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs[i] / rhs[i];
       }
}



// Provide the C function to update Field3D by division with Field3D
void autogen_Field3D_Field3D_division(
  BoutReal * restrict lhs, const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
        {
        lhs[i] /= rhs[i];
      }
}
// Do the actual addition of Field3D and Field3D
void autogen_Field3D_Field3D_Field3D_addition(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs[i] + rhs[i];
       }
}



// Provide the C function to update Field3D by addition with Field3D
void autogen_Field3D_Field3D_addition(
  BoutReal * restrict lhs, const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
        {
        lhs[i] += rhs[i];
      }
}
// Do the actual subtraction of Field3D and Field3D
void autogen_Field3D_Field3D_Field3D_subtraction(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs[i] - rhs[i];
       }
}



// Provide the C function to update Field3D by subtraction with Field3D
void autogen_Field3D_Field3D_subtraction(
  BoutReal * restrict lhs, const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
        {
        lhs[i] -= rhs[i];
      }
}
// Do the actual multiplication of Field3D and Field2D
void autogen_Field3D_Field3D_Field2D_multiplication(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal * restrict rhs, int nx,int ny,int nz) {

       for (int x=0; x < nx; ++x)
       for (int y=0; y < ny; ++y)
       for (int z=0; z < nz; ++z)
         {
         result[z + nz*(y + ny*x)] =
           lhs[z + nz*(y + ny*x)] * rhs[y + x*ny];
       }
}



// Provide the C function to update Field3D by multiplication with Field2D
void autogen_Field3D_Field2D_multiplication(
  BoutReal * restrict lhs, const BoutReal * restrict rhs, int nx,int ny,int nz) {

       for (int x=0; x < nx; ++x)
       for (int y=0; y < ny; ++y)
       for (int z=0; z < nz; ++z)
        {
        lhs[z + nz*(y + ny*x)] *= rhs[y + x*ny];
      }
}
// Do the actual division of Field3D and Field2D
void autogen_Field3D_Field3D_Field2D_division(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal * restrict rhs, int nx,int ny,int nz) {

       for (int x=0; x < nx; ++x)
       for (int y=0; y < ny; ++y)
       for (int z=0; z < nz; ++z)
         {
         result[z + nz*(y + ny*x)] =
           lhs[z + nz*(y + ny*x)] / rhs[y + x*ny];
       }
}



// Provide the C function to update Field3D by division with Field2D
void autogen_Field3D_Field2D_division(
  BoutReal * restrict lhs, const BoutReal * restrict rhs, int nx,int ny,int nz) {

       for (int x=0; x < nx; ++x)
       for (int y=0; y < ny; ++y)
       for (int z=0; z < nz; ++z)
        {
        lhs[z + nz*(y + ny*x)] /= rhs[y + x*ny];
      }
}
// Do the actual addition of Field3D and Field2D
void autogen_Field3D_Field3D_Field2D_addition(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal * restrict rhs, int nx,int ny,int nz) {

       for (int x=0; x < nx; ++x)
       for (int y=0; y < ny; ++y)
       for (int z=0; z < nz; ++z)
         {
         result[z + nz*(y + ny*x)] =
           lhs[z + nz*(y + ny*x)] + rhs[y + x*ny];
       }
}



// Provide the C function to update Field3D by addition with Field2D
void autogen_Field3D_Field2D_addition(
  BoutReal * restrict lhs, const BoutReal * restrict rhs, int nx,int ny,int nz) {

       for (int x=0; x < nx; ++x)
       for (int y=0; y < ny; ++y)
       for (int z=0; z < nz; ++z)
        {
        lhs[z + nz*(y + ny*x)] += rhs[y + x*ny];
      }
}
// Do the actual subtraction of Field3D and Field2D
void autogen_Field3D_Field3D_Field2D_subtraction(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal * restrict rhs, int nx,int ny,int nz) {

       for (int x=0; x < nx; ++x)
       for (int y=0; y < ny; ++y)
       for (int z=0; z < nz; ++z)
         {
         result[z + nz*(y + ny*x)] =
           lhs[z + nz*(y + ny*x)] - rhs[y + x*ny];
       }
}



// Provide the C function to update Field3D by subtraction with Field2D
void autogen_Field3D_Field2D_subtraction(
  BoutReal * restrict lhs, const BoutReal * restrict rhs, int nx,int ny,int nz) {

       for (int x=0; x < nx; ++x)
       for (int y=0; y < ny; ++y)
       for (int z=0; z < nz; ++z)
        {
        lhs[z + nz*(y + ny*x)] -= rhs[y + x*ny];
      }
}
// Do the actual multiplication of Field3D and BoutReal
void autogen_Field3D_Field3D_BoutReal_multiplication(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal  rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs[i] * rhs;
       }
}



// Provide the C function to update Field3D by multiplication with BoutReal
void autogen_Field3D_BoutReal_multiplication(
  BoutReal * restrict lhs, const BoutReal  rhs, int len) {

       for (int i=0; i < len; ++i)
        {
        lhs[i] *= rhs;
      }
}
// Do the actual division of Field3D and BoutReal
void autogen_Field3D_Field3D_BoutReal_division(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal  rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs[i] / rhs;
       }
}



// Provide the C function to update Field3D by division with BoutReal
void autogen_Field3D_BoutReal_division(
  BoutReal * restrict lhs, const BoutReal  rhs, int len) {

       for (int i=0; i < len; ++i)
        {
        lhs[i] /= rhs;
      }
}
// Do the actual addition of Field3D and BoutReal
void autogen_Field3D_Field3D_BoutReal_addition(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal  rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs[i] + rhs;
       }
}



// Provide the C function to update Field3D by addition with BoutReal
void autogen_Field3D_BoutReal_addition(
  BoutReal * restrict lhs, const BoutReal  rhs, int len) {

       for (int i=0; i < len; ++i)
        {
        lhs[i] += rhs;
      }
}
// Do the actual subtraction of Field3D and BoutReal
void autogen_Field3D_Field3D_BoutReal_subtraction(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal  rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs[i] - rhs;
       }
}



// Provide the C function to update Field3D by subtraction with BoutReal
void autogen_Field3D_BoutReal_subtraction(
  BoutReal * restrict lhs, const BoutReal  rhs, int len) {

       for (int i=0; i < len; ++i)
        {
        lhs[i] -= rhs;
      }
}
// Do the actual multiplication of Field2D and Field3D
void autogen_Field3D_Field2D_Field3D_multiplication(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal * restrict rhs, int nx,int ny,int nz) {

       for (int x=0; x < nx; ++x)
       for (int y=0; y < ny; ++y)
       for (int z=0; z < nz; ++z)
         {
         result[z + nz*(y + ny*x)] =
           lhs[y + x*ny] * rhs[z + nz*(y + ny*x)];
       }
}



// Do the actual division of Field2D and Field3D
void autogen_Field3D_Field2D_Field3D_division(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal * restrict rhs, int nx,int ny,int nz) {

       for (int x=0; x < nx; ++x)
       for (int y=0; y < ny; ++y)
       for (int z=0; z < nz; ++z)
         {
         result[z + nz*(y + ny*x)] =
           lhs[y + x*ny] / rhs[z + nz*(y + ny*x)];
       }
}



// Do the actual addition of Field2D and Field3D
void autogen_Field3D_Field2D_Field3D_addition(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal * restrict rhs, int nx,int ny,int nz) {

       for (int x=0; x < nx; ++x)
       for (int y=0; y < ny; ++y)
       for (int z=0; z < nz; ++z)
         {
         result[z + nz*(y + ny*x)] =
           lhs[y + x*ny] + rhs[z + nz*(y + ny*x)];
       }
}



// Do the actual subtraction of Field2D and Field3D
void autogen_Field3D_Field2D_Field3D_subtraction(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal * restrict rhs, int nx,int ny,int nz) {

       for (int x=0; x < nx; ++x)
       for (int y=0; y < ny; ++y)
       for (int z=0; z < nz; ++z)
         {
         result[z + nz*(y + ny*x)] =
           lhs[y + x*ny] - rhs[z + nz*(y + ny*x)];
       }
}



// Do the actual multiplication of Field2D and Field2D
void autogen_Field2D_Field2D_Field2D_multiplication(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs[i] * rhs[i];
       }
}



// Provide the C function to update Field2D by multiplication with Field2D
void autogen_Field2D_Field2D_multiplication(
  BoutReal * restrict lhs, const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
        {
        lhs[i] *= rhs[i];
      }
}
// Do the actual division of Field2D and Field2D
void autogen_Field2D_Field2D_Field2D_division(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs[i] / rhs[i];
       }
}



// Provide the C function to update Field2D by division with Field2D
void autogen_Field2D_Field2D_division(
  BoutReal * restrict lhs, const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
        {
        lhs[i] /= rhs[i];
      }
}
// Do the actual addition of Field2D and Field2D
void autogen_Field2D_Field2D_Field2D_addition(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs[i] + rhs[i];
       }
}



// Provide the C function to update Field2D by addition with Field2D
void autogen_Field2D_Field2D_addition(
  BoutReal * restrict lhs, const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
        {
        lhs[i] += rhs[i];
      }
}
// Do the actual subtraction of Field2D and Field2D
void autogen_Field2D_Field2D_Field2D_subtraction(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs[i] - rhs[i];
       }
}



// Provide the C function to update Field2D by subtraction with Field2D
void autogen_Field2D_Field2D_subtraction(
  BoutReal * restrict lhs, const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
        {
        lhs[i] -= rhs[i];
      }
}
// Do the actual multiplication of Field2D and BoutReal
void autogen_Field2D_Field2D_BoutReal_multiplication(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal  rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs[i] * rhs;
       }
}



// Provide the C function to update Field2D by multiplication with BoutReal
void autogen_Field2D_BoutReal_multiplication(
  BoutReal * restrict lhs, const BoutReal  rhs, int len) {

       for (int i=0; i < len; ++i)
        {
        lhs[i] *= rhs;
      }
}
// Do the actual division of Field2D and BoutReal
void autogen_Field2D_Field2D_BoutReal_division(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal  rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs[i] / rhs;
       }
}



// Provide the C function to update Field2D by division with BoutReal
void autogen_Field2D_BoutReal_division(
  BoutReal * restrict lhs, const BoutReal  rhs, int len) {

       for (int i=0; i < len; ++i)
        {
        lhs[i] /= rhs;
      }
}
// Do the actual addition of Field2D and BoutReal
void autogen_Field2D_Field2D_BoutReal_addition(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal  rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs[i] + rhs;
       }
}



// Provide the C function to update Field2D by addition with BoutReal
void autogen_Field2D_BoutReal_addition(
  BoutReal * restrict lhs, const BoutReal  rhs, int len) {

       for (int i=0; i < len; ++i)
        {
        lhs[i] += rhs;
      }
}
// Do the actual subtraction of Field2D and BoutReal
void autogen_Field2D_Field2D_BoutReal_subtraction(
  BoutReal * restrict result, const BoutReal * restrict lhs,
  const BoutReal  rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs[i] - rhs;
       }
}



// Provide the C function to update Field2D by subtraction with BoutReal
void autogen_Field2D_BoutReal_subtraction(
  BoutReal * restrict lhs, const BoutReal  rhs, int len) {

       for (int i=0; i < len; ++i)
        {
        lhs[i] -= rhs;
      }
}
// Do the actual multiplication of BoutReal and Field3D
void autogen_Field3D_BoutReal_Field3D_multiplication(
  BoutReal * restrict result, const BoutReal  lhs,
  const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs * rhs[i];
       }
}



// Do the actual division of BoutReal and Field3D
void autogen_Field3D_BoutReal_Field3D_division(
  BoutReal * restrict result, const BoutReal  lhs,
  const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs / rhs[i];
       }
}



// Do the actual addition of BoutReal and Field3D
void autogen_Field3D_BoutReal_Field3D_addition(
  BoutReal * restrict result, const BoutReal  lhs,
  const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs + rhs[i];
       }
}



// Do the actual subtraction of BoutReal and Field3D
void autogen_Field3D_BoutReal_Field3D_subtraction(
  BoutReal * restrict result, const BoutReal  lhs,
  const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs - rhs[i];
       }
}



// Do the actual multiplication of BoutReal and Field2D
void autogen_Field2D_BoutReal_Field2D_multiplication(
  BoutReal * restrict result, const BoutReal  lhs,
  const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs * rhs[i];
       }
}



// Do the actual division of BoutReal and Field2D
void autogen_Field2D_BoutReal_Field2D_division(
  BoutReal * restrict result, const BoutReal  lhs,
  const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs / rhs[i];
       }
}



// Do the actual addition of BoutReal and Field2D
void autogen_Field2D_BoutReal_Field2D_addition(
  BoutReal * restrict result, const BoutReal  lhs,
  const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs + rhs[i];
       }
}



// Do the actual subtraction of BoutReal and Field2D
void autogen_Field2D_BoutReal_Field2D_subtraction(
  BoutReal * restrict result, const BoutReal  lhs,
  const BoutReal * restrict rhs, int len) {

       for (int i=0; i < len; ++i)
         {
         result[i] =
           lhs - rhs[i];
       }
}



