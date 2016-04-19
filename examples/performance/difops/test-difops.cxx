/*
 * Test case for differential operators on single 
 * index, performance comparison against current operators
 * which contain individual loops over the domain.
 * 
 * The equations tested are a subset of those in the 2D blob
 * simulations (examples/blob2d)
 * 
 * New operators seem to be slightly faster:
 * 
 * Inner loops          : 0.379856
 * Outer loop           : 0.218945
 * 
 */

#include <bout.hxx>
#include <derivs.hxx>

#include <chrono>
#include <iostream>

BoutReal DDX(const Field3D &f, const DataIterator &i) {
  return (f[i.xp()] - f[i.xm()])/(2.*mesh->coordinates()->dx[i]);
}

BoutReal DDZ(const Field3D &f, const DataIterator &i) {
  return (f[i.zp()] - f[i.zm()])/(2.*mesh->coordinates()->dz);
}

BoutReal bracket_arakawa(const Field3D &f, const Field3D &g, const DataIterator &i) {
  Coordinates *metric = mesh->coordinates();
  
  // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
  BoutReal Jpp = 0.25*( (f[i.zp()] - f[i.zm()])*
                        (g[i.xp()] - g[i.xm()]) -
                        (f[i.xp()] - f[i.xm()])*
                        (g[i.zp()] - g[i.zm()]) );
  
  // J+x
  BoutReal Jpx = 0.25*( g[i.xp()]*(f[i.offset(1,0,1)]-f[i.offset(1,0,-1)]) -
                        g[i.xm()]*(f[i.offset(-1,0,1)]-f[i.offset(-1,0,-1)]) -
                        g[i.zp()]*(f[i.offset(1,0,1)]-f[i.offset(-1,0,1)]) +
                        g[i.zm()]*(f[i.offset(1,0,-1)]-f[i.offset(-1,0,-1)]));
  // Jx+
  BoutReal Jxp = 0.25*( g[i.offset(1,0,1)]*(f[i.zp()]-f[i.xp()]) -
                        g[i.offset(-1,0,-1)]*(f[i.xm()]-f[i.zm()]) -
                        g[i.offset(-1,0,1)]*(f[i.zp()]-f[i.xm()]) +
                        g[i.offset(1,0,-1)]*(f[i.xp()]-f[i.zm()]) );
  
  return (Jpp + Jpx + Jxp) / (3. * metric->dx[i] * metric->dz);
}

int main(int argc, char **argv) {
  BoutInitialise(argc, argv);

  typedef std::chrono::time_point<std::chrono::steady_clock> SteadyClock;
  typedef std::chrono::duration<double> Duration;
  using namespace std::chrono;
  
  Field3D phi, n, vort;
  GRID_LOAD3(phi, n, vort);
  
  BoutReal R_c, L_par;
  GRID_LOAD2(R_c, L_par);
  
  /////////////////////////////////////////////////
  // Single loop over domain

  SteadyClock start2 = steady_clock::now();
  // Note: Need to allocate ddt() fields before [] access
  ddt(n).allocate();
  ddt(vort).allocate();
  for(auto i : n.region(RGN_NOBNDRY)) {
    
    ddt(n)[i] = 
      - bracket_arakawa(phi, n, i) 
      + (2./R_c)*DDZ(n, i) 
      + n[i]*phi[i]/L_par
      ;
    
    ddt(vort)[i] = 
      - bracket_arakawa(phi, vort, i)
      + (2./R_c)*DDZ(n, i)
      + phi[i]/L_par;
      
  }
  Duration elapsed2 = steady_clock::now() - start2;
  // Save ddt values for comparison later
  Field3D dn_2 = ddt(n);
  Field3D dvort_2 = ddt(vort);

  /////////////////////////////////////////////////
  // Using several internal loops over domain 
  
  SteadyClock start1 = steady_clock::now();
  ddt(n) = 
    - bracket(phi, n, BRACKET_ARAKAWA)
    + (2./R_c)*DDZ(n)
    + n*phi/L_par
    ;
  
  ddt(vort) = 
      - bracket(phi, vort, BRACKET_ARAKAWA)
      + (2./R_c)*DDZ(n)
      + phi/L_par;
  Duration elapsed1 = steady_clock::now() - start1;
  
  /////////////////////////////////////////////////
  // Check results
  
  for(auto i : n.region(RGN_NOBNDRY)) {
    if(abs(ddt(n)[i] - dn_2[i]) > 1e-5) {
      output.write("Difference in ddt(n) at (%d,%d,%d): %e, %e\n", 
                   i.x, i.y, i.z, ddt(n)[i], dn_2[i]);
    }
    if(abs(ddt(vort)[i] - dvort_2[i]) > 1e-5) {
      output.write("Difference in ddt(vort) at (%d,%d,%d): %e, %e\n", 
                   i.x, i.y, i.z, ddt(vort)[i], dvort_2[i]);
    }
  }
  
  output << "TIMING\n======\n";
  output << "Inner loops          : " << elapsed1.count() << std::endl;
  output << "Outer loop           : " << elapsed2.count() << std::endl;
  
  BoutFinalise();
  return 0;
}
