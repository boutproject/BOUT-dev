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
 * Inner loops          : 0.353403
 * Outer loop           : 0.186931
 * 
 */

#include <bout.hxx>
#include <bout/derivs.hxx>

#include <chrono>
#include <iostream>

/////////////////////////////////////////////////////////////////////////////
// 2nd order central differencing

/*!
 * \brief 2nd order central differencing in X
 * 
 * Performs calculation at a single point. The input field
 * must be defined at i.xp() and i.xm()
 * 
 * @param f  The field to be differentiated
 * @param i  The point where the result is calculated
 */
inline BoutReal DDX_C2(const Field3D &f, const DataIterator &i) {
  return (f[i.xp()] - f[i.xm()])/(2.*mesh->coordinates()->dx[i]);
}

inline BoutReal DDY_C2(const Field3D &f, const DataIterator &i) {
  return (f.yup()[i.yp()] - f.ydown()[i.ym()])/(2.*mesh->coordinates()->dy[i]);
}


inline BoutReal DDZ_C2(const Field3D &f, const DataIterator &i) {
  return (f[i.zp()] - f[i.zm()])/(2.*mesh->coordinates()->dz);
}

/*!
 * 2nd order central differencing in X
 * Performs operation over RGN_NOBNDRY
 */
const Field3D DDX_C2(const Field3D &f) {
  Field3D result;
  result.allocate();
  for(auto i : result.region(RGN_NOBNDRY)) {
    result[i] = DDX_C2(f, i);
  }
  return result;
}

/*!
 * 2nd order central differencing in Y
 * Performs operation over RGN_NOBNDRY
 */
const Field3D DDY_C2(const Field3D &f) {
  Field3D result;
  result.allocate();
  for(auto i : result.region(RGN_NOBNDRY)) {
    result[i] = DDY_C2(f, i);
  }
  return result;
}

/*!
 * 2nd order central differencing in Z
 * Performs operation over RGN_NOBNDRY
 */
const Field3D DDZ_C2(const Field3D &f) {
  Field3D result;
  result.allocate();
  for(auto i : result.region(RGN_NOBNDRY)) {
    result[i] = DDZ_C2(f, i);
  }
  return result;
}

/////////////////////////////////////////////////////////////////////////////
// 4th order central differencing

/*!
 * 4th order central differencing in X
 */
inline BoutReal DDX_C4(const Field3D &f, const DataIterator &i) {
  return (8.*f[i.xp()] - 8.*f[i.xm()] + f[i.offset(-2,0,0)] - f[i.offset(2,0,0)])/(12.*mesh->coordinates()->dx[i]);
}

/*!
 * 4th order central differencing in Z
 */
inline BoutReal DDZ_C4(const Field3D &f, const DataIterator &i) {
  return (8.*f[i.zp()] - 8.*f[i.zm()] + f[i.offset(0,0,-2)] - f[i.offset(0,0,2)])/(12.*mesh->coordinates()->dz);
}

/*!
 * 4th order central differencing in X
 * Performs operation over RGN_NOBNDRY
 */
const Field3D DDX_C4(const Field3D &f) {
  Field3D result;
  result.allocate();
  for(auto i : result.region(RGN_NOBNDRY)) {
    result[i] = DDX_C4(f, i);
  }
  return result;
}

/*!
 * 4th order central differencing in X
 * Performs operation over RGN_NOBNDRY
 */
const Field3D DDZ_C4(const Field3D &f) {
  Field3D result;
  result.allocate();
  for(auto i : result.region(RGN_NOBNDRY)) {
    result[i] = DDZ_C4(f, i);
  }
  return result;
}

/////////////////////////////////////////////////////////////////////////////
// First derivatives

const Field3D DDX_test(const Field3D &f) {
  // This could be done with a static function
  // but a) this has another "if" to check the function, and 
  //     b) is not thread safe (OpenMP)
  static const Field3D (*func)(const Field3D &f) = nullptr; 
  
  if(!func) {
    // Not defined yet
    string setting;
    Options::getRoot()->getSection("operators")->get("ddx", setting, "c2");
    setting = lowercase(setting);
    
    if(setting == "c2") {
      func = &DDX_C2;
    }else if(setting == "c4") {
      func = &DDX_C4;
    }else {
      throw BoutException("Unrecognised option for DDX: '%s'", setting.c_str());
    }
  }
  // Call the function
  return (*func)(f);
}



/////////////////////////////////////////////////////////////////////////////
// Arakawa brackets

BoutReal bracket_arakawa(const Field3D &f, const Field3D &g, const DataIterator &i) {
  Coordinates *metric = mesh->coordinates();
  
  BoutReal fxp = f[i.xp()];
  BoutReal fxm = f[i.xm()];
  BoutReal fzp = f[i.zp()];
  BoutReal fzm = f[i.zm()];

  BoutReal fpp = f[i.offset(1,0,1)];
  BoutReal fpm = f[i.offset(1,0,-1)];
  BoutReal fmp = f[i.offset(-1,0,1)];
  BoutReal fmm = f[i.offset(-1,0,-1)];
  
  BoutReal gxp = g[i.xp()];
  BoutReal gxm = g[i.xm()];
  BoutReal gzp = g[i.zp()];
  BoutReal gzm = g[i.zm()];

  // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
  BoutReal Jpp = 0.25*( (fzp - fzm)*
                        (gxp - gxm) -
                        (fxp - fxm)*
                        (gzp - gzm) );
  
  // J+x
  BoutReal Jpx = 0.25*( gxp*(fpp-fpm) -
                        gxm*(fmp-fmm) -
                        gzp*(fpp-fmp) +
                        gzm*(fpm-fmm));
  // Jx+
  BoutReal Jxp = 0.25*( g[i.offset(1,0,1)]*(fzp-fxp) -
                        g[i.offset(-1,0,-1)]*(fxm-fzm) -
                        g[i.offset(-1,0,1)]*(fzp-fxm) +
                        g[i.offset(1,0,-1)]*(fxp-fzm) );
  
  return (Jpp + Jpx + Jxp) / (3. * metric->dx[i] * metric->dz);
}

const Field3D bracket_arakawa(const Field3D &f, const Field3D &g) {
  Field3D result;
  result.allocate();
  for(auto i : result.region(RGN_NOBNDRY)) {
    result[i] = bracket_arakawa(f, g, i);
  }
  return result;
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
      + (2./R_c)*DDZ_C2(n, i) 
      + n[i]*phi[i]/L_par
      ;
    
    ddt(vort)[i] = 
      - bracket_arakawa(phi, vort, i)
      + (2./R_c)*DDZ_C2(n, i)
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
    - bracket_arakawa(phi, n) //bracket(phi, n, BRACKET_ARAKAWA)
    + (2./R_c)*DDZ(n)
    + n*phi/L_par
    ;
  
  ddt(vort) = 
    - bracket_arakawa(phi, vort)//bracket(phi, vort, BRACKET_ARAKAWA)
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
