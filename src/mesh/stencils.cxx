
/*******************************************************************************
 * Basic derivative methods.
 * All expect to have an input grid cell at the same location as the output
 * Hence convert cell centred values -> centred values, or left -> left
 *******************************************************************************/

const BoutReal WENO_SMALL = 1.0e-8; // Small number for WENO schemes

////////////////////// FIRST DERIVATIVES /////////////////////

/// central, 2nd order
BoutReal DDX_C2(stencil &f) {
  return 0.5*(f.p - f.m);
}

/// forward, 2nd order
Mesh::boundary_derivs_pair DDX_F2(forward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = -1.5*f.c+2.*f.p-0.5*f.p2;
  result.outer = -1.5*f.m+2.*f.c-0.5*f.p;
  return result;
}

/// backward 2nd order
Mesh::boundary_derivs_pair DDX_B2(backward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 0.5*f.m2-2.*f.m+1.5*f.c;
  result.outer = 0.5*f.m-2.*f.c+1.5*f.p;
  return result;
}

/// central, 4th order
BoutReal DDX_C4(stencil &f) {
  return (8.*f.p - 8.*f.m + f.mm - f.pp)/12.;
}

/// forward, 4th order
Mesh::boundary_derivs_pair DDX_F4(forward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = -1./4.*f.m-5./6.*f.c+3./2.*f.p-1./2.*f.p2+1./12.*f.p3; // uncentred (forward-biased) derivative
  result.outer = -25./12.*f.m+4.*f.c-3.*f.p+4./3.*f.p2-1./4.*f.p3; // forward derivative
  return result;
}

/// backward, 4th order
Mesh::boundary_derivs_pair DDX_B4(backward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 1./4.*f.p+5./6.*f.c-3./2.*f.m+1./2.*f.m2-1./12.*f.m3; // uncentred (forward-biased) derivative
  result.outer = 25./12.*f.p-4.*f.c+3.*f.m-4./3.*f.m2+1./4.*f.m3; // forward derivative
  return result;
}

/// Central WENO method, 2nd order (reverts to 1st order near shocks)
BoutReal DDX_CWENO2(stencil &f) {
  BoutReal isl, isr, isc; // Smoothness indicators
  BoutReal al, ar, ac, sa; // Un-normalised weights
  BoutReal dl, dr, dc; // Derivatives using different stencils
  
  dc = 0.5*(f.p - f.m);
  dl = f.c - f.m;
  dr = f.p - f.c;
  
  isl = SQ(dl);
  isr = SQ(dr);
  isc = (13./3.)*SQ(f.p - 2.*f.c + f.m) + 0.25*SQ(f.p-f.m);
  
  al = 0.25/SQ(WENO_SMALL + isl);
  ar = 0.25/SQ(WENO_SMALL + isr);
  ac = 0.5/SQ(WENO_SMALL + isc);
  sa = al + ar + ac;

  return (al*dl + ar*dr + ac*dc)/sa;
}

// Smoothing 2nd order derivative
BoutReal DDX_S2(stencil &f) {
  
  // 4th-order differencing
  BoutReal result = (8.*f.p - 8.*f.m + f.mm - f.pp)/12.;

  result += SIGN(f.c)*(f.pp - 4.*f.p + 6.*f.c - 4.*f.m + f.mm)/12.;
  
  return result;
}

///////////////////// SECOND DERIVATIVES ////////////////////

/// Second derivative: Central, 2nd order
BoutReal D2DX2_C2(stencil &f) {
  return f.p + f.m - 2.*f.c;
}

/// Second derivative: Forward, 2nd order
Mesh::boundary_derivs_pair D2DX2_F2(forward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 2.*f.c-5.*f.p+4.*f.p2-3.*f.p3;
  result.outer = 2.*f.m-5.*f.c+4.*f.p-3.*f.p2;
  return result;
}

/// Second derivative: Backward, 2nd order
Mesh::boundary_derivs_pair D2DX2_B2(backward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = -2.*f.c+5.*f.m-4.*f.m2+3.*f.m3;
  result.outer = -2.*f.p+5.*f.c-4.*f.m+3.*f.m2;
  return result;
}

/// Second derivative: Central, 4th order
BoutReal D2DX2_C4(stencil &f) {
  return (-f.pp + 16.*f.p - 30.*f.c + 16.*f.m - f.mm)/12.;
}

/// Second derivatives: Forward, 4th order
Mesh::boundary_derivs_pair D2DX2_F4(forward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 5./6.*f.m-5./4.*f.c-1./3.*f.p+7./6.*f.p2-1./2.*f.p3+1./12.*f.p4; // uncentred (forward-biased) derivative
  result.outer = 15./4.*f.m-77./6.*f.c+107./6.*f.p-13.*f.p2+61./12.*f.p3-5./6.*f.p4; // forward derivative
  return result;
}

/// Second derivatives: Backward, 4th order
Mesh::boundary_derivs_pair D2DX2_B4(backward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 5./6.*f.p-5./4.*f.c-1./3.*f.m+7./6.*f.m2-1./2.*f.m3+1./12.*f.m4; // uncentred (backward-biased) derivative
  result.outer = 15./4.*f.p-77./6.*f.c+107./6.*f.m-13.*f.m2+61./12.*f.m3-5./6.*f.m4; // backward derivative
  return result;
}

//////////////////////// UPWIND METHODS ///////////////////////

/// Upwinding: Central, 2nd order
BoutReal VDDX_C2(stencil &v, stencil &f) {
  return v.c*0.5*(f.p - f.m);
}

/// Upwinding: Central, 4th order
BoutReal VDDX_C4(stencil &v, stencil &f) {
  return v.c*(8.*f.p - 8.*f.m + f.mm - f.pp)/12.;
}

/// upwind, 1st order
BoutReal VDDX_U1(stencil &v, stencil &f) {
  return v.c>=0.0 ? v.c*(f.c - f.m): v.c*(f.p - f.c);
}

/// upwind, 2nd order
BoutReal VDDX_U2(stencil &v, stencil &f) {
  return v.c>=0.0 ? v.c*(1.5*f.c - 2.0*f.m + 0.5*f.mm): v.c*(-0.5*f.pp + 2.0*f.p - 1.5*f.c);
}

/// upwind, 4th order
BoutReal VDDX_U4(stencil &v, stencil &f) {
  return v.c >= 0.0 ? v.c*(4.*f.p - 12.*f.m + 2.*f.mm + 6.*f.c)/12.
    : v.c*(-4.*f.m + 12.*f.p - 2.*f.pp - 6.*f.c)/12.;
}

/// TVD upwinding (2nd order)
/// WARNING WARNING : THIS TVD IMPLEMENTATION DOESN'T WORK PROPERLY
BoutReal VDDX_TVD(BoutReal vc, BoutReal vm, BoutReal vp, BoutReal fc, BoutReal fm, BoutReal fp, BoutReal fmm, BoutReal fpp) {
  BoutReal denom, res, ri, ri1, ri_1, fluxLeft, fluxRight;

  if (vc>=0.0){

    /*-smoothness indicator*/
    denom=fp-fc;
    if (fabs(denom) < 1e-20) denom=1e-20;
    ri=(fc-fm)/denom;

    denom=fc-fm;
    if (fabs(denom) < 1e-20) denom=1e-20;
    ri_1=(fm-fmm)/denom;

    /*;-nonlinear TVD flux*/
    fluxRight=fc + 0.5*VANLEER(ri)*(fp-fc);
    fluxLeft=fm + 0.5*VANLEER(ri_1)*(fc-fm);
    res = vc*(fluxRight-fluxLeft); /*divide by dx outside*/
    
  }
  else{

    /*;-smoothness indicator*/
    denom=fm-fc;
    if (fabs(denom) < 1e-20) denom=1e-20;
    ri=(fc-fp)/denom;
    
    denom=(fc-fp);
    if (fabs(denom) < 1e-20) denom=1e-20;
    ri1=(fp-fpp)/denom;

    /*;-nonlinear TVD flux*/
    fluxRight=(fp - 0.5*VANLEER(ri1)*(fp-fc));
    fluxLeft=(fc - 0.5*VANLEER(ri)*(fc-fm));
    res = vc*(fluxRight-fluxLeft); /*divide by dx outside*/    
  }
  return res;
}

/// 3rd-order WENO scheme
BoutReal VDDX_WENO3(stencil &v, stencil &f) {
  BoutReal deriv, w, r;

  if(v.c > 0.0) {
    // Left-biased stencil
    
    r = (WENO_SMALL + SQ(f.c - 2.0*f.m + f.mm)) / (WENO_SMALL + SQ(f.p - 2.0*f.c + f.m));
    w = 1.0 / (1.0 + 2.0*r*r);

    deriv = 0.5*(f.p - f.m) - 0.5*w*(-f.mm + 3.*f.m - 3.*f.c + f.p);
    
  }else {
    // Right-biased
    
    r = (WENO_SMALL + SQ(f.pp - 2.0*f.p + f.c)) / (WENO_SMALL + SQ(f.p - 2.0*f.c + f.m));
    w = 1.0 / (1.0 + 2.0*r*r);
    
    deriv = 0.5*(f.p - f.m) - 0.5*w*( -f.m + 3.*f.c - 3.*f.p + f.pp );
  }

  return v.c*deriv;
}

/// 3rd-order CWENO. Uses the upwinding code and split flux
BoutReal DDX_CWENO3(stencil &f) {
  BoutReal a, ma = fabs(f.c);
  // Split flux
  a = fabs(f.m); if(a > ma) ma = a;
  a = fabs(f.p); if(a > ma) ma = a;
  a = fabs(f.mm); if(a > ma) ma = a;
  a = fabs(f.pp); if(a > ma) ma = a;
  
  stencil sp, vp, sm, vm;
  
  vp.c = 0.5; vm.c = -0.5;
  
  sp = f + ma;
  sm = ma - f;
  
  return VDDX_WENO3(vp, sp) + VDDX_WENO3(vm, sm);
}

//////////////////////// FLUX METHODS ///////////////////////

BoutReal FDDX_U1(stencil &v, stencil &f) {
  // Velocity at lower end
  BoutReal vs = 0.5*(v.m + v.c);
  BoutReal result = (vs >= 0.0) ? vs * f.m : vs * f.c;
  // and at upper 
  vs = 0.5*(v.c + v.p);
  result -= (vs >= 0.0) ? vs * f.c : vs * f.p;

  return result;
}

BoutReal FDDX_C2(stencil &v, stencil &f) {
  return 0.5*(v.p*f.p - v.m*f.m);
}

BoutReal FDDX_C4(stencil &v, stencil &f) {
  return (8.*v.p*f.p - 8.*v.m*f.m + v.mm*f.mm - v.pp*f.pp)/12.;
}

/// Non-oscillatory, containing No free parameters and Dissipative (NND) scheme
/// http://arxiv.org/abs/1010.4135v1
BoutReal FDDX_NND(stencil &v, stencil &f) {
  // f{+-} i
  BoutReal fp = 0.5*(v.c + fabs(v.c))*f.c;
  BoutReal fm = 0.5*(v.c - fabs(v.c))*f.c;
  
  // f{+-} i+1
  BoutReal fp1 = 0.5*(v.p + fabs(v.p))*f.p;
  BoutReal fm1 = 0.5*(v.p - fabs(v.p))*f.p;
  
  // f{+-} i+2
  BoutReal fm2 = 0.5*(v.pp - fabs(v.pp))*f.pp;

  // f{+-} i-1
  BoutReal fp_1 = 0.5*(v.m + fabs(v.m))*f.m;
  BoutReal fm_1 = 0.5*(v.m - fabs(v.m))*f.m;
  
  // f{+-} i-2
  BoutReal fp_2 = 0.5*(v.mm + fabs(v.mm))*f.mm;

  // f^{LR} {i+1/2}
  BoutReal flp = fp  + 0.5*MINMOD(fp1 - fp, fp - fp_1);
  BoutReal frp = fm1 - 0.5*MINMOD(fm1 - fm, fm2 - fm1);
  
  // f^{LR} {i-1/2}
  BoutReal flm = fp_1  + 0.5*MINMOD(fp - fp_1, fp_1 - fp_2);
  BoutReal frm = fm - 0.5*MINMOD(fm - fm_1, fm1 - fm);
    
  // h{+-}
  BoutReal hp = flp + frp;
  BoutReal hm = flm + frm; 
  
  return hp - hm;
}

//////////////////////// MUSCL scheme ///////////////////////

void DDX_KT_LR(const stencil &f, BoutReal &fLp, BoutReal &fRp, BoutReal &fLm, BoutReal &fRm) {
  // Limiter functions
  BoutReal phi   = SUPERBEE( (f.c - f.m) / (f.p - f.c) );
  BoutReal phi_m = SUPERBEE( (f.m - f.mm) / (f.c - f.m) );
  BoutReal phi_p = SUPERBEE( (f.p - f.c) / (f.pp - f.p) );

  fLp = f.c + 0.5*phi*(f.p - f.c);
  fRp = f.p - 0.5*phi_p*(f.pp - f.p);

  fLm = f.m + 0.5*phi_m*(f.c - f.m);
  fRm = f.c - 0.5*phi*(f.p - f.c);
}

// du/dt = d/dx(f)  with maximum local velocity Vmax
BoutReal DDX_KT(const stencil &f, const stencil &u, const BoutReal Vmax) {
  BoutReal uLp, uRp, uLm, uRm;
  BoutReal fLp, fRp, fLm, fRm;
  
  DDX_KT_LR(u, uLp, uRp, uLm, uRm);
  DDX_KT_LR(f, fLp, fRp, fLm, fRm);
  
  BoutReal Fm = 0.5*( fRm + fLm - Vmax*(uRm - uLm));
  BoutReal Fp = 0.5*( fRp + fLp - Vmax*(uRp - uLp));
  
  return Fm - Fp;
}

/*******************************************************************************
 * Staggered differencing methods
 * These expect the output grid cell to be at a different location to the input
 * 
 * The stencil no longer has a value in 'C' (centre)
 * instead, points are shifted as follows:
 * 
 * mm  -> -3/2 h
 * m   -> -1/2 h
 * p   -> +1/2 h
 * pp  -? +3/2 h
 *
 * NOTE: Cell widths (dx, dy, dz) are currently defined as centre->centre
 * for the methods above. This is currently not taken account of, so large
 * variations in cell size will cause issues.
 *******************************************************************************/

/////////////////////// FIRST DERIVATIVES //////////////////////
// Map Centre -> Low or Low -> Centre

// Second order differencing (staggered)
BoutReal DDX_C2_stag(stencil &f) {
  return f.p - f.m;
}

Mesh::boundary_derivs_pair DDX_F2_stag(forward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = -2.*f.c+3*f.p-f.p2;
  result.outer = -2.*f.m+3*f.c-f.p;
  return result;
}

Mesh::boundary_derivs_pair DDX_B2_stag(backward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 2.*f.c-3*f.m+f.m2;
  result.outer = 2.*f.p-3*f.c+f.m;
  return result;
}

BoutReal DDX_C4_stag(stencil &f) {
  return ( 27.*(f.p - f.m) - (f.pp - f.mm) ) / 24.;
}

Mesh::boundary_derivs_pair DDX_F4_stag(forward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = -11./12.*f.m+17./24.*f.c+3./8.*f.p-5./24.*f.p2+1./24.*f.p3; // uncentred (forward-biased) derivative
  result.outer = -31./8*f.m+229./24.*f.c-75./8.*f.p+37./8.*f.p2-11./12.*f.p3; // forward derivative
  return result;
}

Mesh::boundary_derivs_pair DDX_B4_stag(backward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 11./12.*f.p-17./24.*f.c-3./8.*f.m+5./24.*f.m2-1./24.*f.m3; // uncentred (backward-biased) derivative
  result.outer = 31./8*f.p-229./24.*f.c+75./8.*f.m-37./8.*f.m2+11./12.*f.m3; // backward derivative
  return result;
}

/////////////////////// SECOND DERIVATIVES //////////////////////
// Map Centre -> Low or Low -> Centre

Mesh::boundary_derivs_pair D2DX2_F2_stag(forward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 5./2.*f.c-13./2.*f.p+11./2.*f.p2-3./2.*f.p3;
  result.outer = 5./2.*f.m-13./2.*f.c+11./2.*f.p-3./2.*f.p2;
  return result;
}

Mesh::boundary_derivs_pair D2DX2_B2_stag(backward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 5./2.*f.c-13./2.*f.m+11./2.*f.m2-3./2.*f.m3;
  result.outer = 5./2.*f.p-13./2.*f.c+11./2.*f.m-3./2.*f.m2;
  return result;
}

BoutReal D2DX2_C4_stag(stencil &f) {
  return ( f.pp + f.mm - f.p - f.m ) / 2.;
}

Mesh::boundary_derivs_pair D2DX2_F4_stag(forward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 95./48.*f.m-269./48.*f.c+49./8.*f.p-85./24.*f.p2+59./48.*f.p3-3./16.*f.p4;
  result.outer = 301./48.*f.m-377./16.*f.c+865./24.*f.p-683./24.*f.p2+187./16.*f.p3-95./48.*f.p4;
  return result;
}

Mesh::boundary_derivs_pair D2DX2_B4_stag(backward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 95./48.*f.p-269./48.*f.c+49./8.*f.m-85./24.*f.m2+59./48.*f.m3-3./16.*f.m4;
  result.outer = 301./48.*f.p-377./16.*f.c+865./24.*f.m-683./24.*f.m2+187./16.*f.m3-95./48.*f.m4;
  return result;
}

/////////////////////////// UPWINDING ///////////////////////////
// Map (Low, Centre) -> Centre  or (Centre, Low) -> Low
// Hence v contains only (mm, m, p, pp) fields whilst f has 'c' too
//
// v.p is v at +1/2, v.m is at -1/2

BoutReal VDDX_U1_stag(stencil &v, stencil &f) {
  // Lower cell boundary
  BoutReal result = (v.m >= 0) ? v.m * f.m : v.m * f.c;
  
  // Upper cell boundary
  result -= (v.p >= 0) ? v.p * f.c : v.p * f.p;
  
  result *= -1;
  
  // result is now d/dx(v*f), but want v*d/dx(f) so subtract f*d/dx(v)
  result -= f.c*(v.p - v.m);
  
  return result;
}

BoutReal VDDX_U2_stag(stencil &v, stencil &f) {
  BoutReal result;
  
  if (v.p>0 && v.m>0) {
    // Extrapolate v to centre from below, use 2nd order backward difference on f
    result = (1.5*v.m - .5*v.mm) * (.5*f.mm - 2.*f.m + 1.5*f.c);
  }
  else if (v.p<0 && v.m<0) {
    // Extrapolate v to centre from above, use 2nd order forward difference on f
    result = (1.5*v.p - .5*v.pp) * (-1.5*f.c + 2.*f.p - .5*f.pp);
  }
  else {
    // Velocity changes sign, hence is almost zero: use centred interpolation/differencing
    result = .25 * (v.p + v.m) * (f.p - f.m);
  }
  
  return result;
}

BoutReal VDDX_C2_stag(stencil &v, stencil &f) {
  // Result is needed at location of f: interpolate v to f's location and take an unstaggered derivative of f
  return 0.5*(v.p+v.m) * 0.5*(f.p - f.m);
}

BoutReal VDDX_C4_stag(stencil &v, stencil &f) {
  // Result is needed at location of f: interpolate v to f's location and take an unstaggered derivative of f
  return (9.*(v.m + v.p) - v.mm - v.pp)/16. * (8.*f.p - 8.*f.m + f.mm - f.pp)/12.;
}

/////////////////////////// FLUX ///////////////////////////
// Map (Low, Centre) -> Centre  or (Centre, Low) -> Low
// Hence v contains only (mm, m, p, pp) fields whilst f has 'c' too
//
// v.p is v at +1/2, v.m is at -1/2

BoutReal FDDX_U1_stag(stencil &v, stencil &f) {
  // Lower cell boundary
  BoutReal result = (v.m >= 0) ? v.m * f.m : v.m * f.c;
  
  // Upper cell boundary
  result -= (v.p >= 0) ? v.p * f.c : v.p * f.p;

  return result;
}
