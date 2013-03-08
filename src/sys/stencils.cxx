/**************************************************************************
 * Stencil class for differencing, and handles indexing
 * NOTE: bstencil now depreciated
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 * 
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#include <cmath>

typedef double BoutReal;

#include <globals.hxx>
#include <stencils.hxx>
#include <output.hxx>

/**************************************************************************
 * bvalue class
 **************************************************************************/

// Operators
  
bvalue & bvalue::operator=(const bvalue &rhs)
{
  val = rhs.val;
  jx = rhs.jx;
  jy = rhs.jy;
  jz = rhs.jz;
  return(*this);
}

bvalue & bvalue::operator+=(const bvalue &rhs)
{
  // Should check that jx, jy and jz are the same
  val += rhs.val;
  return(*this);
}

bvalue & bvalue::operator-=(const bvalue &rhs)
{
  val -= rhs.val;
  return(*this);
}

bvalue & bvalue::operator*=(const bvalue &rhs)
{
  val *= rhs.val;
  return(*this);
}

bvalue & bvalue::operator*=(const BoutReal rhs)
{
  val *= rhs;
  return(*this);
}

bvalue & bvalue::operator/=(const bvalue &rhs)
{
  val /= rhs.val;
  return(*this);
}

bvalue & bvalue::operator/=(const BoutReal rhs)
{
  val /= rhs;
  return(*this);
}

// Binary operators
  
const bvalue bvalue::operator+(const bvalue &other) const
{
  bvalue result = *this;
  result += other;
  return(result);
}

const bvalue bvalue::operator-(const bvalue &other) const
{
  bvalue result = *this;
  result -= other;
  return(result);
}

const bvalue bvalue::operator*(const bvalue &other) const
{
  bvalue result = *this;
  result *= other;
  return(result);
}

const bvalue bvalue::operator*(const BoutReal rhs) const
{
  bvalue result = *this;
  result *= rhs;
  return(result);
}

const bvalue bvalue::operator/(const bvalue &other) const
{
  bvalue result = *this;
  result /= other;
  return(result);
}

const bvalue bvalue::operator/(const BoutReal rhs) const
{
  bvalue result = *this;
  result /= rhs;
  return(result);
}

const bvalue operator*(const BoutReal lhs, const bvalue &rhs)
{
  return(rhs * lhs);
}

const bvalue operator/(const BoutReal lhs, const bvalue &rhs)
{
  bvalue result = rhs;
  result.val = lhs / result.val;
  return(result);
}

/*******************************************************************************
 * stencil class
 *******************************************************************************/

stencil::stencil()
{
  jx = jy = jz = 0;
  c = p = m = pp = mm = 0.0;
}

stencil::stencil(const stencil &s)
{
  *this = s;
}

stencil::stencil(BoutReal fc)
{
  c = p = m = pp = mm = fc;
}

stencil::stencil(BoutReal fc, BoutReal fm, BoutReal fp, BoutReal fmm, BoutReal fpp)
{
  c  = fc;
  m  = fm;
  p  = fp;
  mm = fmm;
  pp = fpp;
}

stencil & stencil::operator=(const stencil &s)
{
  // Check for self-assignment
  if(this == &s)
    return(*this); // skip this assignment
  
  jx = s.jx;
  jy = s.jy;
  jz = s.jz;

  c = s.c;
  p = s.p;
  m = s.m;
  pp = s.pp;
  mm = s.mm;
  
  return *this;
}

stencil & stencil::operator=(const BoutReal rhs)
{
  jx = jy = jz = 0;

  c = p = m = pp = mm = rhs;

  return *this;
}

stencil & stencil::operator+=(const stencil &s)
{
#ifdef CHECK
  if((jx != s.jx) || (jy != s.jy) || (jz != s.jz)) {
    output.write("Error: Stencil += not at same location\n");
    exit(1);
  }
#endif

  c += s.c;
  p += s.p;
  m += s.m;
  pp += s.pp;
  mm += s.mm;

  return *this;
}

stencil & stencil::operator+=(const BoutReal &rhs)
{
  c += rhs;
  p += rhs;
  m += rhs;
  pp += rhs;
  mm += rhs;

  return *this;
}

stencil & stencil::operator-=(const stencil &s)
{
#ifdef CHECK
  if((jx != s.jx) || (jy != s.jy) || (jz != s.jz)) {
    output.write("Error: Stencil -= not at same location\n");
    exit(1);
  }
#endif

  c -= s.c;
  p -= s.p;
  m -= s.m;
  pp -= s.pp;
  mm -= s.mm;

  return *this;
}

stencil & stencil::operator-=(const BoutReal &rhs)
{
  c -= rhs;
  p -= rhs;
  m -= rhs;
  pp -= rhs;
  mm -= rhs;

  return *this;
}

stencil & stencil::operator*=(const stencil &s)
{
#ifdef CHECK
  if((jx != s.jx) || (jy != s.jy) || (jz != s.jz)) {
    output.write("Error: Stencil *= not at same location\n");
    exit(1);
  }
#endif

  c *= s.c;
  p *= s.p;
  m *= s.m;
  pp *= s.pp;
  mm *= s.mm;

  return *this;
}

stencil & stencil::operator*=(const BoutReal &rhs)
{
  c *= rhs;
  p *= rhs;
  m *= rhs;
  pp *= rhs;
  mm *= rhs;

  return *this;
}

stencil & stencil::operator/=(const stencil &s)
{
#ifdef CHECK
  if((jx != s.jx) || (jy != s.jy) || (jz != s.jz)) {
    output.write("Error: Stencil /= not at same location\n");
    exit(1);
  }
#endif

  c /= s.c;
  p /= s.p;
  m /= s.m;
  pp /= s.pp;
  mm /= s.mm;

  return *this;
}

stencil & stencil::operator/=(const BoutReal &rhs)
{
  c /= rhs;
  p /= rhs;
  m /= rhs;
  pp /= rhs;
  mm /= rhs;

  return *this;
}

const stencil stencil::operator+(const stencil &other) const
{
  stencil result = *this;
  result += other;
  return result;
}

const stencil stencil::operator+(const BoutReal &other) const
{
  stencil result = *this;
  result += other;
  return result;
}

const stencil stencil::operator-(const stencil &other) const
{
  stencil result = *this;
  result -= other;
  return result;
}

const stencil stencil::operator-(const BoutReal &other) const
{
  stencil result = *this;
  result -= other;
  return result;
}

const stencil stencil::operator*(const stencil &other) const
{
  stencil result = *this;
  result *= other;
  return result;
}

const stencil stencil::operator*(const BoutReal &other) const
{
  stencil result = *this;
  result *= other;
  return result;
}

const stencil stencil::operator/(const stencil &other) const
{
  stencil result = *this;
  result /= other;
  return result;
}

const stencil stencil::operator/(const BoutReal &other) const
{
  stencil result = *this;
  result /= other;
  return result;
}

const stencil operator+(const BoutReal &lhs, const stencil &rhs)
{
  return rhs + lhs;
}

const stencil operator-(const BoutReal &lhs, const stencil &rhs)
{
  stencil result;
  
  result.c = lhs - rhs.c;
  result.p = lhs - rhs.p;
  result.m = lhs - rhs.m;
  result.pp = lhs - rhs.pp;
  result.mm = lhs - rhs.mm;

  return result;
}

const stencil operator*(const BoutReal &lhs, const stencil &rhs)
{
  return rhs * lhs;
}

const stencil operator/(const BoutReal &lhs, const stencil &rhs)
{
  stencil result;
  
  result.c = lhs / rhs.c;
  result.p = lhs / rhs.p;
  result.m = lhs / rhs.m;
  result.pp = lhs / rhs.pp;
  result.mm = lhs / rhs.mm;

  return result;
}

BoutReal stencil::min() const
{
  BoutReal r;

  r = c;
  if(p < r)
    r = p;
  if(m < r)
    r = m;
  if(pp < r)
    r = pp;
  if(mm < r)
    r = mm;
  return r;
}

BoutReal stencil::max() const
{
  BoutReal r;

  r = c;
  if(p > r)
    r = p;
  if(m > r)
    r = m;
  if(pp > r)
    r = pp;
  if(mm > r)
    r = mm;
  return r;
}

const stencil stencil::abs() const
{
  stencil result;
  
  result.c = fabs(c);
  result.p = fabs(p);
  result.m = fabs(m);
  result.pp = fabs(pp);
  result.mm = fabs(mm);

  return result;
}

BoutReal min(const stencil &s)
{
  return s.min();
}

BoutReal max(const stencil &s)
{
  return s.max();
}

const stencil abs(const stencil &s)
{
  return s.abs();
}

/*******************************************************************************
 * bstencil class
 *******************************************************************************/

////////////////// Operators ///////////////////////

bstencil & bstencil::operator=(const bstencil &rhs)
{
  jx = rhs.jx;
  jy = rhs.jy;
  jz = rhs.jz;

  cc = rhs.cc;

  xp = rhs.xp; xm = rhs.xm;
  yp = rhs.yp; ym = rhs.ym;
  zp = rhs.zp; zm = rhs.zm;

  x2p = rhs.x2p; x2m = rhs.x2m;
  y2p = rhs.y2p; y2m = rhs.y2m;
  z2p = rhs.z2p; z2m = rhs.z2m;

  return(*this);
}

bstencil & bstencil::operator+=(const bstencil &rhs)
{
  // Should check that jx,jy and jz are the same

  cc += rhs.cc;
  
  xp += rhs.xp; xm += rhs.xm;
  yp += rhs.yp; ym += rhs.ym;
  zp += rhs.zp; zm += rhs.zm;

  x2p += rhs.x2p; x2m += rhs.x2m;
  y2p += rhs.y2p; y2m += rhs.y2m;
  z2p += rhs.z2p; z2m += rhs.z2m;

  return(*this);
}

bstencil & bstencil::operator-=(const bstencil &rhs)
{
  cc -= rhs.cc;
  
  xp -= rhs.xp; xm -= rhs.xm;
  yp -= rhs.yp; ym -= rhs.ym;
  zp -= rhs.zp; zm -= rhs.zm;

  x2p -= rhs.x2p; x2m -= rhs.x2m;
  y2p -= rhs.y2p; y2m -= rhs.y2m;
  z2p -= rhs.z2p; z2m -= rhs.z2m;

  return(*this);
}

bstencil & bstencil::operator*=(const bstencil &rhs)
{
  cc *= rhs.cc;
  
  xp *= rhs.xp; xm *= rhs.xm;
  yp *= rhs.yp; ym *= rhs.ym;
  zp *= rhs.zp; zm *= rhs.zm;

  x2p *= rhs.x2p; x2m *= rhs.x2m;
  y2p *= rhs.y2p; y2m *= rhs.y2m;
  z2p *= rhs.z2p; z2m *= rhs.z2m;

  return(*this);
}

bstencil & bstencil::operator*=(const BoutReal rhs)
{
  cc *= rhs;
  
  xp *= rhs; xm *= rhs;
  yp *= rhs; ym *= rhs;
  zp *= rhs; zm *= rhs;

  x2p *= rhs; x2m *= rhs;
  y2p *= rhs; y2m *= rhs;
  z2p *= rhs; z2m *= rhs;

  return(*this);
}

bstencil & bstencil::operator/=(const bstencil &rhs)
{
  cc /= rhs.cc;
  
  xp /= rhs.xp; xm *= rhs.xm;
  yp /= rhs.yp; ym *= rhs.ym;
  zp /= rhs.zp; zm *= rhs.zm;

  x2p /= rhs.x2p; x2m *= rhs.x2m;
  y2p /= rhs.y2p; y2m *= rhs.y2m;
  z2p /= rhs.z2p; z2m *= rhs.z2m;

  return(*this);
}

bstencil & bstencil::operator/=(const BoutReal rhs)
{
  cc /= rhs;
  
  xp /= rhs; xm /= rhs;
  yp /= rhs; ym /= rhs;
  zp /= rhs; zm /= rhs;

  x2p /= rhs; x2m /= rhs;
  y2p /= rhs; y2m /= rhs;
  z2p /= rhs; z2m /= rhs;
  
  return(*this);
}

bstencil & bstencil::operator^=(const bstencil &rhs)
{
  cc = pow(cc, rhs.cc);
  
  xp = pow(xp, rhs.xp); xm = pow(xm, rhs.xm);
  yp = pow(yp, rhs.yp); ym = pow(ym, rhs.ym);
  zp = pow(zp, rhs.zp); zm = pow(zm, rhs.zm);

  x2p = pow(x2p, rhs.x2p); x2m = pow(x2m, rhs.x2m);
  y2p = pow(y2p, rhs.y2p); y2m = pow(y2m, rhs.y2m);
  z2p = pow(z2p, rhs.z2p); z2m = pow(z2m, rhs.z2m);

  return(*this);
}

bstencil & bstencil::operator^=(const BoutReal rhs)
{
  cc = pow(cc, rhs);
  
  xp = pow(xp, rhs); xm = pow(xm, rhs);
  yp = pow(yp, rhs); ym = pow(ym, rhs);
  zp = pow(zp, rhs); zm = pow(zm, rhs);

  x2p = pow(x2p, rhs); x2m = pow(x2m, rhs);
  y2p = pow(y2p, rhs); y2m = pow(y2m, rhs);
  z2p = pow(z2p, rhs); z2m = pow(z2m, rhs);
  
  return(*this);
}

//////////////// Binary operators ///////////////////////

const bstencil bstencil::operator+(const bstencil &other) const
{
  bstencil result = *this;
  result += other;
  return(result);
}

const bstencil bstencil::operator-(const bstencil &other) const
{
  bstencil result = *this;
  result -= other;
  return(result);
}

const bstencil bstencil::operator*(const bstencil &other) const
{
  bstencil result = *this;
  result *= other;
  return(result);
}

const bstencil bstencil::operator*(const BoutReal rhs) const
{
  bstencil result = *this;
  result *= rhs;
  return(result);
}
const bstencil bstencil::operator/(const bstencil &other) const
{
  bstencil result = *this;
  result /= other;
  return(result);
}

const bstencil bstencil::operator/(const BoutReal rhs) const
{
  bstencil result = *this;
  result /= rhs;
  return(result);
}

const bstencil bstencil::operator^(const bstencil &other) const
{
  bstencil result = *this;
  result ^= other;
  return(result);
}

const bstencil bstencil::operator^(const BoutReal rhs) const
{
  bstencil result = *this;
  result ^= rhs;
  return(result);
}

/////////////// NON-MEMBER OPERATORS /////////////////////

const bstencil operator*(const BoutReal lhs, const bstencil &rhs)
{
  return(rhs * lhs);
}

const bstencil operator/(const BoutReal lhs, const bstencil &rhs)
{
  bstencil result = rhs;
  result.cc = lhs / rhs.cc;
  
  result.xp = lhs / rhs.xp; result.xm = lhs / rhs.xm;
  result.yp = lhs / rhs.yp; result.ym = lhs / rhs.ym;
  result.zp = lhs / rhs.zp; result.zm = lhs / rhs.zm;
  
  result.x2p = lhs / rhs.x2p; result.x2m = lhs / rhs.x2m;
  result.y2p = lhs / rhs.y2p; result.y2m = lhs / rhs.y2m;
  result.z2p = lhs / rhs.z2p; result.z2m = lhs / rhs.z2m;

  return(result);
}

const bstencil operator^(const BoutReal lhs, const bstencil &rhs)
{
  bstencil result = rhs;

  result.cc = pow(lhs, rhs.cc);

  result.xp = pow(lhs, rhs.xp); result.xm = pow(lhs, rhs.xm);
  result.yp = pow(lhs, rhs.yp); result.ym = pow(lhs, rhs.ym);
  result.zp = pow(lhs, rhs.zp); result.zm = pow(lhs, rhs.zm);

  result.x2p = pow(lhs, rhs.x2p); result.x2m = pow(lhs, rhs.x2m);
  result.y2p = pow(lhs, rhs.y2p); result.y2m = pow(lhs, rhs.y2m);
  result.z2p = pow(lhs, rhs.z2p); result.z2m = pow(lhs, rhs.z2m);

  return(result);
}

/*******************************************************************************
 * Functions to handle indexing
 *******************************************************************************/

/* Given index centre, calculates shifted indices */
void calc_index(bindex *bx)
{
  bx->jxp = bx->jx+1;
  if(bx->jxp >= mesh->ngx)
    bx->jxp = mesh->ngx-1;
  
  bx->jxm = bx->jx-1;
  if(bx->jxm < 0)
    bx->jxm = 0;

  if (bx->jx>1) bx->jx2m=bx->jx-2; else bx->jx2m=bx->jxm;
  if (bx->jx<(mesh->ngx-2)) bx->jx2p=bx->jx+2; else bx->jx2p=bx->jxp;

  bx->jyp  = bx->jy+1;
  bx->jym  = bx->jy-1;
  if (bx->jy<mesh->yend || mesh->ystart>1) bx->jy2p = bx->jy+2; else bx->jy2p=bx->jy+1;
  if (bx->jy>mesh->ystart || mesh->ystart>1) bx->jy2m = bx->jy-2; else bx->jy2m=bx->jy-1;

  int ncz = mesh->ngz-1;
  
  bx->jzp  = (bx->jz+1)%ncz;
  bx->jzm  = (bx->jz+ncz-1)%ncz;
  bx->jz2p = (bx->jzp+1)%ncz;
  bx->jz2m = (bx->jzm+ncz-1)%ncz;

  bx->yp_shift = bx->y2p_shift = false;
  bx->ym_shift = bx->y2m_shift = false;

  /* Twist-Shift boundary condition */
  /*
  if(TwistShift && (mesh->TwistOrder != 0)) {
    if( (TS_down_in  && (DDATA_INDEST  != -1) && (bx->jx <  DDATA_XSPLIT)) ||
	(TS_down_out && (DDATA_OUTDEST != -1) && (bx->jx >= DDATA_XSPLIT)) ) {
      
      bx->ym_offset = -ShiftAngle[bx->jx] / mesh->dz;
      if (bx->jy==mesh->ystart) {
	bx->ym_shift = bx->y2m_shift = true;
	
      }else if (bx->jy==mesh->ystart+1) {
	bx->y2m_shift = true;
      }
    }
    if( (TS_up_in  && (UDATA_INDEST  != -1) && (bx->jx <  UDATA_XSPLIT)) ||
	(TS_up_out && (UDATA_OUTDEST != -1) && (bx->jx >= UDATA_XSPLIT)) ) {
      
      bx->yp_offset = ShiftAngle[bx->jx] / mesh->dz;
      if (bx->jy==mesh->yend) {
	bx->yp_shift = bx->y2p_shift = true;
      } else if (bx->jy==mesh->yend-1) {
	bx->y2m_shift = true;
      }
    }
  }
  */

  /* shifted z-indices for x differencing */
  if(mesh->ShiftXderivs && (mesh->ShiftOrder != 0)) {
    bx->xp_offset = (-mesh->zShift[bx->jxp][bx->jy] + mesh->zShift[bx->jx][bx->jy]) / mesh->dz;
    bx->x2p_offset = (-mesh->zShift[bx->jx2p][bx->jy] + mesh->zShift[bx->jx][bx->jy]) / mesh->dz;
    bx->xm_offset =  (-mesh->zShift[bx->jxm][bx->jy] + mesh->zShift[bx->jx][bx->jy]) / mesh->dz;
    bx->x2m_offset = (-mesh->zShift[bx->jx2m][bx->jy] + mesh->zShift[bx->jx][bx->jy]) / mesh->dz;
  }else {
    bx->xp_offset = 0.0;
    bx->x2p_offset = 0.0;
    bx->xm_offset = 0.0;
    bx->x2m_offset = 0.0;
  }
}

/* Resets the index bx */
void start_index(bindex *bx, REGION region)
{
	// Initialize it to something
  bx->jx = 0;

  if((region == RGN_NOBNDRY) || (region == RGN_NOX))
    bx->jx = mesh->xstart;
  
  bx->jy = mesh->ystart;
  bx->jz = 0;

  bx->region = region;

  calc_index(bx);
}

/* Loops the index over all points. Returns 0 when no more */
int next_index3(bindex *bx)
{
  bx->jz++;
  if(bx->jz >= mesh->ngz-1) {
    
    bx->jz = 0;
    bx->jy++;
    
    if(bx->jy > mesh->yend) {
      bx->jy =mesh->ystart;
      bx->jx++;
      
      if((bx->region == RGN_NOBNDRY) || (bx->region == RGN_NOX)) {
	// Missing out X boundary
	if(bx->jx > mesh->xend) {
	  bx->jx = mesh->xstart;
	  return(0);
	}
      }else {
	// Including X boundary regions
	if(bx->jx >= mesh->ngx) {
	  bx->jx = 0;
	  return(0);
	}
      }
    }
  }
  
  calc_index(bx);

  return(1);
}

/* Loops the index over 2D domain (no Z) */
int next_index2(bindex *bx)
{
  bx->jy++;
    
  if(bx->jy > mesh->yend) {
    bx->jy =mesh->ystart;
    bx->jx++;
    
    if((bx->region == RGN_NOBNDRY) || (bx->region == RGN_NOX)) {
      if(bx->jx > mesh->xend) {
	bx->jx = mesh->xstart;
	return(0);
      }
    }else
      if(bx->jx >= mesh->ngx) {
	bx->jx = 0;
	return(0);
      }
  }
  
  calc_index(bx);

  return(1);
}

/* Loops over all perpendicular indices (no Y) */
int next_indexperp(bindex *bx)
{
  bx->jz++;
  if(bx->jz >= mesh->ngz-1) {
    
    bx->jz = 0;
    bx->jx++;
      
    if(bx->jx > mesh->xend) {
      bx->jx = mesh->xstart;
      return(0);
    }
  }
  
  calc_index(bx);

  return(1);
}

/* Resets the index bx to the end of the region*/
void reverse_start_index(bindex *bx, REGION region)
{
	// Initialize it to something
  bx->jx = mesh->ngx-1;

  if((region == RGN_NOBNDRY) || (region == RGN_NOX))
    bx->jx = mesh->xend;
  
  bx->jy = mesh->yend;
  bx->jz = mesh->ngz-2;

  bx->region = region;

  calc_index(bx);
}

/* Resets the index bx to the first x and z but last y values */
void start_index_lasty(bindex *bx, REGION region)
{
	// Initialize it to something
  bx->jx = 0;

  if((region == RGN_NOBNDRY) || (region == RGN_NOX))
    bx->jx = mesh->xstart;
  
  bx->jy = mesh->yend;
  bx->jz = 0;

  bx->region = region;

  calc_index(bx);
}

/* Loops the index backwards over all points. Returns 0 when no more */
int reverse_next_index3(bindex *bx)
{
  bx->jz--;
  if(bx->jz < 0) {
    
    bx->jz = mesh->ngz-2;
    bx->jy--;
    
    if(bx->jy < mesh->ystart) {
      bx->jy =mesh->yend;
      bx->jx--;
      
      if((bx->region == RGN_NOBNDRY) || (bx->region == RGN_NOX)) {
	// Missing out X boundary
	if(bx->jx < mesh->xstart) {
	  bx->jx = mesh->xend;
	  return(0);
	}
      }else {
	// Including X boundary regions
	if(bx->jx < 0) {
	  bx->jx = mesh->ngx-1;
	  return(0);
	}
      }
    }
  }
  
  calc_index(bx);

  return(1);
}

/* Loops the index along y points. Returns 0 when no more */
int next_index_y(bindex *bx)
{
  bx->jy++;
  if(bx->jy > mesh->yend) {
    bx->jy--;
    return(0);
  }
  
  calc_index(bx);

  return(1);
}

/* Loops the index backwards along y points. Returns 0 when no more */
int previous_index_y(bindex *bx)
{
  bx->jy--;
  if(bx->jy < mesh->ystart) {
    bx->jy++;
    return(0);
  }
  
  calc_index(bx);

  return(1);
}
