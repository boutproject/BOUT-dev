/************************************************************************
 * Inversion of parallel derivatives
 * Intended for use in preconditioner for reduced MHD
 * 
 * Inverts a matrix of the form 
 *
 * A + B * Grad2_par2 + C*D2DYDZ + D*D2DZ2 + E*DDY
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
 ************************************************************************/

#ifndef __INV_PAR_H__
#define __INV_PAR_H__

#include "field3d.hxx"
#include "field2d.hxx"
#include "options.hxx"

// Parderiv implementations
#define PARDERIVSERIAL "serial"
#define PARDERIVCYCLIC "cyclic"

/// Base class for parallel inversion solvers
class InvertPar {
public:
  InvertPar(Options *opt) {}
  virtual ~InvertPar() {}
  
  static InvertPar* Create();
  
  virtual const Field2D solve(const Field2D &f); ///< Warning: Default implementation very inefficient
  virtual const Field3D solve(const Field3D &f) = 0;  ///< This method must be implemented
  
  virtual const Field3D solve(const Field2D &f, const Field2D &start) {return solve(f);}
  virtual const Field3D solve(const Field3D &f, const Field3D &start) {return solve(f);}
  
  virtual void setCoefA(const Field2D &f) = 0;
  virtual void setCoefA(const Field3D &f) {setCoefA(f.DC());}
  virtual void setCoefA(const BoutReal f) {setCoefA(Field2D(f));}
  
  virtual void setCoefB(const Field2D &f) = 0;
  virtual void setCoefB(const Field3D &f) {setCoefB(f.DC());}
  virtual void setCoefB(const BoutReal f) {setCoefB(Field2D(f));}
  
  virtual void setCoefC(const Field2D &f) = 0;
  virtual void setCoefC(const Field3D &f) {setCoefB(f.DC());}
  virtual void setCoefC(const BoutReal f) {setCoefB(Field2D(f));}
  
  virtual void setCoefD(const Field2D &f) = 0;
  virtual void setCoefD(const Field3D &f) {setCoefB(f.DC());}
  virtual void setCoefD(const BoutReal f) {setCoefB(Field2D(f));}
  
  virtual void setCoefE(const Field2D &f) = 0;
  virtual void setCoefE(const Field3D &f) {setCoefB(f.DC());}
  virtual void setCoefE(const BoutReal f) {setCoefB(Field2D(f));}
  
private:
};


#endif // __INV_PAR_H__

