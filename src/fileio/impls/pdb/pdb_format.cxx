/**************************************************************************
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

#include "pdb_format.hxx"

#ifdef PDBF

#include <utils.hxx>

static char* REALSTR = (char*)"double";
static char* LOWPRECSTR = (char*)"float";

PdbFormat::PdbFormat()
{
  fp = NULL;
  x0 = y0 = z0 = t0 = 0;
  lowPrecision = false;
}

PdbFormat::PdbFormat(const char *name)
{
  fp = NULL;
  x0 = y0 = z0 = t0 = 0;
  lowPrecision = false;
  openr(name);
}

PdbFormat::PdbFormat(const string &name)
{
  fp = NULL;
  x0 = y0 = z0 = t0 = 0;
  lowPrecision = false;
  openr(name);
}

PdbFormat::~PdbFormat()
{
  close();
}


bool PdbFormat::openr(const string &name)
{
  return openr(name.c_str());
}

bool PdbFormat::openr(const char *name)
{
  if(fp != NULL) // Already open. Close then re-open
    close(); 

  if((fp = PD_open((char*) name, (char*)"r")) == NULL)
    return false;

  fname = copy_string(name);

  return true;
}

bool PdbFormat::openw(const string &name, bool append)
{
  return openw(name.c_str(), append);
}

bool PdbFormat::openw(const char *name, bool append)
{
  if(fp != NULL) // Already open. Close then re-open
    close(); 
  
  appending = append;

  if(append) {
    if((fp = PD_open((char*) name, (char*)"a")) == NULL)
      return false;
  }else
    if((fp = PD_open((char*) name, (char*)"w")) == NULL)
      return false;

  fname = copy_string(name);

  return true;
}

bool PdbFormat::is_valid() {
  return (fp != NULL);
}

void PdbFormat::close() {
  if(fp == NULL)
    return;

  PD_close(fp);

  fp = NULL;
}

void PdbFormat::flush() {
  if(!is_valid())
    return;
  
  PD_flush(fp);
}

const char* PdbFormat::filename()
{
  return fname;
}

const vector<int> PdbFormat::getSize(const char *var)
{
  vector<int> size;
  
  syment* ep = PD_query_entry(fp, (char*) var, NULL);

  if(ep == NULL)
    return size;

  dimdes* dims = PD_entry_dimensions(ep);

  if(dims == NULL) {
    size.push_back(1);
    return size;
  }

  while(dims != NULL) {
    size.push_back(dims->index_max - dims->index_min + 1);
    
    dims = dims->next;
  }

  return size;
}

const vector<int> PdbFormat::getSize(const string &var)
{
  return getSize(var.c_str());
}

bool PdbFormat::setGlobalOrigin(int x, int y, int z)
{
  x0 = x;
  y0 = y;
  z0 = z;
  
  return true;
}

bool PdbFormat::setRecord(int t)
{
  t0 = t;

  return true;
}

bool PdbFormat::read(int *var, const char *name, int lx, int ly, int lz)
{
  long inds[9];

  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  // Indices in the format [min, max, step, ...]
  inds[0] = x0; inds[1] = x0 + lx - 1; inds[2] = 1L;
  inds[3] = y0; inds[4] = y0 + ly - 1; inds[5] = 1L;
  inds[6] = z0; inds[7] = z0 + lz - 1; inds[8] = 1L;

  if(PD_read_as_alt(fp, (char*) name, (char*)"integer", var, inds) == 0) {
    return false;
  }
  
  return true;
}

bool PdbFormat::read(int *var, const string &name, int lx, int ly, int lz)
{
  return read(var, name.c_str(), lx, ly, lz);
}

bool PdbFormat::read(BoutReal *var, const char *name, int lx, int ly, int lz)
{
  long inds[9];

  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  // Indices in the format [min, max, step, ...]
  inds[0] = x0; inds[1] = x0 + lx - 1; inds[2] = 1L;
  inds[3] = y0; inds[4] = y0 + ly - 1; inds[5] = 1L;
  inds[6] = z0; inds[7] = z0 + lz - 1; inds[8] = 1L;

  if(PD_read_as_alt(fp, (char*) name, REALSTR, var, inds) == 0) {
    return false;
  }
  
  return true;
}

bool PdbFormat::read(BoutReal *var, const string &name, int lx, int ly, int lz)
{
  return read(var, name.c_str(), lx, ly, lz);
}

bool PdbFormat::write(int *var, const char *name, int lx, int ly, int lz)
{
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  if((lx == 0) && (ly == 0) && (lz == 0)) {
    // Writing a scalar integer

    return (PD_write(fp, (char*) name, (char*)"integer", var) == TRUE);
  }

  int nd = 1; // Number of dimensions
  if(ly != 0) nd = 2;
  if(lz != 0) nd = 3;

  long inds[9];
  
  // Indices in the format [min, max, step, ...]
  inds[0] = x0; inds[1] = x0 + lx - 1; inds[2] = 1L;
  inds[3] = y0; inds[4] = y0 + ly - 1; inds[5] = 1L;
  inds[6] = z0; inds[7] = z0 + lz - 1; inds[8] = 1L;

  return (PD_write_alt(fp, (char*) name, (char*)"integer", var, nd, inds) == TRUE);
}

bool PdbFormat::write(int *var, const string &name, int lx, int ly, int lz)
{
  return write(var, name.c_str(), lx, ly, lz);
}

bool PdbFormat::write(BoutReal *var, const char *name, int lx, int ly, int lz)
{
  char *vartype = REALSTR;
  if(lowPrecision)
    vartype = LOWPRECSTR;

  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  if((lx == 0) && (ly == 0) && (lz == 0)) {
    // Writing a scalar integer

    return (PD_write_as(fp, (char*) name, REALSTR, vartype, var) == TRUE);
  }

  int nd = 1; // Number of dimensions
  if(ly != 0) nd = 2;
  if(lz != 0) nd = 3;

  long inds[9];
  
  // Indices in the format [min, max, step, ...]
  inds[0] = x0; inds[1] = x0 + lx - 1; inds[2] = 1L;
  inds[3] = y0; inds[4] = y0 + ly - 1; inds[5] = 1L;
  inds[6] = z0; inds[7] = z0 + lz - 1; inds[8] = 1L;

  return (PD_write_as_alt(fp, (char*) name, REALSTR, vartype, var, nd, inds) == TRUE);
}

bool PdbFormat::write(BoutReal *var, const string &name, int lx, int ly, int lz)
{
  return write(var, name.c_str(), lx, ly, lz);
}

/***************************************************************************
 * Record-based (time-dependent) data
 ***************************************************************************/

bool PdbFormat::read_rec(int *var, const char *name, int lx, int ly, int lz)
{
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  int nt, t;
  
  if((nt = nrecs(name)) < 1) // Number of time-points 
    return false;

  t = t0;
  if((t < 0) || (t > (nt-1)))
    t = nt-1; // Reading the last time-point
  
  long inds[12];
  
  // Indices in the format [min, max, step, ...]
  inds[0] = t;  inds[1] = t;           inds[2] = 1L; // Just reading one time-slice
  inds[3] = x0; inds[4] = x0 + lx - 1; inds[5] = 1L;
  inds[6] = y0; inds[7] = y0 + ly - 1; inds[8] = 1L;
  inds[9] = z0; inds[10] = z0 + lz - 1; inds[11] = 1L;

  if(PD_read_as_alt(fp, (char*) name, (char*)"integer", var, inds) == 0) {
    return false;
  }
  
  return true;
}

bool PdbFormat::read_rec(int *var, const string &name, int lx, int ly, int lz)
{
  return read_rec(var, name.c_str(), lx, ly, lz);
}

bool PdbFormat::read_rec(BoutReal *var, const char *name, int lx, int ly, int lz)
{
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  int nt, t;
  
  if((nt = nrecs(name)) < 1) // Number of time-points 
    return false;

  t = t0;
  if((t < 0) || (t > (nt-1)))
    t = nt-1; // Reading the last time-point
  
  long inds[12];
  
  // Indices in the format [min, max, step, ...]
  inds[0] = t;  inds[1] = t;           inds[2] = 1L; // Just reading one time-slice
  inds[3] = x0; inds[4] = x0 + lx - 1; inds[5] = 1L;
  inds[6] = y0; inds[7] = y0 + ly - 1; inds[8] = 1L;
  inds[9] = z0; inds[10] = z0 + lz - 1; inds[11] = 1L;

  if(PD_read_as_alt(fp, (char*) name, REALSTR, var, inds) == 0) {
    return false;
  }
  
  return true;
}

bool PdbFormat::read_rec(BoutReal *var, const string &name, int lx, int ly, int lz)
{
  return read_rec(var, name.c_str(), lx, ly, lz);
}

bool PdbFormat::write_rec(int *var, const char *name, int lx, int ly, int lz)
{
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  int t = 0;
  long inds[12];
  
  // Indices in the format [min, max, step, ...]
  inds[0] = t;  inds[1] = t;           inds[2] = 1L;
  inds[3] = x0; inds[4] = x0 + lx - 1; inds[5] = 1L;
  inds[6] = y0; inds[7] = y0 + ly - 1; inds[8] = 1L;
  inds[9] = z0; inds[10] = z0 + lz - 1; inds[11] = 1L;

  int nd = 1; // Number of dimensions
  if(lx != 0) nd = 2;
  if(ly != 0) nd = 3;
  if(lz != 0) nd = 4;

  if(appending) {
    return (PD_append_alt(fp, (char*) name, var, nd, inds) == TRUE);
  }
  
  return (PD_write_alt(fp, (char*) name, (char*)"integer", var, nd, inds) == TRUE);
}

bool PdbFormat::write_rec(int *var, const string &name, int lx, int ly, int lz)
{
  return write_rec(var, name.c_str(), lx, ly, lz);
}

bool PdbFormat::write_rec(BoutReal *var, const char *name, int lx, int ly, int lz)
{
  char *vartype = REALSTR;
  if(lowPrecision)
    vartype = LOWPRECSTR;
  
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  int t = 0;
  long inds[12];
  
  // Indices in the format [min, max, step, ...]
  inds[0] = t;  inds[1] = t;           inds[2] = 1L;
  inds[3] = x0; inds[4] = x0 + lx - 1; inds[5] = 1L;
  inds[6] = y0; inds[7] = y0 + ly - 1; inds[8] = 1L;
  inds[9] = z0; inds[10] = z0 + lz - 1; inds[11] = 1L;

  int nd = 1; // Number of dimensions
  if(lx != 0) nd = 2;
  if(ly != 0) nd = 3;
  if(lz != 0) nd = 4;

  if(appending) {
    return (PD_append_as_alt(fp, (char*) name, REALSTR, var, nd, inds) == TRUE);
  }
  
  return (PD_write_as_alt(fp, (char*) name, REALSTR, vartype, var, nd, inds) == TRUE);
}

bool PdbFormat::write_rec(BoutReal *var, const string &name, int lx, int ly, int lz)
{
  return write_rec(var, name.c_str(), lx, ly, lz);
}

/***************************************************************************
 * Private functions
 ***************************************************************************/

int PdbFormat::nrecs(const char *name)
{
  syment *ep;
  char *type;
  long size;
  int nd;
  long *dims;
  int nt;
  
  if((ep = PD_query_entry(fp, (char*) name, NULL)) == NULL)
    return 0;
  
  if(!PD_get_entry_info(ep, &type, &size, &nd, &dims))
    return 0;
  
  if(nd < 1) {
    nt = 0;
  }else
    nt = dims[1] - dims[0] + 1;
  
  //PD_free_entry_info(type, dims);
  
  return nt;
}

#endif // PDBF

