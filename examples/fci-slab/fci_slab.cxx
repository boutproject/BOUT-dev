#include <bout/physicsmodel.hxx>
#include <derivs.hxx>
#include <utils.hxx>
#include <interpolation.hxx>
#include <bout/assert.hxx>

Field3D fx;
Field3D fz;
Field3D fxz;

class FCIMap {
public:
  FCIMap() {};

  void genCoeffs(BoutReal* xt_prime, BoutReal* zt_prime, int nx, int ny, int nz);
  void genCoeffs(Field3D xt_prime, Field3D zt_prime, int nx, int ny, int nz);

  int*** i_corner;
  int*** k_corner;
  Field3D a_x;
  Field3D a_z;
  Field3D a_1mx;
  Field3D a_1mz;
  Field3D b_x;
  Field3D b_z;
  Field3D b_1mx;
  Field3D b_1mz;
};

void FCIMap::genCoeffs(Field3D xt_prime, Field3D zt_prime, int nx, int ny, int nz) {

  BoutReal t_x, t_z, temp;

  i_corner = i3tensor(nx, ny, nz);
  k_corner = i3tensor(nx, ny, nz);

  int ncz = mesh->ngz-1;

  for(int x=mesh->xstart;x<=mesh->xend;x++) {
	for(int y=mesh->ystart; y<=mesh->yend;y++) {
	  for(int z=0;z<ncz;z++) {
		i_corner[x][y][z] = (int)(xt_prime[x][y][z]);

                
		zt_prime[x][y][z] = zt_prime[x][y][z] - ncz * ( (int) (zt_prime[x][y][z] / ((BoutReal) ncz)) );
		if(zt_prime[x][y][z] < 0.0)
		  zt_prime[x][y][z] += ncz;
                  
		k_corner[x][y][z] = (int)(zt_prime[x][y][z]);

		t_x = xt_prime[x][y][z] - (BoutReal)i_corner[x][y][z];
		t_z = zt_prime[x][y][z] - (BoutReal)k_corner[x][y][z];
                
		// Check that t_x and t_z are in range
		if( (t_x < 0.0) || (t_x > 1.0) )
		  throw BoutException("t_x=%e out of range at (%d,%d,%d)", t_x, x,y,z);
                
		if( (t_z < 0.0) || (t_z > 1.0) )
		  throw BoutException("t_z=%e out of range at (%d,%d,%d)", t_z, x,y,z);
                

		// NOTE: A (small) hack to avoid one-sided differences
		if( i_corner[x][y][z] == mesh->xend ) {
		  i_corner[x][y][z] -= 1;
		  t_x = 1.0;
		}

		// This is h00 in Wikipedia page
		temp = 2.*t_x*t_x*t_x - 3.*t_x*t_x + 1.;
		a_x.setData(x, y, z, &temp);
		temp = 2.*t_z*t_z*t_z - 3.*t_z*t_z + 1.;
		a_z.setData(x, y, z, &temp);
                
		// h01 in Wikipedia page
		temp = -2.*t_x*t_x*t_x + 3.*t_x*t_x;
		a_1mx.setData(x, y, z, &temp);
		temp = -2.*t_z*t_z*t_z + 3.*t_z*t_z;
		a_1mz.setData(x, y, z, &temp);

		// h10
		temp = t_x*(1.-t_x)*(1.-t_x);
		b_x.setData(x, y, z, &temp);
		temp = t_z*(1.-t_z)*(1.-t_z);
		b_z.setData(x, y, z, &temp);

		// h11
		temp = t_x*t_x*t_x - t_x*t_x;
		b_1mx.setData(x, y, z, &temp);
		temp = t_z*t_z*t_z - t_z*t_z;
		b_1mz.setData(x, y, z, &temp);

	  }
	}
  }

}

class FCISlab : public PhysicsModel {
public:

  FCIMap forward_map;
  FCIMap backward_map;

  int init(bool restarting) {

	Field3D forward_xt_prime; 
	Field3D forward_zt_prime; 
	Field3D backward_xt_prime;
	Field3D backward_zt_prime;

	// mesh->dx = 0.1;
	// mesh->dy = 0.6283185307179586;

	// SAVE_ONCE(mesh->dx);
	// SAVE_ONCE(mesh->dy);

	GRID_LOAD(forward_xt_prime); 
	GRID_LOAD(forward_zt_prime); 
	GRID_LOAD(backward_xt_prime);
	GRID_LOAD(backward_zt_prime);

	forward_map.genCoeffs(forward_xt_prime, forward_zt_prime, mesh->ngx, mesh->ngy, mesh->ngz-1);
	backward_map.genCoeffs(backward_xt_prime, backward_zt_prime, mesh->ngx, mesh->ngy, mesh->ngz-1);

	dump.add(*(f.yup()), "yup", 0);
	dump.add(*(f.ydown()), "ydown", 0);

	solver->add(f, "f");
	solver->add(g, "g");

	f.applyBoundary("dirichlet");
	g.applyBoundary("dirichlet");

	SAVE_REPEAT2(ddt(f), ddt(g));

	return 0;
  }
  int rhs(BoutReal time) {
	mesh->communicate(f,g);
	ddt(f) = Grad_par_FCI(g);
	ddt(g) = Grad_par_FCI(f);
	return 0;
  }

  void interpolate_FCI(const Field3D &f, Field3D &f_next, FCIMap &fcimap, int dir);
  const Field3D Grad_par_FCI(Field3D &f);
private:
  Field3D f, g;
};

BOUTMAIN(FCISlab);

void FCISlab::interpolate_FCI(const Field3D &f, Field3D &f_next, FCIMap &fcimap, int dir) {

  if(!mesh->FCI)
    return; // Not using FCI method. Print error / warning?

  fx = DDX(f) * mesh->dx;
  mesh->communicate(fx);
  fz = DDZ(f) * mesh->dz;
  mesh->communicate(fz);
  fxz = D2DXDZ(f) * mesh->dx * mesh->dz;
  mesh->communicate(fxz);

  f_next = 0;

  for(int x=mesh->xstart;x<=mesh->xend;x++) {
	for(int y=mesh->ystart; y<=mesh->yend;y++) {
	  for(int z=0;z<mesh->ngz-1;z++) {

		int ncz = mesh->ngz-1;
		int z_mod = ((fcimap.k_corner[x][y][z] % ncz) + ncz) % ncz;
		int z_mod_p1 = (z_mod + 1) % ncz;

		// Interpolate f in X at Z
		BoutReal f_z = f(fcimap.i_corner[x][y][z], y + dir, z_mod)*fcimap.a_x[x][y][z]
		  + f(fcimap.i_corner[x][y][z]+1, y + dir, z_mod)*fcimap.a_1mx[x][y][z]
		  + fx( fcimap.i_corner[x][y][z], y + dir, z_mod)*fcimap.b_x[x][y][z]
		  + fx( fcimap.i_corner[x][y][z]+1, y + dir, z_mod)*fcimap.b_1mx[x][y][z];
                   
		// Interpolate f in X at Z+1
		BoutReal f_zp1 = f( fcimap.i_corner[x][y][z], y + dir, z_mod_p1)*fcimap.a_x[x][y][z]
		  + f( fcimap.i_corner[x][y][z]+1, y + dir, z_mod_p1)*fcimap.a_1mx[x][y][z]
		  + fx( fcimap.i_corner[x][y][z], y + dir, z_mod_p1)*fcimap.b_x[x][y][z]
		  + fx( fcimap.i_corner[x][y][z]+1, y + dir, z_mod_p1)*fcimap.b_1mx[x][y][z];
                   
		// Interpolate fz in X at Z
		BoutReal fz_z = fz(fcimap.i_corner[x][y][z], y + dir, z_mod)*fcimap.a_x[x][y][z]
		  + fz( fcimap.i_corner[x][y][z]+1, y + dir, z_mod)*fcimap.a_1mx[x][y][z]
		  + fxz(fcimap.i_corner[x][y][z], y + dir, z_mod)*fcimap.b_x[x][y][z]
		  + fxz(fcimap.i_corner[x][y][z]+1, y + dir, z_mod)*fcimap.b_1mx[x][y][z];
                   
		// Interpolate fz in X at Z+1
		BoutReal fz_zp1 = fz(fcimap.i_corner[x][y][z], y + dir, z_mod_p1)*fcimap.a_x[x][y][z]
		  + fz( fcimap.i_corner[x][y][z]+1, y + dir, z_mod_p1)*fcimap.a_1mx[x][y][z]
		  + fxz(fcimap.i_corner[x][y][z], y + dir, z_mod_p1)*fcimap.b_x[x][y][z]
		  + fxz(fcimap.i_corner[x][y][z]+1, y + dir, z_mod_p1)*fcimap.b_1mx[x][y][z];
                  
		// Interpolate in Z
		f_next(x,y + dir,z) = 
		  + f_z    * fcimap.a_z[x][y][z]
		  + f_zp1  * fcimap.a_1mz[x][y][z]
		  + fz_z   * fcimap.b_z[x][y][z]
		  + fz_zp1 * fcimap.b_1mz[x][y][z];
	  }
	}
  }
}

const Field3D FCISlab::Grad_par_FCI(Field3D &f) {

  Field3D result;
  Field3D *yup, *ydown;

  result.allocate();

  yup = f.yup();
  ydown = f.ydown();

  interpolate_FCI(f, *yup, forward_map, +1);
  interpolate_FCI(f, *ydown, backward_map, -1);

  for (int x=mesh->xstart;x<=mesh->xend;++x) {
	for (int y=mesh->ystart;y<=mesh->yend;++y) {
	  for (int z=0;z<mesh->ngz-1;++z) {
		result(x,y,z) = ((*yup)(x,y+1,z) - (*ydown)(x,y-1,z))/(2*mesh->dy(x,y));
	  }
	}
  }
  return result;
}
