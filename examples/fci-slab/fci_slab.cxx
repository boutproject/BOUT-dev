#include <bout/physicsmodel.hxx>
#include <derivs.hxx>
#include <utils.hxx>

class FCIMap {
public:
  FCIMap(Field3D xt_prime, Field3D zt_prime, int nx, int ny, int nz);
  FCIMap(BoutReal* xt_prime, BoutReal* zt_prime, int nx, int ny, int nz);
  // ~FCIMap();

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

FCIMap::FCIMap(Field3D xt_prime, Field3D zt_prime, int nx, int ny, int nz) {

  BoutReal t_x, t_z, temp;

  i_corner = i3tensor(nx, ny, nz);
  k_corner = i3tensor(nx, ny, nz);

  for(int x=0;x<nx;x++) {
	for(int y=0; y<ny;y++) {
	  for(int z=0;z<nz;z++) {
		i_corner[x][y][z] = (int)xt_prime[x][y][z];
		k_corner[x][y][z] = (int)zt_prime[x][y][z];

		t_x = xt_prime[x][y][z] - (BoutReal)i_corner[x][y][z];
		t_z = zt_prime[x][y][z] - (BoutReal)k_corner[x][y][z];

		std::cout << k_corner[x][y][z] << " " << zt_prime[x][y][z] << " " << t_z << " " << 2.*t_z*t_z*t_z - 3.*t_z*t_z + 1. << std::endl;

		temp = 2.*t_x*t_x*t_x - 3.*t_x*t_x + 1.;
        a_x.setData(x, y, z, &temp);
		temp = 2.*t_z*t_z*t_z - 3.*t_z*t_z + 1.;
        a_z.setData(x, y, z, &temp);

		temp = -2.*t_x*t_x*t_x + 3.*t_x*t_x;
        a_1mx.setData(x, y, z, &temp);
		temp = -2.*t_z*t_z*t_z + 3.*t_z*t_z;
        a_1mz.setData(x, y, z, &temp);

		temp = t_x*(1.-t_x)*(1.-t_x);
        b_x.setData(x, y, z, &temp);
		temp = t_z*(1.-t_z)*(1.-t_z);
		b_z.setData(x, y, z, &temp);

		temp = t_x*t_x*t_x - t_x*t_x;
		b_1mx.setData(x, y, z, &temp);
		temp = t_z*t_z*t_z - t_z*t_z;
		b_1mz.setData(x, y, z, &temp);

	  }
	}
  }

}

FCIMap::FCIMap(BoutReal* xt_prime, BoutReal* zt_prime, int nx, int ny, int nz) {

  BoutReal t_x, t_z, temp;

  i_corner = i3tensor(nx, ny, nz);
  k_corner = i3tensor(nx, ny, nz);

  for(int x=0;x<nx;x++) {
	for(int y=0; y<ny;y++) {
	  for(int z=0;z<nz;z++) {
		i_corner[x][y][z] = (int)xt_prime[x + y*nx + z*ny*nx];
		k_corner[x][y][z] = (int)zt_prime[x + y*nx + z*ny*nx];

		t_x = xt_prime[x + y*nx + z*ny*nx] - (BoutReal)i_corner[x][y][z];
		t_z = zt_prime[x + y*nx + z*ny*nx] - (BoutReal)k_corner[x][y][z];

		temp = 2.*t_x*t_x*t_x - 3.*t_x*t_x + 1.;
        a_x.setData(x, y, z, &temp);
		temp = 2.*t_z*t_z*t_z - 3.*t_z*t_z + 1.;
        a_z.setData(x, y, z, &temp);

		temp = -2.*t_x*t_x*t_x + 3.*t_x*t_x;
        a_1mx.setData(x, y, z, &temp);
		temp = -2.*t_z*t_z*t_z + 3.*t_z*t_z;
        a_1mz.setData(x, y, z, &temp);

		temp = t_x*(1.-t_x)*(1.-t_x);
        b_x.setData(x, y, z, &temp);
		temp = t_z*(1.-t_z)*(1.-t_z);
		b_z.setData(x, y, z, &temp);

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
	int init(bool restarting) {
		const string fwdfilename = "forward_coefs.nc";

		// Create a file format handler
		DataFormat *fwdfile = data_format(fwdfilename.c_str());

		fwdfile->openr(fwdfilename);

		if(!fwdfile->is_valid()) {
			output << "\tERROR: Could not open file " << fwdfilename << endl;
		}

		BoutReal forward_xt[5][5][5];
		BoutReal forward_zt[5][5][5];
		Field3D forward_xt_f;
		Field3D forward_zt_f;

		const string xt_prime = "xt_prime";
		const string zt_prime = "zt_prime";

		fwdfile->read(**forward_xt, xt_prime, 5, 5, 5);
		fwdfile->read(**forward_zt, zt_prime, 5, 5, 5);

		fwdfile->close();

		const string bkdfilename = "backward_coefs.nc";

		// Create a file format handler
		DataFormat *bkdfile = data_format(bkdfilename.c_str());

		bkdfile->openr(bkdfilename);

		if(!bkdfile->is_valid()) {
			output << "\tERROR: Could not open file " << bkdfilename << endl;
		}

		BoutReal backward_xt[5][5][5];
		BoutReal backward_zt[5][5][5];
		Field3D backward_xt_f;
		Field3D backward_zt_f;

		bkdfile->read(**backward_xt, xt_prime, 5, 5, 5);
		bkdfile->read(**backward_zt, zt_prime, 5, 5, 5);

		bkdfile->close();

		// Turn array of BoutReals into Field3D
		// for (int x=0;x<5;++x) {
		// 	for (int y=0;y<5;++y) {
		// 		for (int z =0;z<5;++z) {
		// 		  forward_xt_f.setData(x, y, z, &forward_xt[x][y][z]);
		// 		  forward_zt_f.setData(x, y, z, &forward_zt[x][y][z]);
		// 		  backward_xt_f.setData(x, y, z, &backward_xt[x][y][z]);
		// 		  backward_zt_f.setData(x, y, z, &backward_zt[x][y][z]);
		// 		}
		// 	}
		// }

		// FCIMap forward_map(forward_xt_f, forward_zt_f, 5, 5, 5);
		// FCIMap backward_map(backward_xt_f, backward_zt_f, 5, 5, 5);

		FCIMap forward_map(**forward_xt, **forward_zt, 5, 5, 5);
		FCIMap backward_map(**backward_xt, **backward_zt, 5, 5, 5);

		// for (int i=0;i<5;++i) {
		// 	for (int j=0;j<5;++j) {
		// 		for (int k =0;k<5;++k) {
		// 			std::cout << forward_zt[i][j][k] << "\t";
		// 		}
		// 	}
		// 	std::cout << "\n";
		// }

		// std::cout << "\n\n";

		// for (int i=0;i<5;++i) {
		// 	for (int j=0;j<5;++j) {
		// 		for (int k =0;k<5;++k) {
		// 			std::cout << forward_map.k_corner[i][j][k] << "\t";
		// 		}
		// 	}
		// 	std::cout << "\n";
		// }



		solver->add(f, "f");
		return 0;
	}
	int rhs(BoutReal time) {
		ddt(f) = f;
		return 0;
	}

	void interpolate_FCI(Field3D &f, Field3D &f_next, FCIMap &fcimap, int dir);
	const Field3D& Grad_par_FCI(const Field3D &f);
private:
	Field3D f;
};

BOUTMAIN(FCISlab);

void FCISlab::interpolate_FCI(Field3D &f, Field3D &f_next, FCIMap &fcimap, int dir) {
    // Field3D *fup = f.yup();

    Field3D fx = DDX(f);
    mesh->communicate(fx);
    Field3D fz = DDZ(f);
    mesh->communicate(fz);
    Field3D fxz = D2DXDZ(f);
    mesh->communicate(fxz);

    for(int x=mesh->xstart;x<=mesh->xend;x++) {
	  for(int y=mesh->ystart; y<=mesh->yend;y++) {
		for(int z=0;z<mesh->ngz-1;z++) {
		  f_next(x,y,z) = (f(fcimap.i_corner[x][y][z], y + dir, fcimap.k_corner[x][y][z])*fcimap.a_x[x][y][z]
						   + f(fcimap.i_corner[x][y][z]+1, y + dir, fcimap.k_corner[x][y][z])*fcimap.a_1mx[x][y][z]
						   + fx( fcimap.i_corner[x][y][z], y + dir, fcimap.k_corner[x][y][z])*fcimap.b_x[x][y][z]
						   - fx( fcimap.i_corner[x][y][z]+1, y + dir, fcimap.k_corner[x][y][z])*fcimap.b_1mx[x][y][z])*fcimap.a_z[x][y][z]
			+ (f( fcimap.i_corner[x][y][z], y + dir, fcimap.k_corner[x][y][z]+1)*fcimap.a_x[x][y][z]
			   + f( fcimap.i_corner[x][y][z]+1, y + dir, fcimap.k_corner[x][y][z]+1)*fcimap.a_1mx[x][y][z]
			   + fx( fcimap.i_corner[x][y][z], y + dir, fcimap.k_corner[x][y][z]+1)*fcimap.b_x[x][y][z]
			   - fx( fcimap.i_corner[x][y][z]+1, y + dir, fcimap.k_corner[x][y][z]+1)*fcimap.b_1mx[x][y][z])*fcimap.a_1mz[x][y][z]
			+ (fz(fcimap.i_corner[x][y][z], y + dir, fcimap.k_corner[x][y][z])*fcimap.a_x[x][y][z]
			   + fz( fcimap.i_corner[x][y][z]+1, y + dir, fcimap.k_corner[x][y][z])*fcimap.a_1mx[x][y][z]
			   + fxz(fcimap.i_corner[x][y][z], y + dir, fcimap.k_corner[x][y][z])*fcimap.b_x[x][y][z]
			   - fxz(fcimap.i_corner[x][y][z]+1, y + dir, fcimap.k_corner[x][y][z])*fcimap.b_1mx[x][y][z])*fcimap.b_z[x][y][z]
			- (fz(fcimap.i_corner[x][y][z], y + dir, fcimap.k_corner[x][y][z]+1)*fcimap.a_x[x][y][z]
			   + fz( fcimap.i_corner[x][y][z]+1, y + dir, fcimap.k_corner[x][y][z]+1)*fcimap.a_1mx[x][y][z]
			   + fxz(fcimap.i_corner[x][y][z], y + dir, fcimap.k_corner[x][y][z]+1)*fcimap.b_x[x][y][z]
			   - fxz(fcimap.i_corner[x][y][z]+1, y + dir, fcimap.k_corner[x][y][z]+1)*fcimap.b_1mx[x][y][z])*fcimap.b_1mz[x][y][z];
		}
	  }
    }
}

const Field3D& FCISlab::Grad_par_FCI(const Field3D &f) {

}
