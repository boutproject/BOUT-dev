#include <boutmain.hxx>
#include <derivs.hxx>
#include <field_factory.hxx>
#include <bout/mesh.hxx>
#include <assert.h>

#include <interpolation.hxx>

Field3D n;

void save_diff(Field3D a, Field3D b,CELL_LOC in,CELL_LOC out=CELL_CENTRE);


#define test(x,xl) {                                                    \
    if (Options::getRoot()->                                            \
        getSection("diff_" #x)->isSet("function")){                     \
      for (auto in : {CELL_CENTRE, CELL_ ##x ## LOW} ){                 \
        auto inv = ff->create3D("v:func",oo,mesh,in,0);                 \
        for (auto out : {CELL_CENTRE, CELL_ ##x ## LOW} ){              \
          auto inf = ff->create3D("f:func",oo,mesh,out,0);              \
          save_diff(/*mesh->indexVDD##x(inv,inf,out,DIFF_DEFAULT)       \
                      /mesh->coordinates()->d ##xl                  */  \
                    VDD##x(inv,inf,out)                                 \
                    ,                                                   \
                    ff->create3D(                                       \
                                 "diff_" #x ":function",                \
                                 oo,                                    \
                                 mesh, out, 0),in,out);                 \
        }                                                               \
      }                                                                 \
    }                                                                   \
  }


Field3D error;
Field3D diff;
Field3D exact;
int physics_init(bool restart) {
  SOLVE_FOR(n);
  auto ff = FieldFactory::get();
  auto oo = Options::getRoot();
  n=0;
  error=0;
  SAVE_REPEAT(error);
  SAVE_REPEAT(diff);
  SAVE_REPEAT(exact);
  test(X,x);
  
  test(Y,y);

  test(Z,z);
  
  return 0;
}


int physics_run(BoutReal time) {
  ddt(n)=0;
  return 0;
}
void save_diff(Field3D a, Field3D b,CELL_LOC in,CELL_LOC out){
  exact=b;
  diff=a;
  error= a-b;
  
  checkData(error);
  // std::string name="error_";
  // name+=in;//strLocation(in);
  // name+="_";
  // name+=out;//strLocation(out);
  //dump.addOnce(error,name.c_str());
  dump.write();
}

void compare( Field3D diff, Field3D exp){
  static int counter=0;
  int print=-1;
  //OPTION(Options::getRoot(),print,-1);
  if (counter++==print){
    //for (int x=mesh->xstart;x<=mesh->xend;++x)
      for (int x=0;x<mesh->LocalNx;++x)
      for (int y=mesh->ystart;y<=mesh->yend;++y)
	//for (int y=0;y<mesh->LocalNy;++y)
	for (int z=0;z<mesh->LocalNz;++z)
	  output.write("\t%2d %2d %2d   %8g\t%8g\t%8g\n",x,y,z,exp(x,y,z)
		       ,diff(x,y,z),n(x,y,z));
  
  }
  auto error=max(abs(diff-exp)); 
  output.write("\nerror: %g\n",error);
  if (error > 1e-8){
    //int x=mesh->xend-1;
    //int y=mesh->ystart;
    //int z=mesh->LocalNz-1;
    auto err=exp-diff;
    
    //for (int x=mesh->xstart;x<=mesh->xend;++x)
    for (int x=0;x<mesh->LocalNx;++x)
      //for (int y=mesh->ystart;y<=mesh->yend;++y)
      for (int y=0;y<mesh->LocalNy;++y)
	for (int z=0;z<mesh->LocalNz;++z)
	  output.write("\t%2d %2d %2d   %8g\t%8g\t%8g\n",x,y,z,exp(x,y,z)
		       ,diff(x,y,z),n(x,y,z));
    
    throw BoutException("Error is to large: %g",error);
  }
}

