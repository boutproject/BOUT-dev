/*
 * Global fields for gather/scatter
 * 
 * 
 */

#include <bout.hxx>
#include <boutmain.hxx>
#include <bout/globalfield.hxx>

int physics_init(bool UNUSED(restarting)) {
  
  /////////////////////////////////////////////////////////////
  // 2D fields

  // Create local variables, fill with data
  Field2D localX, localY;

  localX.allocate();
  localY.allocate();
  
  for(int x=0;x<mesh->LocalNx;x++) {
    for(int y=0;y<mesh->LocalNy;y++) {
      localX(x,y) = mesh->getGlobalXIndex(x);
      localY(x,y) = mesh->getGlobalYIndex(y - mesh->ystart);
    }
  }
  
  // Gather onto one processor (0 by default)
  GlobalField2D gX(mesh), gY(mesh);

  gX.gather(localX);
  gY.gather(localY);
  
  if(gX.dataIsLocal()) {
    // Data is on this processor
    bool gather_pass = true;
    for(int x=0;x<gX.xSize();x++)
      for(int y=0;y<gX.ySize();y++) {
        if( (ROUND(gX(x,y)) != x) || (ROUND(gY(x,y)) != y) ) {
          output.write("{:d}, {:d} :  {:e}, {:e}\n", x,y, gX(x,y), gY(x,y));
          gather_pass = false;
        }
      }
    output << "2D GATHER TEST: " << gather_pass << endl;
  }
  
  // Scatter back and check
  Field2D scatX = gX.scatter();
  Field2D scatY = gY.scatter();
  
  bool scatter_pass = true;
  for(int x=mesh->xstart;x<=mesh->xend;x++)
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      if( (localX(x,y) != scatX(x,y)) || (localY(x,y) != scatY(x,y)) ) {
        output.write("{:d}, {:d} :  ({:e}, {:e}) ({:e}, {:e})", x, y,
                   localX(x,y), localY(x,y), scatX(x,y), scatY(x,y));
        scatter_pass = false;
      }
    }
  output << "2D SCATTER TEST: " << scatter_pass << endl;
  
  /////////////////////////////////////////////////////////////
  // 3D fields

  
  // Create local variables, fill with data
  Field3D localX3D, localY3D;

  localX3D.allocate();
  localY3D.allocate();
  
  for(int x=0;x<mesh->LocalNx;x++)
    for(int y=0;y<mesh->LocalNy;y++)
      for(int z=0;z<mesh->LocalNz;z++) {
        localX3D(x,y,z) = mesh->getGlobalXIndex(x) + z;
        localY3D(x,y,z) = mesh->getGlobalYIndex(y - mesh->ystart) + z;
      }
  
  // Gather onto one processor (0 by default)
  GlobalField3D gX3D(mesh), gY3D(mesh);

  gX3D.gather(localX3D);
  gY3D.gather(localY3D);
  
  if(gX3D.dataIsLocal()) {
    // Data is on this processor
    bool gather_pass3D = true;
    for(int x=0;x<gX3D.xSize();x++)
      for(int y=0;y<gX3D.ySize();y++) 
        for(int z=0;z<gX3D.zSize();z++) {
          if( (ROUND(gX3D(x,y,z)) != x + z) || (ROUND(gY3D(x,y,z)) != y + z) ) {
            output.write("{:d}, {:d}, {:d} :  {:e}, {:e}\n", x,y,z, gX3D(x,y,z), gY3D(x,y,z));
            gather_pass3D = false;
          }
      }
    output << "3D GATHER TEST: " << gather_pass3D << endl;
  }
  
  // Scatter back and check
  Field3D scatX3D = gX3D.scatter();
  Field3D scatY3D = gY3D.scatter();
  
  bool scatter_pass3D = true;
  for(int x=mesh->xstart;x<=mesh->xend;x++)
    for(int y=mesh->ystart;y<=mesh->yend;y++)
      for(int z=0;z<mesh->LocalNz;z++) {
        if( (localX3D(x,y,z) != scatX3D(x,y,z)) || (localY3D(x,y,z) != scatY3D(x,y,z)) ) {
          output.write("{:d}, {:d}, {:d} :  ({:e}, {:e}) ({:e}, {:e})", x, y, z,
                       localX3D(x,y,z), localY3D(x,y,z), scatX3D(x,y,z), scatY3D(x,y,z));
          scatter_pass3D = false;
        }
      }
  output << "2D SCATTER TEST: " << scatter_pass3D << endl;


  return 1; // Signal an error, so quits
}

int physics_run(BoutReal UNUSED(t)) {
  // Doesn't do anything
  return 1;
}
