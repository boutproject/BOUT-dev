
#include "boundary_region.h"
#include "utils.h"
#include "globals.h"

BoundaryRegionXIn::BoundaryRegionXIn(int ymin, int ymax)
{
  location = BNDRY_XIN;
  x = mesh->xstart-1; // First point inside the boundary
  ys = ymin;
  ye = ymax;
  if(ye < ys)
    SWAP(ys, ye);
 
  // Unit vector out of the domain
  bx = -1;
  by =  0;
}

void BoundaryRegionXIn::first()
{
  x = mesh->xstart-1;
  y = ys;
}

void BoundaryRegionXIn::next()
{
  // Loop over all points, from inside out
  y++;
  if(y > ye) {
    y = ys;
    x--; // Going from inside out
  }
}

void BoundaryRegionXIn::nextX()
{
  x--;
  if(y > ye)
    y = ys;
}

void BoundaryRegionXIn::nextY()
{
  y++;
  if(x < 0)
    x = mesh->xstart-1;
}

bool BoundaryRegionXIn::isDone()
{
  return (x < 0) || (y > ye); // Return true if gone out of the boundary
}

///////////////////////////////////////////////////////////////


BoundaryRegionXOut::BoundaryRegionXOut(int ymin, int ymax)
{
  location = BNDRY_XOUT;
  x = mesh->xend+1; // First point inside the boundary
  ys = ymin;
  ye = ymax;
  if(ye < ys)
    SWAP(ys, ye);
 
  // Unit vector out of the domain
  bx = 1;
  by = 0;
}

void BoundaryRegionXOut::first()
{
  x = mesh->xend+1;
  y = ys;
}

void BoundaryRegionXOut::next()
{
  // Loop over all points, from inside out
  y++;
  if(y > ye) {
    y = ys;
    x++; // Going from inside out
  }
}

void BoundaryRegionXOut::nextX()
{
  x++;
  if(y > ye)
    y = ys;
}

void BoundaryRegionXOut::nextY()
{
  y++;
  if(x >= mesh->ngx)
    x = mesh->xend+1;
}

bool BoundaryRegionXOut::isDone()
{
  return (x >= mesh->ngx) || (y > ye); // Return true if gone out of the boundary
}

///////////////////////////////////////////////////////////////


BoundaryRegionYDown::BoundaryRegionYDown(int xmin, int xmax)
{
  location = BNDRY_XOUT;
  y = mesh->ystart-1; // First point inside the boundary
  xs = xmin;
  xe = xmax;
  
  if(xe < xs)
    SWAP(xs, xe);
 
  // Unit vector out of the domain
  bx =  0;
  by = -1;
}

void BoundaryRegionYDown::first()
{
  x = xs;
  y = mesh->ystart-1;
}

void BoundaryRegionYDown::next()
{
  // Loop over all points, from inside out
  y--;
  if(y < 0) {
    y = mesh->ystart-1;
    x++;
  }
}

void BoundaryRegionYDown::nextX()
{
  x++;
  if(y < 0)
    y = mesh->ystart-1;
}

void BoundaryRegionYDown::nextY()
{
  y--;
  if(x > xe)
    x = xs;
}

bool BoundaryRegionYDown::isDone()
{
  return (x > xe) || (y < 0); // Return true if gone out of the boundary
}


///////////////////////////////////////////////////////////////


BoundaryRegionYUp::BoundaryRegionYUp(int xmin, int xmax)
{
  location = BNDRY_XOUT;
  y = mesh->ystart-1; // First point inside the boundary
  xs = xmin;
  xe = xmax;
  
  if(xe < xs)
    SWAP(xs, xe);
 
  // Unit vector out of the domain
  bx = 0;
  by = 1;
}

void BoundaryRegionYUp::first()
{
  x = xs;
  y = mesh->yend+1;
}

void BoundaryRegionYUp::next()
{
  // Loop over all points, from inside out
  y++;
  if(y >= mesh->ngy) {
    y = mesh->yend+1;
    x++;
  }
}

void BoundaryRegionYUp::nextX()
{
  x++;
  if(y >= mesh->ngy)
    y = mesh->yend+1;
}

void BoundaryRegionYUp::nextY()
{
  y++;
  if(x > xe)
    x = xs;
}

bool BoundaryRegionYUp::isDone()
{
  return (x > xe) || (y >= mesh->ngy); // Return true if gone out of the boundary
}
