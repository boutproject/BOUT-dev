
#include "boundary_region.h"
#include "utils.h"
#include "globals.h"

BoundaryRegionXIn::BoundaryRegionXIn(const string &name, int ymin, int ymax)
  : ys(ymin), ye(ymax), BoundaryRegion(name, -1, 0)
{
  location = BNDRY_XIN;
  x = mesh->xstart-1; // First point inside the boundary
  if(ye < ys)
    SWAP(ys, ye);
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


BoundaryRegionXOut::BoundaryRegionXOut(const string &name, int ymin, int ymax)
  : ys(ymin), ye(ymax), BoundaryRegion(name, 1, 0)
{
  location = BNDRY_XOUT;
  x = mesh->xend+1; // First point inside the boundary
  if(ye < ys)
    SWAP(ys, ye);
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


BoundaryRegionYDown::BoundaryRegionYDown(const string &name, int xmin, int xmax)
  : xs(xmin), xe(xmax), BoundaryRegion(name, 0, -1)
{
  location = BNDRY_XOUT;
  y = mesh->ystart-1; // First point inside the boundary
  if(xe < xs)
    SWAP(xs, xe);
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


BoundaryRegionYUp::BoundaryRegionYUp(const string &name, int xmin, int xmax)
  : xs(xmin), xe(xmax), BoundaryRegion(name, 0, 1)
{
  location = BNDRY_XOUT;
  y = mesh->ystart-1; // First point inside the boundary
  if(xe < xs)
    SWAP(xs, xe);
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
