
#include <globals.hxx>
#include <boundary_region.hxx>
#include <utils.hxx>

BoundaryRegionXIn::BoundaryRegionXIn(const string &name, int ymin, int ymax)
  : ys(ymin), ye(ymax), BoundaryRegion(name, -1, 0)
{
  location = BNDRY_XIN;
  width = mesh->xstart;
  x = width-1; // First point inside the boundary
  if(ye < ys)
    swap(ys, ye);
}

void BoundaryRegionXIn::first()
{
  x = width-1;
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
    x = width-1;
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
  width = mesh->ngx - mesh->xend - 1;
  x = mesh->ngx - width; // First point inside the boundary
  if(ye < ys)
    swap(ys, ye);
}

void BoundaryRegionXOut::first()
{
  x = mesh->ngx - width;
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
    x = mesh->ngx - width;
}

bool BoundaryRegionXOut::isDone()
{
  return (x >= mesh->ngx) || (y > ye); // Return true if gone out of the boundary
}

///////////////////////////////////////////////////////////////


BoundaryRegionYDown::BoundaryRegionYDown(const string &name, int xmin, int xmax)
  : xs(xmin), xe(xmax), BoundaryRegion(name, 0, -1)
{
  location = BNDRY_YDOWN;
  width = mesh->ystart;
  y = width-1; // First point inside the boundary
  if(xe < xs)
    swap(xs, xe);
}

void BoundaryRegionYDown::first()
{
  x = xs;
  y = width-1;
}

void BoundaryRegionYDown::next()
{
  // Loop over all points, from inside out
  y--;
  if(y < 0) {
    y = width-1;
    x++;
  }
}

void BoundaryRegionYDown::nextX()
{
  x++;
  if(y < 0)
    y = width-1;
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
  location = BNDRY_YUP;
  width = mesh->ngy - mesh->yend - 1;
  y = mesh->ngy - width; // First point inside the boundary
  if(xe < xs)
    swap(xs, xe);
}

void BoundaryRegionYUp::first()
{
  x = xs;
  y = mesh->ngy - width;
}

void BoundaryRegionYUp::next()
{
  // Loop over all points, from inside out
  y++;
  if(y >= mesh->ngy) {
    y = mesh->ngy - width;
    x++;
  }
}

void BoundaryRegionYUp::nextX()
{
  x++;
  if(y >= mesh->ngy)
    y = mesh->ngy - width;
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
