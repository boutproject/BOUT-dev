
#include <globals.hxx>
#include <boundary_region.hxx>
#include <utils.hxx>

BoundaryRegionXIn::BoundaryRegionXIn(const string &name, int ymin, int ymax)
  : BoundaryRegion(name, -1, 0), ys(ymin), ye(ymax)
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

void BoundaryRegionXIn::next1d()
{
  // Loop over the innermost points
  y++;
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
  : BoundaryRegion(name, 1, 0), ys(ymin), ye(ymax)
{
  location = BNDRY_XOUT;
  width = mesh->LocalNx - mesh->xend - 1;
  x = mesh->LocalNx - width; // First point inside the boundary
  if(ye < ys)
    swap(ys, ye);
}

void BoundaryRegionXOut::first()
{
  x = mesh->LocalNx - width;
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

void BoundaryRegionXOut::next1d()
{
  // Loop over the innermost points
  y++;
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
  if(x >= mesh->LocalNx)
    x = mesh->LocalNx - width;
}

bool BoundaryRegionXOut::isDone()
{
  return (x >= mesh->LocalNx) || (y > ye); // Return true if gone out of the boundary
}

///////////////////////////////////////////////////////////////


BoundaryRegionYDown::BoundaryRegionYDown(const string &name, int xmin, int xmax)
  : BoundaryRegion(name, 0, -1), xs(xmin), xe(xmax)
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

void BoundaryRegionYDown::next1d()
{
  // Loop over the innermost points
  x++;
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
  : BoundaryRegion(name, 0, 1), xs(xmin), xe(xmax)
{
  location = BNDRY_YUP;
  width = mesh->LocalNy - mesh->yend - 1;
  y = mesh->LocalNy - width; // First point inside the boundary
  if(xe < xs)
    swap(xs, xe);
}

void BoundaryRegionYUp::first()
{
  x = xs;
  y = mesh->LocalNy - width;
}

void BoundaryRegionYUp::next()
{
  // Loop over all points, from inside out
  y++;
  if(y >= mesh->LocalNy) {
    y = mesh->LocalNy - width;
    x++;
  }
}

void BoundaryRegionYUp::next1d()
{
  // Loop over the innermost points
  x++;
}

void BoundaryRegionYUp::nextX()
{
  x++;
  if(y >= mesh->LocalNy)
    y = mesh->LocalNy - width;
}

void BoundaryRegionYUp::nextY()
{
  y++;
  if(x > xe)
    x = xs;
}

bool BoundaryRegionYUp::isDone()
{
  return (x > xe) || (y >= mesh->LocalNy); // Return true if gone out of the boundary
}
