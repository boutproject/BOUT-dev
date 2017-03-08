
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
  width = mesh->local_nx - mesh->xend - 1;
  x = mesh->local_nx - width; // First point inside the boundary
  if(ye < ys)
    swap(ys, ye);
}

void BoundaryRegionXOut::first()
{
  x = mesh->local_nx - width;
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
  if(x >= mesh->local_nx)
    x = mesh->local_nx - width;
}

bool BoundaryRegionXOut::isDone()
{
  return (x >= mesh->local_nx) || (y > ye); // Return true if gone out of the boundary
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
  width = mesh->local_ny - mesh->yend - 1;
  y = mesh->local_ny - width; // First point inside the boundary
  if(xe < xs)
    swap(xs, xe);
}

void BoundaryRegionYUp::first()
{
  x = xs;
  y = mesh->local_ny - width;
}

void BoundaryRegionYUp::next()
{
  // Loop over all points, from inside out
  y++;
  if(y >= mesh->local_ny) {
    y = mesh->local_ny - width;
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
  if(y >= mesh->local_ny)
    y = mesh->local_ny - width;
}

void BoundaryRegionYUp::nextY()
{
  y++;
  if(x > xe)
    x = xs;
}

bool BoundaryRegionYUp::isDone()
{
  return (x > xe) || (y >= mesh->local_ny); // Return true if gone out of the boundary
}
