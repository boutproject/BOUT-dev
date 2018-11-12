
#include <globals.hxx>
#include <boundary_region.hxx>
#include <utils.hxx>

BoundaryRegionXIn::BoundaryRegionXIn(std::string name, int ymin, int ymax, Mesh* passmesh, int wid)
  : BoundaryRegion(name, -1, 0, BNDRY_XIN,
      wid < 0 ? (passmesh == nullptr ? mesh : passmesh)->xstart : wid,
      passmesh),
    ys(ymin), ye(ymax)
{
  x = width-1; // First point inside the boundary
  if(ye < ys)
    swap(ys, ye);
}

BoundaryRegion* BoundaryRegionXIn::copy(int wid) {
  return new BoundaryRegionXIn(label, ys, ye, localmesh, wid);
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


BoundaryRegionXOut::BoundaryRegionXOut(std::string name, int ymin, int ymax, Mesh* passmesh, int wid)
  : BoundaryRegion(name, 1, 0, BNDRY_XOUT,
      wid < 0 ? (passmesh == nullptr ? mesh : passmesh)->LocalNx - (passmesh == nullptr ? mesh : passmesh)->xend - 1 : wid,
      passmesh),
    ys(ymin), ye(ymax)
{
  x = localmesh->LocalNx - width; // First point inside the boundary
  if(ye < ys)
    swap(ys, ye);
}

BoundaryRegion* BoundaryRegionXOut::copy(int wid) {
  return new BoundaryRegionXOut(label, ys, ye, localmesh, wid);
}

void BoundaryRegionXOut::first()
{
  x = localmesh->LocalNx - width;
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
  if(x >= localmesh->LocalNx)
    x = localmesh->LocalNx - width;
}

bool BoundaryRegionXOut::isDone()
{
  return (x >= localmesh->LocalNx) || (y > ye); // Return true if gone out of the boundary
}

///////////////////////////////////////////////////////////////


BoundaryRegionYDown::BoundaryRegionYDown(std::string name, int xmin, int xmax, Mesh* passmesh, int wid)
  : BoundaryRegion(name, 0, -1, BNDRY_YDOWN,
      wid < 0 ? (passmesh == nullptr ? mesh : passmesh)->ystart : wid,
      passmesh),
    xs(xmin), xe(xmax)
{
  y = width-1; // First point inside the boundary
  if(xe < xs)
    swap(xs, xe);
}

BoundaryRegion* BoundaryRegionYDown::copy(int wid) {
  return new BoundaryRegionYDown(label, xs, xe, localmesh, wid);
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


BoundaryRegionYUp::BoundaryRegionYUp(std::string name, int xmin, int xmax, Mesh* passmesh, int wid)
  : BoundaryRegion(name, 0, 1, BNDRY_YUP,
      wid < 0 ? (passmesh == nullptr ? mesh : passmesh)->LocalNy - (passmesh == nullptr ? mesh : passmesh)->yend - 1 : wid,
      passmesh),
    xs(xmin), xe(xmax)
{
  y = localmesh->LocalNy - width; // First point inside the boundary
  if(xe < xs)
    swap(xs, xe);
}

BoundaryRegion* BoundaryRegionYUp::copy(int wid) {
  return new BoundaryRegionYUp(label, xs, xe, localmesh, wid);
}

void BoundaryRegionYUp::first()
{
  x = xs;
  y = localmesh->LocalNy - width;
}

void BoundaryRegionYUp::next()
{
  // Loop over all points, from inside out
  y++;
  if(y >= localmesh->LocalNy) {
    y = localmesh->LocalNy - width;
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
  if(y >= localmesh->LocalNy)
    y = localmesh->LocalNy - width;
}

void BoundaryRegionYUp::nextY()
{
  y++;
  if(x > xe)
    x = xs;
}

bool BoundaryRegionYUp::isDone()
{
  return (x > xe) || (y >= localmesh->LocalNy); // Return true if gone out of the boundary
}
