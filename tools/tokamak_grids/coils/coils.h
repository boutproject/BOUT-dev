
#ifndef __COILS_H__
#define __COILS_H__

class Point {
 public:
  Point() : x(0.), y(0.), z(0.) {}
  Point(double px, double py, double pz) : x(px), y(py), z(pz) {}
  
  double norm() const;
  double distanceTo(const Point &p);
  
  const Point operator+() const;
  const Point operator+(const Point &p) const;
  Point & operator+=(const Point &p);
  
  const Point operator-() const;
  const Point operator-(const Point &p) const;
  Point & operator-=(const Point &p);
  
  const Point operator/(double val) const;
  Point & operator/=(const Point &p);
  Point & operator/=(double val);
  
  const Point operator*(double val) const;
  Point & operator*=(const Point &p);
  Point & operator*=(double val);
  
  double x, y, z;
};

const Point operator*(const double lhs, const Point &rhs);

class Vector {
 public:
  Vector() {}
  Vector(Point o, Point d) : origin(o), direction(d) {}
  
  Vector & operator+=(const Vector &v);

  Point origin;    // Starting point for the vector 
  Point direction; // Direction the vector's pointing in (from origin at point)
};

#endif // __COILS_H__
