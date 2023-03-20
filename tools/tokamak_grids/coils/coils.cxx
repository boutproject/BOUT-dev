/************************************************************
 * Coil calculation to replace coils.pro
 ************************************************************/

#include "coils.h"
#include <cmath>
#include <iostream>
#include <netcdfcpp.h>
#include <vector>

using namespace std;

const double PI = 3.1415926;

/////////////////////////////////////////////////////////////
// Point implementation

double Point::norm() const { return sqrt(x * x + y * y + z * z); }

double Point::distanceTo(const Point& p) { return (*this - p).norm(); }

const Point Point::operator+() const { return *this; }

const Point Point::operator+(const Point& p) const {
  return Point(x + p.x, y + p.y, z + p.z);
}

Point& Point::operator+=(const Point& p) {
  x += p.x;
  y += p.y;
  z += p.z;
  return *this;
}

const Point Point::operator-() const { return Point(-x, -y, -z); }

const Point Point::operator-(const Point& p) const {
  return Point(x + p.x, y + p.y, z + p.z);
}

Point& Point::operator-=(const Point& p) {
  x -= p.x;
  y -= p.y;
  z -= p.z;
  return *this;
}

const Point Point::operator/(double val) const {
  return Point(x / val, y / val, z / val);
}

Point& Point::operator/=(double val) {
  x /= val;
  y /= val;
  z /= val;
  return *this;
}

const Point Point::operator*(double val) const {
  return Point(x * val, y * val, z * val);
}

Point& Point::operator*=(double val) {
  x *= val;
  y *= val;
  z *= val;
  return *this;
}

const Point operator*(const double lhs, const Point& rhs) { return rhs * lhs; }

/////////////////////////////////////////////////////////////
// Vector implementation

Vector& Vector::operator+=(const Vector& v) {
  direction += v.direction;
  return *this;
}

/////////////////////////////////////////////////////////////

// Calculate vector potential from a line at given position
const Point AfromLine(Point start, Point end, double current, Point pos) {
  double len = start.distanceTo(end);

  Point ivec = current * (end - start) / len; // Current vector

  // Integrate ivec * 1/d over wire
  double integral = 0., last;
  int n = 1;
  do {
    last = integral;
    n *= 2;

    // Use Simpson's rule to integrate
    double h = len / ((double)n);

    integral = 0.;
    for (int i = 0; i <= n; i++) {
      double frac = ((double)i) / ((double)(n - 1));
      double d = pos.distanceTo(frac * start + (1. - frac) * end);

      if ((i == 0) || (i == (n - 1))) {
        integral += 1. / d;
      } else if (i % 2 == 1) {
        integral += 4. / d;
      } else
        integral += 2. / d;
      integral *= h / 3.;
    }
  } while (abs((integral - last) / (integral + last)) > 1.e-3);

  return ivec * integral;
}

const Point AfromCoil(vector<Point> corners, double current, Point pos) {
  Point A;

  for (int i = 0; i < corners.size(); i++)
    A += AfromLine(corners[i], corners[(i + 1) % corners.size()], current, pos);

  return A;
}

/////////////////////////////////////////////////////////////
// Main

double** matrix(int nx, int ny) {
  double** m;
  m = new double*[nx];
  m[0] = new double[nx * ny];
  for (int i = 1; i < nx; i++)
    m[i] = m[i - 1] + ny;
  return m;
}

int readInteger(NcFile* dataFile, const char* name) {
  NcVar* var = dataFile->get_var(name);
  if (!var) {
    cout << "ERROR: Couldn't read " << string(name) << endl;
    return NULL;
  }
  if (!var->is_valid()) {
    cout << "ERROR: Couldn't read " << string(name) << endl;
    return NULL;
  }

  int result;
  var->get(&result, 1);
  return result;
}

double** readMatrix(NcFile* dataFile, const char* name, int& nx, int& ny) {
  NcVar* var = dataFile->get_var(name);
  if (!var) {
    cout << "ERROR: Couldn't read " << string(name) << endl;
    return NULL;
  }
  if (!var->is_valid()) {
    cout << "ERROR: Couldn't read " << string(name) << endl;
    return NULL;
  }
  if (var->num_dims() != 2) {
    cout << "ERROR: " << string(name) << " must be 2D" << endl;
    return NULL;
  }

  NcDim *xDim, *yDim;
  xDim = var->get_dim(0);
  yDim = var->get_dim(1);

  nx = xDim->size();
  ny = yDim->size();

  double** m = matrix(nx, ny);
  var->get(m[0], nx, ny);

  return m;
}

int main(int argc, char** argv) {
  char* name;

  // Check command-line arguments

  if (argc == 1) {
    cout << "Useage: " << string(argv[0]) << " <grid file>" << endl;
    return 1;
  }
  name = argv[1];

  NcFile* dataFile;
  NcDim *xDim, *yDim;

  NcError err(NcError::verbose_nonfatal);
  dataFile = new NcFile(name, NcFile::ReadOnly);
  if (!dataFile->is_valid()) {
    delete dataFile;
    cout << "ERROR: Couldn't open grid file" << endl;
    return 1;
  }
  if (!(xDim = dataFile->get_dim("x"))) {
    cout << "ERROR: Grid file doesn't have an x dimension" << endl;
    return 1;
  }
  if (!(yDim = dataFile->get_dim("y"))) {
    cout << "ERROR: Grid file doesn't have a y dimension" << endl;
    return 1;
  }

  int nx, ny;

  double** Rxy = readMatrix(dataFile, "Rxy", nx, ny);
  if (!Rxy)
    return 1;
  double** Zxy = readMatrix(dataFile, "Zxy", nx, ny);
  if (!Zxy)
    return 1;

  int nz = 16; // Number of points to use in Z

  // Loop over the grid points
  for (int x = 0; x < nx; x++)
    for (int y = 0; y < ny; y++)
      for (int z = 0; z < nz; z++) {
        double phi = 2. * PI * ((double)z) / ((double)nz);

        // Convert to cartesian coordinates
        Point p(Rxy[x][y] * cos(phi), Rxy[x][y] * sin(phi), Zxy[x][y]);

        // Calculate A
      }

  return 0;
}
