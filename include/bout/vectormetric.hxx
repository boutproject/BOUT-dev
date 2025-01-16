#pragma once

#include "bout/build_defines.hxx"

#if BOUT_USE_METRIC_3D
#include "bout/vector3d.hxx"
using VectorMetric = Vector3D;
#else
#include "bout/vector2d.hxx"
using VectorMetric = Vector2D;
#endif
