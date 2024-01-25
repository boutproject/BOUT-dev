/// \file
/// Header to simplify using `Vector2D`/`Vector3D` when switching
/// between 2D and 3D metrics in user code

#pragma once

#include "bout/build_config.hxx"

#if BOUT_USE_METRIC_3D
#include "bout/vector3d.hxx"
using VectorMetric = Vector3D;
#else
#include "bout/vector2d.hxx"
using VectorMetric = Vector2D;
#endif
