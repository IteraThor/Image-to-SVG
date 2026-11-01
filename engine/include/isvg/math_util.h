// libm entry points as the engine calls them (logf = (float)log((double)x)).
// Results depend on the C library's exp/log/pow/atan2; this port relies on the Emscripten libm.
#pragma once

#include <cmath>

namespace isvg::math {

inline double log(double x) { return std::log(x); }
inline double exp(double x) { return std::exp(x); }
inline double pow(double x, double y) { return std::pow(x, y); }
inline float logf(float x) { return static_cast<float>(log(static_cast<double>(x))); }

}  // namespace isvg::math
