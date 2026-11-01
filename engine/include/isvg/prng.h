// Global PRNG: the engine's only random source, a process-global Lehmer generator
// (s = s*2003 mod 9973, s0 = 1). The application never reseeds it, so only the first vectorization in a
// session starts from s0; the port resets it per run (resetPrng) to match a fresh launch.
#pragma once

#include <cstdint>

namespace isvg {

inline int32_t& prngState() {
    static int32_t s = 1;
    return s;
}
inline void resetPrng() { prngState() = 1; }
inline double prngNext() {
    int32_t& s = prngState();
    s = (s * 2003) % 9973;
    return static_cast<double>(s) / 9973.0;
}

}  // namespace isvg
