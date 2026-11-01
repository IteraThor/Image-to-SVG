// Image: 8-bit, 4-channel pixels in the byte order the engine receives them (B,G,R,A — Qt ARGB32).
#pragma once

#include <cstdint>
#include <vector>

namespace isvg {

struct Image {
    int w = 0, h = 0;
    std::vector<uint8_t> px;  // w*h*4

    const uint8_t* at(int index) const { return &px[static_cast<size_t>(index) * 4]; }
    int numPixels() const { return w * h; }
};

// Loading: the engine stores pixels premultiplied by alpha,
// truncating c * (a / 255.0) for the three colour channels.
inline void premultiplyOnLoad(Image& img) {
    for (size_t i = 0; i + 3 < img.px.size(); i += 4) {
        const double a = static_cast<double>(img.px[i + 3]) / 255.0;
        for (int k = 0; k < 3; ++k) img.px[i + k] = static_cast<uint8_t>(static_cast<int>(static_cast<double>(img.px[i + k]) * a));
    }
}

constexpr float kInv255f = 0.003921569f;  // 0x3b808081, the engine's 1/255 in float

}  // namespace isvg
