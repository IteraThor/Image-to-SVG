// Shadow detection support: Gaussian-blurred image, signal images and boundary mask
// (blurred image, SignalImage, BoundaryMask) and the trained shadow classifier.
#pragma once

#include <array>
#include <vector>

#include "isvg/image.h"
#include "isvg/pixel_map.h"

namespace isvg {

// Blur: 5×5 Gaussian (σ = 0.75), renormalised at borders.
Image blurPixels(const Image& src);

// SignalImage::initialize. s1 = squared colour gradient to the 4 neighbours of the
// original image; s2 = squared Laplacian of the blurred image. Borders are mirrored.
struct SignalImages {
    int w = 0, h = 0;
    std::vector<float> s1, s2;
};
SignalImages computeSignals(const Image& src);

// BoundaryMask::update8 over raw map ids: `level` on every pixel of a 2×2 block with
// differing ids, then rings of level-1 ... 1 around it.
std::vector<uint8_t> boundaryMask8(const PixelMap& map, int level);

using FeatureVector = std::array<double, 42>;

// isProbablyShadowingSegment (decision tree with trained thresholds).
bool isProbablyShadowingSegment(const FeatureVector& f);
// isProbablyNotFalsePositive: accept a shadow dissolution (uses features 33-38).
bool isProbablyNotFalsePositive(const FeatureVector& f, int configT);

}  // namespace isvg
