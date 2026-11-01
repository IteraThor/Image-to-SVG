// PaletteFinder: k-means palettes for 2..12 colours
// on the (downsampled) image and the 62 features classifyPaletteMode decides the colour count from.
#pragma once

#include <array>
#include <string>
#include <vector>

#include "isvg/image.h"
#include "isvg/parameter_classifier.h"
#include "isvg/pixel_map.h"

namespace isvg {

struct PaletteFinderParams {
    float costPerColor = 0.15f;  // PaletteFinder::cost_per_color (EP_CLASSIFY)
    bool use16HistBins = true;  // PaletteFinder::use_16_hist_bins
    int maskLevel = 2;  // PaletteFinder
};

// Downsampling: nearest-neighbour downsampling to ~150000
// pixels when the image has more than 200000.
Image downsamplePixels(const Image& src);

class PaletteFinder {
public:
    // img: the (downsampled) pixels; map: the engine's pixel map (reused from the classify
    // segmentation, redirect table and all); edges/edgeW: the full-resolution EdgeImage, which the
    // original keeps indexing with downsampled coordinates.
    PaletteFinder(const Image& img, PixelMap& map, const std::vector<uint8_t>& edges, int edgeW,
                  bool antiAliased, const PaletteFinderParams& p, const EdgeParams& ep)
        : img_(img), map_(map), edges_(edges), edgeW_(edgeW), aa_(antiAliased), p_(p), ep_(ep) {}

    // findPalettes for k = 2..12, then findOptimalNumColors.
    // Returns classifyPaletteMode's answer; features holds its 62 inputs.
    std::string findPalettes(const std::string& ext, const std::string& typeString);

    // The k-means part alone (findPalettes stage): fresh histogram, then
    // findPalette(k) for k = kMin..kMax.
    void findPalettesRange(int kMin, int kMax);

    const double* features() const { return features_; }
    const FoundPalette& palette(int k) const { return palettes_[k]; }

private:
    void doHistogram();
    void findPalette(int k);
    void doKMeans(int k);
    void initializeK(int k);
    void assignK(bool final);
    void computeMeanColors(bool weighted);
    void computeCost();
    float computePixelCost(int x, int y, const std::vector<uint8_t>& mask, class WeightsFinder& wf);
    void computeFeatures();

    double edgeWeightLinear(int i) const { return edgeWeight(edges_[i], ep_.scale, ep_.floor); }
    double edgeWeightXY(int x, int y) const { return edgeWeight(edges_[y * edgeW_ + x], ep_.scale, ep_.floor); }
    double edgeIntensity(int i) const { return static_cast<double>(edges_[i]) / 255.0; }

    const Image& img_;
    PixelMap& map_;
    const std::vector<uint8_t>& edges_;
    int edgeW_;
    bool aa_;
    PaletteFinderParams p_;
    EdgeParams ep_;

    std::vector<std::array<uint8_t, 4>> peaks_;
    std::vector<double> peakWeights_;
    FoundPalette cur_;  // current palette
    std::vector<std::array<float, 4>> means_;
    std::vector<float> meanWeights_;
    FoundPalette palettes_[13];
    double features_[62] = {};
};

}  // namespace isvg
