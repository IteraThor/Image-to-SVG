// ParameterClassifier: picks the image type, quality and palette size that
// classifyImage reports and that select the engine parameter file.
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "isvg/image.h"
#include "isvg/pixel_map.h"
#include "isvg/pixel_segmentation.h"

namespace isvg {

// Decision trees (classifier_trees.cpp, generated).
// features: Vec<float,20> from ParameterClassifier::classifyParameters; ext: image file extension
// (".png", ...). Returns "<type> <quality>", e.g. "logo_aa 2", "logo 2", "photo 1".
std::string classifyImageType(const float* features, const std::string& ext);
// features: Vec<double,62> from PaletteFinder::findPalettes; typeString: classifyImageType's result.
// Returns the palette mode as a number string ("2".."9", "13").
std::string classifyPaletteMode(const double* features, const std::string& typeString);

// ColorHistogram with 16 bins per channel (setNumBins(16): shift 4), indexed by the colour bytes
// in memory order.
class ColorHistogram {
public:
    static constexpr int kBins = 16, kShift = 4;
    void setNumBins16();
    void addColor(const uint8_t c[4], double w);
    double totalWeight() const;
    // Zero the colour's bin, or (neighborhood) the 3^4 bins around it.
    void zeroOutColor(const uint8_t c[4], bool neighborhood);
    // Local maxima over the 80-neighbourhood with weight >= 2, as bin-centre colours, sorted by
    // descending weight with connectivity 1).
    void findPeaks(std::vector<std::array<uint8_t, 4>>& colors, std::vector<double>& weights) const;
    std::vector<double>& bins() { return bins_; }

private:
    static int index(int a, int b, int c, int d) { return a + kBins * (b + kBins * (c + kBins * d)); }
    std::vector<double> bins_;
};

// EdgeImage::update: byte edge strength from the 4-neighbour colour differences.
std::vector<uint8_t> edgeImage(const Image& img);
// Edge weight of one pixel.
double edgeWeight(uint8_t e, double scale, double floor);

struct EdgeParams {
    double scale = 42.0;  // EdgeImage::edge_weight_scale
    double floor = 0.0032;  // EdgeImage::edge_weight_floor
};

// The 20 features ParameterClassifier::classifyParameters feeds to classifyImageType, computed from
// the preprocessed pixels and the classify-mode segmentation.
void computeImageTypeFeatures(const Image& img, PixelMap& map, PixelSegmentation& seg,
                              const std::vector<uint8_t>& edges, const EdgeParams& ep, float features[20]);

struct FoundPalette {
    std::vector<std::array<float, 4>> colors;
    std::vector<double> weights;  // per colour: pixel count, or summed edge intensity (anti-aliased)
    float cost = 0, logCost = 0, colorCost = 0;  // Palette,,
};

// classifyImage's result: the four ints it stores at.
struct Classification {
    int imageType = 0;  // 0 photo, 1 logo_aa, 2 logo
    int quality = 0;  // 0 high, 1 medium, 2 low
    int paletteMode = 0;  // 0 two colours, 1 up to maxColors, 2 full colour
    int maxColors = 0;
    std::string typeString, paletteString;
    float typeFeatures[20];
    double paletteFeatures[62];
    FoundPalette palettes[13];  // PaletteFinder palettes for k = 2..12
};
// the engine::classifyImage on preprocessed pixels; ext is the image file extension (".png").
// Advances the global PRNG deterministically.
// `map` is the engine's pixel map; it carries state (redirect table) into findPalettesStage.
Classification classifyImage(const Image& img, const std::string& ext, PixelMap& map);

// findPalettes(settings, 2, 12) as the pipeline calls it : k-means palettes for k = 2..12 on the downsampled image with the vectorization parameters.
// Its only effect on the output is the global PRNG it advances.
void findPalettesStage(const Image& img, PixelMap& map, bool antiAliased, float costPerColor);

}  // namespace isvg
