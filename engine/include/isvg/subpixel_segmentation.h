// SubPixelSegmentation: segments whose boundary pixels are modelled as anti-aliased mixtures of the
// neighbouring segment colours .
#pragma once

#include <set>
#include <vector>

#include "isvg/pixel_segmentation.h"
#include "isvg/weights_finder.h"

namespace isvg {

// SubPixelSegment.
struct SubPixelSegment {
    int n = 0;
    float sum[4] = {};
    int paletteIndex = -1;
    std::set<int> band;  // boundary-band pixel indices (inner and outer)
    int nInterior = 0;
    float sumSq[4] = {};
    float interiorSum[4] = {};
    float interiorSumSq[4] = {};
    bool flag = false;
    float color[4] = {};  // model colour (interior median or mean)
    std::vector<uint8_t> interiorBytes[4];  // per-channel interior samples for the median

    void clearStats();
    void addPixelF(const float c[4]);
    void interiorPushBack(const uint8_t* px);
    void interiorMedianPixelF(float out[4]);
};

struct SubPixelSegmentationParams {
    float self_weight_eps = 0.1f;  // SubPixelSegmentation::self_weight_eps
    float self_weight_penalty_low = 0.01f;  // SubPixelSegmentation::self_weight_penalty_low
    float self_weight_penalty_high = 1.0f;  // SubPixelSegmentation::self_weight_penalty_high
};

class SubPixelSegmentation : public PixelSegmentation {
public:
    SubPixelSegmentation(PixelMap& map, const Image& image, const SubPixelSegmentationParams& p);

    void initialize(int stage, bool resetPalette);
    void rebuildSubPixelSegments(int stage, bool use24, bool keepFlags);
    void updateAllPixelCosts();
    float computePixelCost(int x, int y, int replaceTo, int replaceFrom);
    // Interior colour error of segment `seg` if it took the colour of `colorSeg` (-1: its own).
    double computeInteriorCost(int seg, int colorSeg);
    // Sum of band pixel costs of `seg`: cached (useCached) or recomputed with `seg` read as `other`,
    // optionally storing the recomputed costs.
    double computeBoundaryCost(bool useCached, bool store, int seg, int other);
    double computeSegmentCost(int seg, int other) {
        return computeInteriorCost(seg, other) + computeBoundaryCost(other < 0, false, seg, other);
    }
    void mergeSegments(int into, int from);
    bool isInterior(int id, int pixel);  // PixelMap::isInterior with this->use24
    // Sum of band pixel costs of `seg` read as `other`, stopping once it exceeds `limit`.
    double computeBoundaryCostUpTo(double limit, int seg, int other);
    // Map pixels to palette colours and renumber `dst` into connected components.
    void flattenTo(PixelSegmentation& dst);
    void finalize();
    void resetPaletteIndices();

    SubPixelSegment& sub(int id) { return subs_.at(id); }
    int numSubSlots() const { return static_cast<int>(subs_.size()); }
    const float* colorOf(int id);  // getColorFRef
    std::vector<float>& pixelCosts() { return cost_; }
    WeightsFinder& weights() { return wf_; }
    bool use24() const { return use24_; }

protected:
    bool hasSegmentRecords() const override { return !subs_.empty(); }
    int& paletteIndexOf(int id) override { return subs_.at(id).paletteIndex; }
    void ownColorF(int id, float out[4]) override {
        for (int k = 0; k < 4; ++k) out[k] = subs_.at(id).color[k];
    }

private:
    SubPixelSegmentationParams p_;
    std::vector<SubPixelSegment> subs_;  // / count
    std::vector<float> cost_;  // per-pixel cost
    WeightsFinder wf_;
    bool use24_ = false;
};

}  // namespace isvg
