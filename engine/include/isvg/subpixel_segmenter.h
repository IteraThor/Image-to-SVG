// SubPixelSegmenter: greedy merging of over-segmented regions under the sub-pixel cost model
// .
#pragma once

#include <climits>
#include <map>

#include "isvg/shadows.h"
#include "isvg/subpixel_segmentation.h"

namespace isvg {

struct SubPixelSegmenterParams {  // parameter-file settings
    int max_iterations = 20;  // (not a registered parameter; constructor default)
    float beach_lambda_f = 1e-5f;  // SubPixelSegmenter::beach_lambda_f
    float color_cluster_margin = 0.1f;  // SubPixelSegmenter::color_cluster_margin
    int use_connected_24_iteration = 0;  // SubPixelSegmenter::use_connected_24_iteration
    float lambda_f = 0.1f;  // SubPixelSegmenter::lambda_f
};

class SubPixelSegmenter {
public:
    SubPixelSegmenter(SubPixelSegmentation& seg, const SubPixelSegmenterParams& p) : seg_(seg), p_(p) {}

    void execute();
    // Runs SubPixelSegmenter::execute up to (not including) splitAndMergeBeaches.
    void executeMergePhases();
    void splitAndMergeBeaches();
    // Small segments (< 3 interior pixels) adopt the colour of one of the 100 largest segments when
    // that lowers their cost by more than color_cluster_margin.
    void clusterColors(int stage);
    void finalize();
    // Per-segment 42-dimensional feature vectors and the shadow flags derived from
    // them.
    void updateAllFeatureVectors();
    void flagShadowSegments();
    // Dissolve shadow-flagged segments into their neighbours pixel by pixel, keep the change only
    // when isProbablyNotFalsePositive accepts it. `configT` is.
    void splitAndMergeFlaggedSegmentsFromBlurryImage(int configT);
    const std::vector<FeatureVector>& features() const { return features_; }
    int mergeLoop(int maxNeighborSize, int maxSize);
    double computeMergeMetric(int nb, int seg);
    bool useInteriorCount() const { return useInterior_; }

private:
    SubPixelSegmentation& seg_;
    SubPixelSegmenterParams p_;
    bool useInterior_ = false;  // (size measured by interior pixels; 24-neighbour pass)
    std::map<int, bool> evaluated_;  // pair cache keyed by nb | seg << 16
    std::vector<FeatureVector> features_;  // / count
};

}  // namespace isvg
