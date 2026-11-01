// PixelSegmenter: quadtree pre-segmentation + annealed merge/swap over-segmentation
// .
#pragma once

#include "isvg/pixel_segmentation.h"

namespace isvg {

struct PixelSegmenterParams {  // parameter-file keys in the comments
    int max_iterations = 50;  // Segmenter::max_iterations
    float lambda_pre = 1e-4f;  // Segmenter::lambda_pre
    float lambda_initial = 1e-4f;  // Segmenter::lambda_initial
    float lambda_final = 0.1f;  // Segmenter::lambda_final
    int min_num_pixels = 0;  // Segmenter::min_num_pixels
    int num_iter_btwn_renum = 5;  // Segmenter::num_iter_btwn_renum
    float rise_fraction = 0.4f;  // Segmenter::rise_fraction
    float size_penalty_weight = 0;  // Segmenter::size_penalty_weight
    float size_penalty_power = 1;  // Segmenter::size_penalty_power
};

enum class MergeConnectedness { Four = 0, Eight = 1, MinPixels = 2 };

class PixelSegmenter {
public:
    PixelSegmenter(PixelSegmentation& seg, const PixelSegmenterParams& p) : seg_(seg), p_(p) {}

    void initialize();
    void overSegmentation(bool diagonal);
    void enforceMinNumPixels();
    // Aliased branch: segments with the same truncated mean byte colour become one palette entry and
    // the pixel map is rewritten to palette indices.
    void clusterColorsAndGeneratePalette();

    // Exposed for stage-level testing.
    void preSegmentation();
    void mergeAndSwapLoop(MergeConnectedness mode);
    void stepLambda();
    float lambda() const { return lambda_; }
    float rate() const { return rate_; }

private:
    float computeMergeMetric(int a, int b);
    int considerMergeAndSwap(int p1, int p2, int a, int b);
    void swapPixel(int p, int from, int to);
    void mergeSegments(int into, int from);

    PixelSegmentation& seg_;
    PixelSegmenterParams p_;
    float lambda_ = 0;
    float rate_ = 1;
    int w_ = 0, h_ = 0;
};

}  // namespace isvg
