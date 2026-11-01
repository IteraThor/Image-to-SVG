// WeightsFinder: explains a pixel as a convex mixture of the colours of the distinct segments in
// its 8- (or 24-) neighbourhood .
#pragma once

#include <array>
#include <functional>

#include "isvg/image.h"
#include "isvg/pixel_map.h"

namespace isvg {

struct Color4f {
    float c[4] = {0, 0, 0, 0};
};

class WeightsFinder {
public:
    // colorOf(id) returns the model colour of segment `id` (SubPixelSegmentation::getColorFRef).
    using ColorSource = std::function<const float*(int)>;

    void attach(PixelMap& map, const Image& img, ColorSource colorOf) {
        map_ = &map;
        img_ = &img;
        colorOf_ = std::move(colorOf);
    }
    // Palette-index mode (setUsePaletteIndexSource: segment ids read from the map
    // are replaced by idMap(id) — the region's palette index — before use.
    void setIdMap(std::function<int(int)> idMap) { idMap_ = std::move(idMap); }
    void clear();

    // Model colour of pixel (x,y). Ids equal to `replaceFrom` are treated as `replaceTo` when
    // replaceTo >= 0. *ownWeight receives the weight of the pixel's own segment.
    Color4f computePixelF(int x, int y, int replaceTo, int replaceFrom, double* ownWeight, bool use24);

    int numNeighbors() const { return k_; }
    const int* neighborIds() const { return ids_.data(); }
    double weight(int i) const { return w_[i]; }
    int heaviestSegmentId() const;

private:
    int idWithReplace(int x, int y, int replaceTo, int replaceFrom);  // paletteOrIdOfWithReplace
    int updateNeighbors(int x, int y, int replaceTo, int replaceFrom, bool use24);
    void setupQP();
    void solveQP();

    PixelMap* map_ = nullptr;
    const Image* img_ = nullptr;
    ColorSource colorOf_;
    std::function<int(int)> idMap_;

    int k_ = 0;  // number of distinct neighbour segments
    std::array<int, 9> ids_{};  // own segment first
    double w_[9] = {};  // mixture weights (w_[0] own, w_[1] other for K=2)
    double A_[9 * 3] = {};  // colours, column i = colour of ids_[i] (lda 3)
    double H_[9 * 9] = {};  // AᵀA + 0.01 I, then its Cholesky factor (ld 9)
    double f_[9] = {};  // Aᵀp
    double E_[9 * 9] = {};  // equality constraints, column t = constraint t
    double Y_[9 * 9] = {};  // H⁻¹E
    double S_[9 * 9] = {};  // EᵀH⁻¹E and its Cholesky factor
    double lambda_[9] = {};
    double r_[9] = {};
    bool pinned_[9] = {};
    double p_[4] = {};  // pixel in, model colour out
    //: only ever set to -1 (clear), so heaviestSegmentId always returns the pixel's own
    // segment — verified: nothing else writes this field.
    int numCandidates_ = -1;
};

}  // namespace isvg
