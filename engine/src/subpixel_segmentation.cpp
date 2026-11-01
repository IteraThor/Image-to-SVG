#include "isvg/subpixel_segmentation.h"

#include <algorithm>
#include <cassert>

namespace isvg {
namespace {
constexpr int kN8i[8] = {1, 0, -1, 0, -1, 1, -1, 1};
constexpr int kN8j[8] = {0, 1, 0, -1, 1, -1, -1, 1};
constexpr int kN24i[24] = {1, 0, -1, 0, -1, 1, -1, 1, -2, -2, -2, -2, -2, 2, 2, 2, 2, 2, -1, 0, 1, -1, 0, 1};
constexpr int kN24j[24] = {0, 1, 0, -1, 1, -1, -1, 1, -2, -1, 0, 1, 2, -2, -1, 0, 1, 2, -2, -2, -2, 2, 2, 2};
}  // namespace

void SubPixelSegment::clearStats() {
    n = 0;
    for (int k = 0; k < 4; ++k) {
        sum[k] = sumSq[k] = interiorSum[k] = interiorSumSq[k] = color[k] = 0.0f;
        interiorBytes[k].clear();
    }
    band.clear();
    nInterior = 0;
    flag = false;
}

void SubPixelSegment::addPixelF(const float c[4]) {
    ++n;
    for (int k = 0; k < 4; ++k) sum[k] = sum[k] + c[k];
    for (int k = 0; k < 4; ++k) sumSq[k] = c[k] * c[k] + sumSq[k];
}

void SubPixelSegment::interiorPushBack(const uint8_t* px) {
    for (int k = 0; k < 4; ++k) interiorBytes[k].push_back(px[k]);
}

// Per-channel nth_element at size/2 (the selected value does not depend on the algorithm).
void SubPixelSegment::interiorMedianPixelF(float out[4]) {
    const size_t mid = static_cast<size_t>(static_cast<int>(interiorBytes[0].size()) / 2);
    for (int k = 0; k < 4; ++k) {
        auto& v = interiorBytes[k];
        std::nth_element(v.begin(), v.begin() + mid, v.end());
        out[k] = static_cast<float>(v[mid]) * kInv255f;
    }
}

SubPixelSegmentation::SubPixelSegmentation(PixelMap& map, const Image& image, const SubPixelSegmentationParams& p)
    : PixelSegmentation(map, image), p_(p) {
    wf_.attach(map, image, [this](int id) { return colorOf(id); });
}

const float* SubPixelSegmentation::colorOf(int id) {
    SubPixelSegment& s = subs_.at(id);
    if (s.paletteIndex < 0) return s.color;
    return palette_->color(s.paletteIndex);
}

double SubPixelSegmentation::computeBoundaryCostUpTo(double limit, int seg, int other) {
    SubPixelSegment& s = subs_.at(seg);
    double sum = 0.0;
    const int w = image_.w;
    for (int q : s.band) {
        sum = sum + static_cast<double>(computePixelCost(q % w, q / w, other, seg));
        if (limit < sum) return sum;
    }
    return sum;
}

void SubPixelSegmentation::flattenTo(PixelSegmentation& dst) {
    if (!(initialized_ && hasSegmentRecords() && paletteIndexOf(0) >= 0)) generatePaletteFromSegments();
    const int n = image_.numPixels();
    for (int idx = 0; idx < n; ++idx) map_.raw(idx) = subs_.at(map_.raw(idx)).paletteIndex;
    int maxId = map_.maxId();
    if (maxId <= palette_->size()) maxId = palette_->size();
    dst.updateMaxId(maxId);
    for (int i = 0; i < palette_->size(); ++i) dst.segment(i).paletteIndex = i;
    dst.updateStatistics();
}

void SubPixelSegmentation::finalize() {
    initialized_ = false;
    segments_.clear();
    subs_.clear();
    cost_.clear();
}

void SubPixelSegmentation::resetPaletteIndices() {
    for (int i = 0; i <= map_.maxId(); ++i) subs_.at(i).paletteIndex = -1;
}

void SubPixelSegmentation::initialize(int stage, bool resetPalette) {
    wf_.clear();
    const int want = map_.maxId() + 1;
    if (want != static_cast<int>(subs_.size())) subs_.assign(std::max(0, want), SubPixelSegment{});
    cost_.assign(static_cast<size_t>(image_.numPixels()), 0.0f);
    if (resetPalette)
        for (int i = 0; i <= map_.maxId(); ++i) subs_[i].paletteIndex = -1;
    rebuildSubPixelSegments(stage, false, false);
    initialized_ = true;
}

void SubPixelSegmentation::rebuildSubPixelSegments(int /*stage*/, bool use24, bool keepFlags) {
    if (!keepFlags) renumber();
    use24_ = use24;
    const bool usingPalette = initialized_ && hasSegmentRecords() && paletteIndexOf(0) >= 0;
    assert(map_.maxId() < static_cast<int>(subs_.size()));
    for (int i = 0; i <= map_.maxId(); ++i) {
        SubPixelSegment& s = subs_[i];
        const bool f = s.flag;
        s.clearStats();
        if (keepFlags) s.flag = f;
    }
    const int w = map_.width(), h = map_.height();
    const int nn = use24_ ? 24 : 8;
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            const int idx = y * w + x;
            const int id = map_.raw(idx);
            SubPixelSegment& s = subs_[id];
            bool boundary = false;
            for (int k = 0; k < nn; ++k) {
                const int nx = x + (use24_ ? kN24i[k] : kN8i[k]);
                const int ny = y + (use24_ ? kN24j[k] : kN8j[k]);
                if (nx < 0 || nx >= image_.w || ny < 0 || ny >= image_.h) continue;
                const int q = nx + w * ny;
                if (map_.raw(q) == id) continue;
                s.band.insert(q);
                boundary = true;
            }
            const uint8_t* px = image_.at(idx);
            const float c[4] = {px[0] * kInv255f, px[1] * kInv255f, px[2] * kInv255f, px[3] * kInv255f};
            s.addPixelF(c);
            if (boundary) {
                s.band.insert(idx);
            } else {
                ++s.nInterior;
                for (int k = 0; k < 4; ++k) s.interiorSum[k] = s.interiorSum[k] + c[k];
                for (int k = 0; k < 4; ++k) s.interiorSumSq[k] = s.interiorSumSq[k] + c[k] * c[k];
                if (!usingPalette) s.interiorPushBack(px);
            }
        }
    }
    for (int i = 0; i <= map_.maxId(); ++i) {
        SubPixelSegment& s = subs_[i];
        if (usingPalette) continue;
        if (s.nInterior <= 2) {
            const float n = static_cast<float>(s.n);
            for (int k = 0; k < 4; ++k) s.color[k] = s.sum[k] / n;
        } else {
            s.interiorMedianPixelF(s.color);
        }
        for (int k = 0; k < 4; ++k) s.interiorBytes[k].clear();
    }
    updateAllPixelCosts();
}

void SubPixelSegmentation::updateAllPixelCosts() {
    const int w = map_.width();
    for (int i = 0; i <= map_.maxId(); ++i) {
        for (int q : subs_[i].band) {
            if (map_.idOf(q) != i) continue;
            cost_[q] = computePixelCost(q % w, q / w, -1, -1);
        }
    }
}

double SubPixelSegmentation::computeInteriorCost(int seg, int colorSeg) {
    const SubPixelSegment& s = subs_.at(seg);
    if (s.nInterior < 1) return 0.0;
    if (colorSeg < 0) colorSeg = seg;
    const float* color = colorOf(colorSeg);
    const double n = static_cast<double>(s.nInterior);
    double mean[4], var[4], t[4];
    for (int k = 0; k < 4; ++k) {
        mean[k] = static_cast<double>(s.interiorSum[k]) / n;
        var[k] = static_cast<double>(s.interiorSumSq[k]) - mean[k] * mean[k] * n;
    }
    for (int k = 0; k < 3; ++k)  // the alpha channel's variance is not clamped
        if (!(var[k] > 0.0)) var[k] = 0.0;
    for (int k = 0; k < 4; ++k) {
        const double d = static_cast<double>(color[k]) - mean[k];
        t[k] = d * d * n + var[k];
    }
    double sum = t[0] + 0.0;
    sum = t[1] + sum;
    sum = sum + t[2];
    sum = t[3] + sum;
    return sum;
}

double SubPixelSegmentation::computeBoundaryCost(bool useCached, bool store, int seg, int other) {
    SubPixelSegment& s = subs_.at(seg);
    double sum = 0.0;
    const int w = image_.w;
    for (int q : s.band) {
        if (useCached) {
            sum = sum + static_cast<double>(cost_[q]);
        } else {
            const float c = computePixelCost(q % w, q / w, other, seg);
            if (store) cost_[q] = c;
            sum = sum + static_cast<double>(c);
        }
    }
    return sum;
}

bool SubPixelSegmentation::isInterior(int id, int pixel) {
    const int w = map_.width(), h = map_.height();
    const int x = pixel % w, y = pixel / w;
    const int own = map_.idOf(x, y);
    if (own != id) return false;
    const int nn = use24_ ? 24 : 8;
    for (int k = 0; k < nn; ++k) {
        const int nx = x + (use24_ ? kN24i[k] : kN8i[k]);
        const int ny = y + (use24_ ? kN24j[k] : kN8j[k]);
        if (nx < 0 || nx >= w || ny < 0 || ny >= h) continue;
        if (map_.idOf(nx, ny) != own) return false;
    }
    return true;
}

void SubPixelSegmentation::mergeSegments(int into, int from) {
    SubPixelSegment& B = subs_.at(from);
    SubPixelSegment& A = subs_.at(into);
    map_.redirect(into, from);
    A.n += B.n;
    for (int k = 0; k < 4; ++k) A.sum[k] = A.sum[k] + B.sum[k];
    for (int k = 0; k < 4; ++k) B.sum[k] = 0.0f;
    B.n = 0;
    A.nInterior += B.nInterior;
    for (int k = 0; k < 4; ++k) A.sumSq[k] = B.sumSq[k] + A.sumSq[k];
    for (int k = 0; k < 4; ++k) A.interiorSum[k] = A.interiorSum[k] + B.interiorSum[k];
    for (int k = 0; k < 4; ++k) A.interiorSumSq[k] = B.interiorSumSq[k] + A.interiorSumSq[k];
    B.nInterior = 0;
    for (int k = 0; k < 4; ++k) B.interiorSum[k] = B.interiorSumSq[k] = 0.0f;
    computeBoundaryCost(false, true, from, into);
    for (int q : B.band) {
        auto [it, inserted] = A.band.insert(q);
        if (inserted || !isInterior(into, q)) continue;
        const uint8_t* px = image_.at(q);
        const float c[4] = {px[0] * kInv255f, px[1] * kInv255f, px[2] * kInv255f, px[3] * kInv255f};
        ++A.nInterior;
        for (int k = 0; k < 4; ++k) A.interiorSum[k] = A.interiorSum[k] + c[k];
        for (int k = 0; k < 4; ++k) A.interiorSumSq[k] = A.interiorSumSq[k] + c[k] * c[k];
        A.band.erase(it);
    }
    B.band.clear();
}

float SubPixelSegmentation::computePixelCost(int x, int y, int replaceTo, int replaceFrom) {
    double wOwn = 0.0;
    const Color4f model = wf_.computePixelF(x, y, replaceTo, replaceFrom, &wOwn, use24_);
    assert(x >= 0 && y >= 0 && x < image_.w && y < image_.h);
    const uint8_t* px = image_.at(y * image_.w + x);
    float r = 0.0f;
    float d[4];
    for (int k = 0; k < 4; ++k) {
        const float t = static_cast<float>(px[k]) * kInv255f - model.c[k];
        d[k] = t * t;
    }
    r = d[0] + r;
    r = d[1] + r;
    r = r + d[2];
    r = r + d[3];
    const double a = 1.0 - wOwn;
    double cost = a > 0.0 ? static_cast<double>(p_.self_weight_penalty_low) * a : 0.0;
    const double b = static_cast<double>(p_.self_weight_eps) - wOwn;
    if (b > 0.0) cost = cost + b * static_cast<double>(p_.self_weight_penalty_high);
    cost = cost + static_cast<double>(r);
    return static_cast<float>(cost);
}

}  // namespace isvg
