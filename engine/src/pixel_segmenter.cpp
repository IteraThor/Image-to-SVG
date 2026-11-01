#include "isvg/pixel_segmenter.h"

#include <algorithm>
#include <array>
#include <vector>

#include <cassert>

#include "isvg/math_util.h"

namespace isvg {
namespace {

// neighbors4 / neighbors8 tables (x offsets i, y offsets j): R, D, L, U, DL, UR, UL, DR.
constexpr int kN8i[8] = {1, 0, -1, 0, -1, 1, -1, 1};
constexpr int kN8j[8] = {0, 1, 0, -1, 1, -1, -1, 1};

// Sequential channel sum ((((0+c0)+c1)+c2)+c3), as the SSE code accumulates it.
inline float seqSum4(const float v[4]) {
    float s = 0.0f;
    s = v[0] + s;
    s = s + v[1];
    s = s + v[2];
    s = s + v[3];
    return s;
}

inline float sqDist(const uint8_t* px, const float mean[4]) {
    float d[4];
    for (int k = 0; k < 4; ++k) {
        const float t = static_cast<float>(px[k]) - mean[k];
        d[k] = t * t;
    }
    return seqSum4(d);
}

// QuadTree::computeVarianceAndSetIfTooHigh / setId / variance helper
//. Blocks whose colour variance is <= threshold get a fresh id; otherwise each
// low-variance child does (high-variance children were handled during their own recursion).
struct QuadTree {
    double threshold;
    PixelMap& map;
    const Image& img;

    void setId(int x, int y, int size, int& id) {
        for (int yy = y; yy < y + size; ++yy)
            for (int xx = x; xx < x + size; ++xx)
                if (xx >= 0 && xx < img.w && yy >= 0 && yy < img.h) map.raw(yy * img.w + xx) = id;
        ++id;
    }

    static double variance(const double s[4], const double q[4], int n) {
        double sumq = 0.0;
        sumq = sumq + q[0];
        sumq = sumq + q[1];
        sumq = sumq + q[2];
        sumq = sumq + q[3];
        double ss = 0.0;
        ss = s[0] * s[0] + ss;
        ss = ss + s[1] * s[1];
        ss = ss + s[2] * s[2];
        ss = s[3] * s[3] + ss;
        const double dn = static_cast<double>(n);
        return (sumq - ss / dn) / dn;
    }

    void compute(int x, int y, int size, int& id, double sum[4], double sumsq[4], int& n, double& var,
                 bool setIfLow) {
        if (x < 0 || x >= img.w || y < 0 || y >= img.h) {
            for (int k = 0; k < 4; ++k) sum[k] = sumsq[k] = 0.0;
            n = 0;
            var = 0.0;
            return;
        }
        if (size == 1) {
            const uint8_t* p = img.at(y * img.w + x);
            for (int k = 0; k < 4; ++k) {
                const double c = static_cast<double>(p[k]) * 0.00392156862745098;
                sum[k] = c;
                sumsq[k] = c * c;
            }
            n = 1;
            var = 0.0;
            return;
        }
        const int half = size / 2;
        const int cx[4] = {x, x + half, x, x + half};
        const int cy[4] = {y, y, y + half, y + half};
        for (int k = 0; k < 4; ++k) sum[k] = sumsq[k] = 0.0;
        n = 0;
        double childVar[4];
        for (int c = 0; c < 4; ++c) {
            double cs[4], cq[4];
            int cn;
            compute(cx[c], cy[c], half, id, cs, cq, cn, childVar[c], false);
            n += cn;
            for (int k = 0; k < 4; ++k) {
                sum[k] = sum[k] + cs[k];
                sumsq[k] = sumsq[k] + cq[k];
            }
        }
        var = n != 0 ? variance(sum, sumsq, n) : 0.0;
        if (var <= threshold) {
            if (setIfLow) setId(x, y, size, id);
        } else {
            for (int c = 0; c < 4; ++c)
                if (childVar[c] <= threshold) setId(cx[c], cy[c], half, id);
        }
    }
};

}  // namespace

void PixelSegmenter::initialize() {
    seg_.initialize();
    lambda_ = p_.lambda_initial;
    w_ = seg_.image().w;
    h_ = seg_.image().h;
    const float r = math::logf(p_.lambda_final / p_.lambda_initial);
    const float denom = static_cast<float>(p_.max_iterations) * p_.rise_fraction;
    rate_ = static_cast<float>(math::exp(static_cast<double>(r / denom)));
}

void PixelSegmenter::preSegmentation() {
    QuadTree qt{static_cast<double>(p_.lambda_pre), seg_.map(), seg_.image()};
    int nextId = w_ * h_;
    for (int y = 0; y < h_; y += 16) {
        for (int x = 0; x < w_; x += 16) {
            double s[4], q[4], var;
            int n;
            int id = nextId;
            qt.compute(x, y, 16, id, s, q, n, var, true);
            nextId = id;
        }
    }
}

float PixelSegmenter::computeMergeMetric(int a, int b) {
    const PixelSegment& A = seg_.segment(a);
    const PixelSegment& B = seg_.segment(b);
    float a2[4], b2[4], ab2[4];
    for (int k = 0; k < 4; ++k) {
        a2[k] = A.sum[k] * A.sum[k];
        b2[k] = B.sum[k] * B.sum[k];
        const float t = A.sum[k] + B.sum[k];
        ab2[k] = t * t;
    }
    const float nA = static_cast<float>(A.n);
    const float nB = static_cast<float>(B.n);
    const float nAB = static_cast<float>(A.n + B.n);
    float m = seqSum4(a2) / nA;
    m = m + seqSum4(b2) / nB;
    m = m - seqSum4(ab2) / nAB;
    if (p_.size_penalty_weight > 0.0f) {
        const double pw = static_cast<double>(p_.size_penalty_power);
        const float pa = static_cast<float>(math::pow(static_cast<double>(nA), pw));
        const float pb = static_cast<float>(math::pow(static_cast<double>(nB), pw));
        const float pab = static_cast<float>(math::pow(static_cast<double>(nAB), pw));
        float t = 1.0f / pa;
        t = t + 1.0f / pb;
        t = t - 1.0f / pab;
        m = m - t * p_.size_penalty_weight;
    }
    return m;
}

void PixelSegmenter::swapPixel(int p, int from, int to) {
    PixelSegment& F = seg_.segment(from);
    PixelSegment& T = seg_.segment(to);
    seg_.map().raw(p) = to;
    const uint8_t* px = seg_.image().at(p);
    const float c[4] = {px[0] * kInv255f, px[1] * kInv255f, px[2] * kInv255f, px[3] * kInv255f};
    T.addPixelF(c);
    --F.n;
    assert(F.n >= 0);
    for (int k = 0; k < 4; ++k) F.sum[k] = F.sum[k] - c[k];
    if (F.n == 0) seg_.map().decrementLive();
}

void PixelSegmenter::mergeSegments(int into, int from) {
    seg_.map().redirect(into, from);
    PixelSegment& I = seg_.segment(into);
    PixelSegment& F = seg_.segment(from);
    I.n += F.n;
    for (int k = 0; k < 4; ++k) I.sum[k] = F.sum[k] + I.sum[k];
    F.n = 0;
    for (int k = 0; k < 4; ++k) F.sum[k] = 0.0f;
}

int PixelSegmenter::considerMergeAndSwap(int p1, int p2, int a, int b) {
    PixelSegment& A = seg_.segment(a);
    PixelSegment& B = seg_.segment(b);
    const float m = computeMergeMetric(a, b);
    bool merge = lambda_ > m;
    if (!merge && lambda_ == p_.lambda_final)
        merge = A.n < p_.min_num_pixels || B.n < p_.min_num_pixels;
    if (merge) {
        mergeSegments(a, b);
        return 1;
    }
    float meanB[4], meanA[4];
    const float fB = 255.0f / static_cast<float>(B.n);
    for (int k = 0; k < 4; ++k) meanB[k] = B.sum[k] * fB;
    const uint8_t* px1 = seg_.image().at(p1);
    const float d1B = sqDist(px1, meanB);
    const float fA = 255.0f / static_cast<float>(A.n);
    for (int k = 0; k < 4; ++k) meanA[k] = A.sum[k] * fA;
    const float d1A = sqDist(px1, meanA);
    if (d1A > d1B) {
        swapPixel(p1, a, b);
        return 2;
    }
    const uint8_t* px2 = seg_.image().at(p2);
    const float d2A = sqDist(px2, meanA);
    const float d2B = sqDist(px2, meanB);
    if (d2B > d2A) {
        swapPixel(p2, b, a);
        return 3;
    }
    return 0;
}

void PixelSegmenter::mergeAndSwapLoop(MergeConnectedness mode) {
    PixelMap& map = seg_.map();
    int idx = 0;
    for (int y = 0; y < h_; ++y) {
        for (int x = 0; x < w_; ++x, ++idx) {
            int id = map.idOf(idx);
            const int nk = mode == MergeConnectedness::Eight ? 8 : 2;
            for (int k = 0; k < nk; ++k) {
                const int nx = x + kN8i[k], ny = y + kN8j[k];
                if (nx >= w_ || ny >= h_) continue;
                if (mode == MergeConnectedness::Eight && (nx < 0 || ny < 0)) continue;
                const int q = idx + kN8i[k] + kN8j[k] * w_;
                const int nid = map.idOf(q);
                if (nid == id) continue;
                if (mode == MergeConnectedness::MinPixels) {
                    if (seg_.segment(id).n < p_.min_num_pixels || seg_.segment(nid).n < p_.min_num_pixels) {
                        mergeSegments(id, nid);
                        id = map.idOf(idx);
                    }
                } else if (considerMergeAndSwap(idx, q, id, nid) != 0) {
                    id = map.idOf(idx);
                }
            }
        }
    }
}

void PixelSegmenter::stepLambda() {
    const float next = lambda_ * rate_;
    lambda_ = next <= p_.lambda_final ? next : p_.lambda_final;
}

void PixelSegmenter::overSegmentation(bool diagonal) {
    preSegmentation();
    seg_.updateStatistics();
    const int savedMin = p_.min_num_pixels;
    if (diagonal) p_.min_num_pixels = 0;
    for (int i = 1; i <= p_.max_iterations; ++i) {
        mergeAndSwapLoop(MergeConnectedness::Four);
        if (i % p_.num_iter_btwn_renum == 0) seg_.updateStatistics();
        stepLambda();
    }
    if (diagonal) {
        mergeAndSwapLoop(MergeConnectedness::Eight);
        seg_.performDiagonalExtrusions();
        p_.min_num_pixels = savedMin;
    }
    seg_.updateStatistics();
    if (p_.min_num_pixels > 0) {
        mergeAndSwapLoop(MergeConnectedness::MinPixels);
        seg_.updateStatistics();
    }
    seg_.compact();
}

namespace {
// CVTTSS2SI: NaN and out-of-range values give 0x80000000 (low byte 0).
inline uint8_t truncToByte(float v) {
    if (!(v >= -2147483648.0f && v < 2147483648.0f)) return 0;
    return static_cast<uint8_t>(static_cast<int>(v));
}
}  // namespace

// Small-segment colour reassignment compares costs against -10, which a non-negative cost never
// undercuts, so it never changes a colour and is omitted.
void PixelSegmenter::clusterColorsAndGeneratePalette() {
    PixelMap& map = seg_.map();
    const int n = map.maxId() + 1;
    std::vector<std::array<uint8_t, 4>> colors(std::max(n, 0));
    for (int i = 0; i < n; ++i) {
        const PixelSegment& sg = seg_.segment(i);
        const float scale = 255.0f / static_cast<float>(sg.n);
        for (int k = 0; k < 4; ++k) colors[i][k] = truncToByte(sg.sum[k] * scale);
    }
    Palette& pal = *seg_.palette();
    pal.clearMap();
    for (int i = 0; i < n; ++i) seg_.segment(i).paletteIndex = pal.tryAddColorToMap(colors[i].data());
    pal.setNewPaletteFromMap();
    const int total = w_ * h_;
    for (int idx = 0; idx < total; ++idx) map.raw(idx) = seg_.segment(map.idOf(idx)).paletteIndex;
    seg_.updateMaxId(pal.size() - 1);
    for (int i = 0; i <= map.maxId(); ++i) seg_.segment(i).paletteIndex = i;
}

void PixelSegmenter::enforceMinNumPixels() {
    if (p_.min_num_pixels < 1) return;
    mergeAndSwapLoop(MergeConnectedness::MinPixels);
    seg_.updateStatistics();
}

}  // namespace isvg
