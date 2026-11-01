// SubPixelSegmenter::updateAllFeatureVectors and flagShadowSegments.
//
// Feature index map (unlisted entries stay 0):
//   0 Σ s1          1 Σ s2          2 √f0          3 √f1          4 √(f0/n)       5 √(f1/n)
//   6 Σ √s1         7 Σ s2 (again)  8 n            9 interior n  10 #unmasked px  11 #neighbour ids
//  12-18 neighbour counts by interior/total size   19/20 brighter / not brighter neighbours
//  21 luminance    22 Σ colour variance            23 interior/band ratio         24 bbox aspect
//  25-32 per-pixel averages of 0,1,2,3,4,5,6,7     39 min colour distance to big segments
//  41 Σ squared residual from the principal colour line
#include <cmath>
#include <cstdlib>
#include <list>
#include <set>

#include "isvg/subpixel_segmenter.h"

namespace isvg {
namespace {

inline float seqSum4(const float v[4]) {
    float s = 0.0f;
    s = v[0] + s;
    s = s + v[1];
    s = s + v[2];
    s = s + v[3];
    return s;
}

inline double luminance(const float c[4]) {
    double l = static_cast<double>(c[2]) * 0.3;
    l = l + static_cast<double>(c[1]) * 0.59;
    l = l + static_cast<double>(c[0]) * 0.11;
    return l;
}

}  // namespace

void SubPixelSegmenter::updateAllFeatureVectors() {
    seg_.resetPaletteIndices();
    seg_.rebuildSubPixelSegments(0x1b, false, false);
    const Image& img = seg_.image();
    const SignalImages sig = computeSignals(img);
    PixelMap& map = seg_.map();
    const std::vector<uint8_t> mask = boundaryMask8(map, 2);

    const int numSeg = map.maxId() + 1;
    std::vector<std::array<float, 16>> outer(numSeg);
    std::vector<std::array<float, 4>> axis(numSeg);
    for (auto& o : outer) o.fill(0.0f);

    // Phase 1: per-pixel sums.
    const int w = img.w, h = img.h, mw = map.width();
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            const int id = map.raw(y * mw + x);
            FeatureVector& f = features_.at(id);
            const float s1 = sig.s1[static_cast<size_t>(y) * w + x];
            f[0] = static_cast<double>(s1) + f[0];
            const float s2 = sig.s2[static_cast<size_t>(y) * w + x];
            f[1] = static_cast<double>(s2) + f[1];
            f[6] = static_cast<double>(std::sqrt(s1)) + f[6];
            f[7] = static_cast<double>(s2) + f[7];
            if (mask[static_cast<size_t>(y) * mw + x] == 0) f[10] = f[10] + 1.0;
            const uint8_t* p = img.at(y * w + x);
            float c[4];
            for (int k = 0; k < 4; ++k) c[k] = static_cast<float>(p[k]) * kInv255f;
            auto& o = outer[id];
            for (int r = 0; r < 4; ++r)
                for (int k = 0; k < 4; ++k) o[r * 4 + k] = o[r * 4 + k] + c[r] * c[k];
        }
    }

    // Phase 2: colour covariance and its principal axis by power iteration.
    for (int i = 0; i < numSeg; ++i) {
        SubPixelSegment& s = seg_.sub(i);
        features_.at(i)[41] = 0.0;
        const float n = static_cast<float>(s.n);
        auto& C = outer[i];
        for (float& v : C) v = v / n;
        float mean[4];
        for (int k = 0; k < 4; ++k) mean[k] = s.sum[k] / n;
        for (int r = 0; r < 4; ++r)
            for (int k = 0; k < 4; ++k) C[r * 4 + k] = C[r * 4 + k] - mean[r] * mean[k];
        float v[4] = {1, 1, 1, 1};
        auto mul = [&](const float in[4], float out[4]) {
            for (int k = 0; k < 4; ++k) out[k] = 0.0f;
            for (int r = 0; r < 4; ++r)
                for (int k = 0; k < 4; ++k) out[k] = out[k] + C[r * 4 + k] * in[r];
        };
        float wv[4];
        mul(v, wv);
        float prevLam = 0.0f;
        float keep[4] = {wv[0], wv[1], wv[2], wv[3]};  // XMM11 (first lane is read back)
        for (int it = 10;;) {
            float x[4] = {wv[0], wv[1], wv[2], wv[3]};
            float sq[4];
            for (int k = 0; k < 4; ++k) sq[k] = x[k] * x[k];
            const float nrm = std::sqrt(seqSum4(sq));
            if (nrm > 0.0f) {
                for (int k = 0; k < 4; ++k) x[k] = x[k] / nrm;
                for (int k = 0; k < 4; ++k) keep[k] = x[k];
            }
            v[0] = keep[0];
            v[1] = x[1];
            v[2] = x[2];
            v[3] = x[3];
            float nw[4];
            mul(v, nw);
            float pr[4];
            for (int k = 0; k < 4; ++k) pr[k] = v[k] * nw[k];
            const float lam = seqSum4(pr);
            const float diff = std::fabs(prevLam - lam);
            for (int k = 0; k < 4; ++k) wv[k] = nw[k];
            if (9.99999997e-07f > diff) break;
            if (--it == 0) break;
            for (int k = 0; k < 4; ++k) keep[k] = nw[k];
            prevLam = lam;
        }
        for (int k = 0; k < 4; ++k) axis[i][k] = v[k];
    }

    // Phase 3: squared residual of each pixel from its segment's principal colour line.
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            const int id = map.raw(y * mw + x);
            const SubPixelSegment& s = seg_.sub(id);
            const uint8_t* p = img.at(y * w + x);
            float c[4], mean[4], d[4];
            for (int k = 0; k < 4; ++k) c[k] = static_cast<float>(p[k]) * kInv255f;
            const float n = static_cast<float>(s.n);
            for (int k = 0; k < 4; ++k) mean[k] = s.sum[k] / n;
            const auto& v = axis[id];
            for (int k = 0; k < 4; ++k) d[k] = (c[k] - mean[k]) * v[k];
            const float t = seqSum4(d);
            float r[4];
            for (int k = 0; k < 4; ++k) {
                const float proj = t * v[k] + mean[k];
                const float e = c[k] - proj;
                r[k] = e * e;
            }
            FeatureVector& f = features_.at(id);
            f[41] = static_cast<double>(seqSum4(r)) + f[41];
        }
    }

    // Phase 4: derived per-segment features.
    for (int i = 0; i <= map.maxId(); ++i) {
        SubPixelSegment& s = seg_.sub(i);
        FeatureVector& f = features_.at(i);
        const double f0 = f[0], f1 = f[1];
        f[2] = std::sqrt(f0);
        f[3] = std::sqrt(f1);
        const double n = static_cast<double>(s.n);
        const double f0n = f0 / n, f1n = f1 / n;
        f[4] = std::sqrt(f0n);
        f[5] = std::sqrt(f1n);
        f[8] = n;
        f[26] = f1n;
        f[25] = f0n;
        f[39] = 1e100;
        f[9] = static_cast<double>(s.nInterior);
        f[27] = f[2] / n;
        f[31] = f[6] / n;
        f[28] = f[3] / n;
        f[32] = f[7] / n;
        f[29] = f[4] / n;
        f[30] = f[5] / n;
        for (int j = 0; j <= map.maxId(); ++j) {
            const SubPixelSegment& t = seg_.sub(j);
            if (t.nInterior <= 10) continue;
            float d[4];
            for (int k = 0; k < 4; ++k) {
                const float e = t.color[k] - s.color[k];
                d[k] = e * e;
            }
            const double dist = static_cast<double>(seqSum4(d));
            if (f[39] > dist) f[39] = dist;
        }
        float c[4];
        seg_.getColorF(i, c);
        const double lum = luminance(c);
        f[21] = lum;
        {
            const float nf = static_cast<float>(s.n);
            float var[4];
            for (int k = 0; k < 4; ++k) {
                const float m = s.sum[k] / nf;
                var[k] = s.sumSq[k] / nf - m * m;
            }
            f[22] = static_cast<double>(seqSum4(var));
        }
        f[23] = static_cast<double>(s.nInterior) / static_cast<double>(s.band.size());
        std::set<int> nbr;
        int maxX = 0, minX = img.w, maxY = 0, minY = img.h;
        for (int q : s.band) {
            nbr.insert(map.raw(q));
            const int x = q % mw, y = q / mw;
            if (maxX < x) maxX = x;
            if (!(minX <= x)) minX = x;
            if (maxY < y) maxY = y;
            if (!(minY <= y)) minY = y;
        }
        if (s.band.empty()) {
            f[11] = -1.0;
            f[24] = static_cast<double>(-minX) / static_cast<double>(-minY);
            continue;
        }
        f[11] = static_cast<double>(nbr.size()) - 1.0;
        f[24] = static_cast<double>(maxX - minX) / static_cast<double>(maxY - minY);
        for (int j : nbr) {
            if (j == i) continue;
            const SubPixelSegment& t = seg_.sub(j);
            const int nI = t.nInterior;
            if (nI == 0) {
                f[12] = f[12] + 1.0;
            } else if (nI > 4) {
                f[16] = f[16] + 1.0;
                if (nI > 9) {
                    f[17] = f[17] + 1.0;
                    if (nI > 14) f[18] = f[18] + 1.0;
                }
            }
            if (t.n > 4) {
                f[13] = f[13] + 1.0;
                if (t.n > 9) {
                    f[14] = f[14] + 1.0;
                    if (t.n > 14) f[15] = f[15] + 1.0;
                }
            }
            float cj[4];
            seg_.getColorF(j, cj);
            if (luminance(cj) > lum)
                f[19] = 1.0 + f[19];
            else
                f[20] = 1.0 + f[20];
        }
    }
}

namespace {

// SubPixelSegmentation::computePixelAntiAliasingSwapCost(InDirection) /:
// cheapest 4-neighbour of pixel p that lies in a non-shadow segment, by squared colour difference.
bool antiAliasingSwapCost(SubPixelSegmentation& seg, int p, float* cost, int* newId) {
    const Image& img = seg.image();
    PixelMap& map = seg.map();
    const int w = map.width();
    const int x = p % w, y = p / w;
    float c[4];
    const uint8_t* px = img.at(p);
    for (int k = 0; k < 4; ++k) c[k] = static_cast<float>(px[k]) * kInv255f;
    *cost = 1e10f;
    auto dir = [&](int nx, int ny) {
        if (nx < 0 || nx >= img.w || ny < 0 || ny >= img.h) return false;
        const int id = map.raw(ny * w + nx);
        if (seg.sub(id).flag) return false;
        const uint8_t* q = img.at(ny * img.w + nx);
        float d[4];
        for (int k = 0; k < 4; ++k) {
            const float t = static_cast<float>(q[k]) * kInv255f - c[k];
            d[k] = t * t;
        }
        const float s = seqSum4(d);
        if (!(*cost > s)) return false;
        *cost = s;
        *newId = id;
        return true;
    };
    const bool a = dir(x - 1, y);
    const bool b = dir(x + 1, y);
    const bool c2 = dir(x, y - 1);
    const bool d = dir(x, y + 1);
    return a || b || c2 || d;
}

}  // namespace

void SubPixelSegmenter::splitAndMergeFlaggedSegmentsFromBlurryImage(int configT) {
    PixelMap& map = seg_.map();
    const int numSeg = map.maxId() + 1;
    std::vector<int> flaggedIdx(numSeg, -1);
    int K = 0;
    for (int s = 0; s <= map.maxId(); ++s)
        if (seg_.sub(s).flag) flaggedIdx[s] = K++;
    if (K == 0) return;

    std::vector<std::list<int>> pixels(K);  // local_148: pixels of each flagged segment
    std::vector<int> segOf(K, 0);  // local_158
    std::vector<std::set<int>> neighbours(K);  // local_178: non-flagged neighbour ids
    std::vector<std::list<int>> outerBand(K);  // local_138: band pixels of other segments
    std::vector<float> before(K, 0.0f), after(K, 0.0f);

    seg_.rebuildSubPixelSegments(0x20, true, true);
    const std::vector<float>& cost = seg_.pixelCosts();
    for (int s = 0; s <= map.maxId(); ++s) {
        const int k = flaggedIdx.at(s);
        if (k < 0) continue;
        segOf[k] = s;
        // As compiled: the interior cost is taken for segment id `k` (the flagged index), not `s`.
        before[k] = static_cast<float>(seg_.computeInteriorCost(k, -1)) + before[k];
        for (int q : seg_.sub(s).band) {
            const int id = map.raw(q);
            before[k] = cost[q] + before[k];
            if (id == s) continue;
            outerBand[k].push_back(q);
            if (!seg_.sub(id).flag) neighbours[k].insert(id);
        }
    }
    const int w = map.width(), h = map.height();
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            const int k = flaggedIdx.at(map.raw(y * w + x));
            if (k >= 0) pixels[k].push_back(seg_.image().w * y + x);
        }

    // Peel each flagged segment, cheapest pixel first.
    for (int k = 0; k < K; ++k) {
        std::list<int> todo(pixels[k]);
        while (!todo.empty()) {
            auto best = todo.end();
            int bestId = -1;
            float bestCost = 0.0f;
            for (auto it = todo.begin(); it != todo.end(); ++it) {
                float c;
                int id;
                if (!antiAliasingSwapCost(seg_, *it, &c, &id)) continue;
                if (best == todo.end() || bestCost > c) {
                    best = it;
                    bestId = id;
                    bestCost = c;
                }
            }
            if (best == todo.end()) std::abort();  // invariant violated: the list would be corrupted
            map.raw(*best) = bestId;
            todo.erase(best);
        }
    }

    const PixelMap backup = map;
    seg_.renumber();
    seg_.initialize(0x24, false);
    seg_.rebuildSubPixelSegments(0x27, true, true);
    for (int k = 0; k < K; ++k) {
        for (int p : pixels[k]) after[k] = seg_.computePixelCost(p % w, p / w, -1, -1) + after[k];
        for (int p : outerBand[k]) after[k] = seg_.computePixelCost(p % w, p / w, -1, -1) + after[k];
    }
    map = backup;

    for (int k = 0; k < K; ++k) {
        const int n = static_cast<int>(pixels[k].size()) + static_cast<int>(outerBand[k].size());
        FeatureVector& f = features_.at(segOf[k]);
        const float b = before[k], a = after[k];
        const float nf = static_cast<float>(n);
        f[33] = static_cast<double>(b);
        f[34] = static_cast<double>(a);
        f[37] = static_cast<double>(a - b);
        f[35] = static_cast<double>(b / nf);
        f[36] = static_cast<double>(a / nf);
        f[38] = static_cast<double>(a - b) / static_cast<double>(n);
        if (isProbablyNotFalsePositive(f, configT)) continue;
        for (int p : pixels[k]) map.raw(p) = segOf[k];
    }
    seg_.renumber();
    seg_.initialize(0x2b, false);
}

void SubPixelSegmenter::flagShadowSegments() {
    const int numSeg = seg_.map().maxId() + 1;
    features_.assign(numSeg, FeatureVector{});
    updateAllFeatureVectors();
    for (int i = 0; i <= seg_.map().maxId(); ++i)
        if (isProbablyShadowingSegment(features_.at(i))) seg_.sub(i).flag = true;
}

}  // namespace isvg
