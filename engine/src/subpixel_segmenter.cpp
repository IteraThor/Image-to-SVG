#include "isvg/subpixel_segmenter.h"

#include <cstdint>
#include <list>
#include <set>
#include <utility>
#include <vector>

namespace isvg {

void SubPixelSegmenter::clusterColors(int /*stage*/) {
    PixelMap& map = seg_.map();
    std::set<std::pair<int, int>> bySize;  // (interior pixel count, id)
    for (int i = 0; i <= map.maxId(); ++i) {
        const int nInt = seg_.sub(i).nInterior;
        if (nInt > 2) bySize.emplace(nInt, i);
    }
    // populateBigSegmentIds: up to 100 ids, largest first.
    std::vector<int> big;
    for (auto it = bySize.rbegin(); it != bySize.rend() && big.size() < 100; ++it) big.push_back(it->second);

    const int maxId = map.maxId();
    for (int i = 0; i <= maxId; ++i) {
        if (seg_.sub(i).nInterior >= 3) continue;
        double best = seg_.computeSegmentCost(i, -1) - static_cast<double>(p_.color_cluster_margin);
        int bestId = i;
        for (int k : big) {
            const double interior = seg_.computeInteriorCost(i, k);
            if (best <= interior) continue;
            const double boundary = seg_.computeBoundaryCostUpTo(best - interior, i, k);
            if (interior + boundary < best) {
                best = interior + boundary;
                bestId = k;
            }
        }
        if (bestId == i) continue;
        SubPixelSegment& s = seg_.sub(i);
        const SubPixelSegment& t = seg_.sub(bestId);
        for (int c = 0; c < 4; ++c) s.color[c] = t.color[c];
        s.paletteIndex = t.paletteIndex;
    }
}

void SubPixelSegmenter::finalize() {
    evaluated_.clear();
    seg_.finalize();
}

void SubPixelSegmenter::execute() {
    executeMergePhases();
    splitAndMergeBeaches();
    seg_.performDiagonalExtrusions();
    seg_.rebuildSubPixelSegments(0x12, useInterior_, false);
    for (int i = 0;; ++i)
        if (mergeLoop(INT_MAX, 4) == 0 || p_.max_iterations <= i) break;
    seg_.rebuildSubPixelSegments(0x16, useInterior_, false);
}

// Segments with no interior pixels are relabelled pixel by pixel. As compiled, a segment is only
// processed if none of its band pixels belongs to another segment, and the relabel target
// (WeightsFinder::getHeaviestSegmentId) is always the pixel's own segment; ported as is.
void SubPixelSegmenter::splitAndMergeBeaches() {
    PixelMap& map = seg_.map();
    const int maxId = map.maxId();
    const int w = map.width();
    std::list<int> pixels;
    for (int i = 0; i <= maxId; ++i) {
        SubPixelSegment& s = seg_.sub(i);
        if (!(map.redirectOf(i) < 0 && s.nInterior == 0)) continue;
        pixels.clear();
        bool skip = false;
        for (int q : s.band) {
            if (map.idOf(q) == i) {
                pixels.push_back(q);
            } else if (s.nInterior == 0) {
                skip = true;
                break;
            }
        }
        if (skip) continue;
        if (!pixels.empty()) {
            int count = 0, adjacent = 0;
            for (int q : pixels) {
                const int x = q % w, y = q / w;
                if (map.idOf(x, y - 1) == i || map.idOf(x, y + 1) == i) ++adjacent;
                if (map.idOf(x - 1, y) == i || map.idOf(x + 1, y) == i) ++adjacent;
                ++count;
            }
            if (static_cast<double>(adjacent) / static_cast<double>(count) > 1.2) continue;
        }
        while (!pixels.empty()) {
            int changes = 0;
            for (auto it = pixels.begin(); it != pixels.end();) {
                const int q = *it;
                const float old = seg_.pixelCosts()[q];
                const float c = seg_.computePixelCost(q % w, q / w, -1, i);
                if (c - old < p_.beach_lambda_f) {
                    map.raw(q) = seg_.weights().heaviestSegmentId();
                    it = pixels.erase(it);
                    ++changes;
                    seg_.pixelCosts()[q] = c;
                } else {
                    ++it;
                }
            }
            if (changes == 0) break;
        }
    }
}

double SubPixelSegmenter::computeMergeMetric(int nb, int seg) {
    const double own = seg_.computeSegmentCost(seg, -1);
    double d = seg_.computeInteriorCost(seg, nb) - own;
    if (d <= static_cast<double>(p_.lambda_f)) d = d + seg_.computeBoundaryCost(false, false, seg, nb);
    return d;
}

int SubPixelSegmenter::mergeLoop(int maxNeighborSize, int maxSize) {
    PixelMap& map = seg_.map();
    int merges = 0;
    for (int i = 0; i <= map.maxId(); ++i) {
        SubPixelSegment& s = seg_.sub(i);
        const int size = useInterior_ ? s.nInterior : s.n;
        if (!(size <= maxSize && map.redirectOf(i) < 0)) continue;
        int best = -1;
        double bestMetric = static_cast<double>(p_.lambda_f);
        for (int q : s.band) {
            const int nb = map.idOf(q);
            if (nb == i) continue;
            const int key = static_cast<int>(static_cast<uint32_t>(nb) | (static_cast<uint32_t>(i) << 16));
            auto it = evaluated_.find(key);
            if (it != evaluated_.end() && it->second) continue;
            if (maxNeighborSize < seg_.sub(nb).n) continue;
            const double m = computeMergeMetric(nb, i);
            if (m < bestMetric) {
                bestMetric = m;
                best = nb;
            }
            evaluated_[key] = true;
        }
        if (best >= 0) {
            seg_.mergeSegments(best, i);
            ++merges;
        }
    }
    for (auto& kv : evaluated_) kv.second = false;
    return merges;
}

void SubPixelSegmenter::executeMergePhases() {
    useInterior_ = false;
    seg_.initialize(0, true);
    const int maxIt = p_.max_iterations;
    for (int i = 0;; ++i)
        if (mergeLoop(i + 1, i + 1) == 0 || maxIt <= i) break;
    seg_.rebuildSubPixelSegments(4, useInterior_, false);
    for (int i = 0;; ++i)
        if (mergeLoop(INT_MAX, i + 1) == 0 || maxIt <= i) break;
    for (int i = 0;; ++i)
        if (mergeLoop(INT_MAX, INT_MAX) == 0 || maxIt <= i) break;
    if (p_.use_connected_24_iteration != 0) {
        useInterior_ = true;
        seg_.rebuildSubPixelSegments(9, true, false);
        for (int i = 0;; ++i)
            if (mergeLoop(INT_MAX, maxIt) == 0 || maxIt <= i) break;
        useInterior_ = false;
        seg_.rebuildSubPixelSegments(0xd, false, false);
    }
}

}  // namespace isvg
