#include "isvg/boundary_tracer.h"

#include <algorithm>
#include <cmath>
#include <list>

namespace isvg {

namespace {

// Pixel corners relative to a pixel (BoundaryTracer.
constexpr int kCornerDx[4] = {0, 1, 1, 0};
constexpr int kCornerDy[4] = {0, 0, 1, 1};
// Pixels around a node, counter-clockwise from the top-left.
constexpr int kPixelDx[4] = {-1, -1, 0, 0};
constexpr int kPixelDy[4] = {-1, 0, 0, -1};
// Step directions from a node: left, down, right, up. Step k runs between
// pixel k and pixel k+1.
constexpr int kStepDx[4] = {-1, 0, 1, 0};
constexpr int kStepDy[4] = {0, 1, 0, -1};

}  // namespace

void buildRegions(PixelSegmentation& seg, std::vector<Region>& regions) {
    const int maxId = seg.map().maxId();
    regions.assign(maxId + 1, Region{});
    for (int i = 0; i <= maxId; ++i) {
        Region& r = regions[i];
        float c[4];
        seg.getColorF(i, c);
        for (int k = 0; k < 4; ++k) r.color[k] = static_cast<uint8_t>(static_cast<int>(c[k] * 255.0f + 0.5f));
        r.n = seg.segment(i).n;
        r.paletteIndex = seg.segment(i).paletteIndex;
    }
}

int BoundaryTracer::enumerateBoundaryNodes() {
    // Count, per corner, the region boundaries that touch it. Only edges between pixels inside
    // [0,w-1)×[0,h-1) are scanned, plus a special case for the bottom-right pixel.
    std::fill(nodeIndex_.begin(), nodeIndex_.end(), 0);
    for (int y = 0; y < h_ - 1; ++y) {
        for (int x = 0; x < w_ - 1; ++x) {
            if (id(x, y) != id(x, y + 1)) {
                ++nodeIndex(x, y + 1);
                ++nodeIndex(x + 1, y + 1);
            }
            if (id(x, y) != id(x + 1, y)) {
                ++nodeIndex(x + 1, y);
                ++nodeIndex(x + 1, y + 1);
            }
        }
    }
    if (id(w_ - 1, h_ - 1) != id(w_ - 2, h_ - 1) || id(w_ - 1, h_ - 1) != id(w_ - 1, h_ - 2))
        ++nodeIndex(w_ - 1, h_ - 1);

    // Border corners always become nodes (left column, right column, top row, bottom row), then
    // every interior corner touched by a boundary.
    int n = 0;
    for (int y = 0; y <= h_; ++y) nodeIndex(0, y) = n++;
    for (int y = 0; y <= h_; ++y) nodeIndex(w_, y) = n++;
    for (int x = 1; x < w_; ++x) nodeIndex(x, 0) = n++;
    for (int x = 1; x < w_; ++x) nodeIndex(x, h_) = n++;
    for (int y = 1; y < h_; ++y)
        for (int x = 1; x < w_; ++x) {
            int& v = nodeIndex(x, y);
            v = v > 0 ? n++ : -1;
        }
    return n;
}

void BoundaryTracer::traceBoundary(int seg, int x, int y) {
    // Walk the boundary of `seg` keeping it on the same side, one pixel edge per step, until a
    // corner already claimed by this region is reached. Corners with a single way on are claimed
    // as they are passed; at a choice the step continuing along the last pixel of `seg` wins.
    std::list<int> nodeList, nbList;
    int prevX = -1, prevY = -1;  // pixel of `seg` along the last unambiguous step
    int count = 0;
    // Candidate buffers live across iterations (stale entries are read back).
    int cx[4] = {0, 0, 0, 0}, cy[4] = {0, 0, 0, 0}, cdir[4] = {0, 0, 0, 0}, cnb[4] = {0, 0, 0, 0};
    int nb = 0;
    for (;;) {
        nodeList.push_back(nodeIndex(x, y));
        int nc = 0;
        for (int k = 0; k < 4; ++k) {
            const int nx = x + kStepDx[k], ny = y + kStepDy[k];
            if (nx < 0 || ny < 0 || nx > w_ || ny > h_) continue;
            const int v = visited(nx, ny);
            const int nn = nodeIndex(nx, ny);
            const int ax = x + kPixelDx[k], ay = y + kPixelDy[k];
            const int bx = x + kPixelDx[(k + 1) & 3], by = y + kPixelDy[(k + 1) & 3];
            const int idA = (ax >= 0 && ax < w_ && ay >= 0 && ay < h_) ? id(ax, ay) : -1;
            const int idB = (bx >= 0 && bx < w_ && by >= 0 && by < h_) ? id(bx, by) : -1;
            if (idB != idA && idA == seg && nn >= 0) {
                cnb[nc] = idB;
                if (v != seg) {
                    cx[nc] = nx;
                    cy[nc] = ny;
                    cdir[nc] = k;
                    ++nc;
                }
            }
        }
        bool takeFirst = true;
        if (nc <= 1) {
            visited(x, y) = seg;
        } else if (prevX != -1 && prevY != -1) {
            takeFirst = false;
            for (int i = 0; i < nc; ++i) {
                const int d = cdir[i];
                if (x + kPixelDx[d] == prevX && y + kPixelDy[d] == prevY) {
                    x = cx[i];
                    y = cy[i];
                    nb = cnb[i];
                    break;
                }
            }
        }
        if (takeFirst) {
            const int d = cdir[0];
            prevX = x + kPixelDx[d];
            prevY = y + kPixelDy[d];
            x = cx[0];
            y = cy[0];
            nb = cnb[0];
        }
        nbList.push_back(nb);
        if (visited(x, y) == seg) break;
        ++count;
    }

    Region& r = regions_[seg];
    const int n = count + 1;
    r.nodes.assign(nodeList.begin(), nodeList.end());
    r.elems.assign(n, ContourElement{});
    auto it = nbList.begin();
    for (int i = 0; i < n; ++i, ++it) r.elems[i].neighbor = *it;
    nbList.sort();
    nbList.unique();
    nbList.remove(-1);
    r.neighbors.assign(nbList.begin(), nbList.end());
}

namespace {

// Staircase patterns (BoundaryTracer ctor, data: length, step bits
// (first step = highest bit) and the keep-mask over the len+1 nodes of the matched run.
struct StairPattern {
    int len, bits, keep;
};
constexpr StairPattern kStairPatterns[6] = {
    {5, 0b00000, 0b100001},
    {9, 0b010101010, 0b1001001001},
    {9, 0b010010010, 0b1001001001},
    {13, 0b1000100010001, 0b10000100100001},
    {10, 0b0010000100, 0b10000100001},
    {10, 0b0101001010, 0b10000100001},
};

// BoundaryTracer::patternMatch: first pattern no longer than `len` whose bits equal the
// low bits of the y- or x-step string.
int patternMatch(int bitsX, int bitsY, int len, int* keep) {
    for (const StairPattern& p : kStairPatterns) {
        if (p.len > len) continue;
        const int mask = (1 << p.len) - 1;
        if ((mask & bitsY) == p.bits || (mask & bitsX) == p.bits) {
            *keep = p.keep;
            return p.len;
        }
    }
    return 0;
}

}  // namespace

void BoundaryTracer::removeRedundantNodes(bool antiAliased) {
    std::vector<Node>& nodes = image_.nodes;
    // Border nodes that are not fixed only run along the image edge: drop them. The image
    // corners (first/last node of the left and right columns) are skipped.
    int removed[4] = {0, 0, 0, 0};
    int idx = 0;
    for (int s = 0; s < 4; ++s) {
        const int c = image_.edgeNodeCounts[s];
        int start, num;
        if (s <= 1) {
            if (c <= 2) {
                idx += 2;
                continue;
            }
            start = idx + 1;
            num = c - 2;
            idx += c;
        } else {
            if (c < 1) continue;
            start = idx;
            num = c;
            idx += c;
        }
        for (int j = start; j < start + num; ++j) {
            if (nodes[j].type != 3) {
                nodes[j].flags |= 0x40;
                ++removed[s];
            }
        }
    }
    for (int j = 0; j < 2 * (w_ + h_); ++j) nodes[j].flags |= 0x80;

    if (!antiAliased) {
        // Staircase removal on contours of >= 20 nodes. Walking the contour (with
        // 30 nodes of wrap-around), each run of steps without a direction reversal is recorded as two
        // bit strings - x-steps and y-steps, newest bit lowest. Once a run has 5..31 steps it is
        // matched against the staircase patterns; on a match the replacement mask decides which of
        // the last len+1 nodes survive (0 bits are removed). Nodes passed after the first 30 steps
        // are marked 0x80, which ends the wrap-around at the first already-visited node.
        for (Region& r : regions_) {
            const int n = static_cast<int>(r.nodes.size());
            if (n <= 19) continue;
            const int limit = n + 30;
            int bitsX = 0, bitsY = 0, count = 0, prevDx = 0, prevDy = 0, j = 0;
            auto reset = [&] { count = bitsX = bitsY = prevDx = prevDy = 0; };
            for (int i = 1; i <= limit; ++i) {
                Node& a = nodes[r.nodes[j]];
                j = i % n;
                Node& b = nodes[r.nodes[j]];
                if ((a.flags & 0x40) || (b.flags & 0x80) || a.type == 3 || b.type == 3) {
                    reset();
                    continue;
                }
                if (i > 30) a.flags |= 0x80;
                const int dx = static_cast<int>(std::round(a.x - b.x));
                const int dy = static_cast<int>(std::round(a.y - b.y));
                if (dx != 0 && dx != prevDx && prevDx != 0) {  // x reversal: new run
                    prevDx = dx;
                    count = 1;
                    bitsY = dy != 0;
                    bitsX = 1;
                    if (dy != 0) prevDy = dy;
                    continue;
                }
                if (prevDy != 0 && dy != 0 && dy != prevDy) {  // y reversal: new run
                    prevDy = dy;
                    count = 1;
                    bitsY = 1;
                    bitsX = dx != 0;
                    if (dx != 0) prevDx = dx;
                    continue;
                }
                bitsX <<= 1;
                bitsY <<= 1;
                ++count;
                if (dx != 0) {
                    bitsX |= 1;
                    prevDx = dx;
                }
                if (dy != 0) {
                    bitsY |= 1;
                    prevDy = dy;
                }
                if (count == 32) count = 31;
                if (count <= 4) continue;
                int keep = 0;
                const int len = patternMatch(bitsX, bitsY, count, &keep);
                if (len <= 0) continue;
                for (int c = 0; c <= len; ++c, keep >>= 1)
                    if (!(keep & 1)) nodes[r.nodes[(i + n - c) % n]].flags |= 0xc0;
                reset();
            }
        }
    }

    for (int s = 0; s < 4; ++s) image_.edgeNodeCounts[s] -= removed[s];

    const int total = static_cast<int>(nodes.size());
    std::vector<int> remap(total, -1);
    int kept = 0;
    for (int i = 0; i < total; ++i) {
        if (nodes[i].flags & 0x40) continue;
        remap[i] = kept;
        if (kept != i) nodes[kept] = std::move(nodes[i]);
        ++kept;
    }
    nodes.resize(kept);

    for (Region& r : regions_) {
        const int n = static_cast<int>(r.nodes.size());
        int k = 0, last = 0;
        for (int i = 0; i < n; ++i) {
            const int m = remap[r.nodes[i]];
            if (m < 0) continue;
            r.nodes[k] = m;
            r.elems[k].flag = r.elems[i].flag;
            r.elems[k].neighbor = r.elems[i].neighbor;
            r.elems[k].steps = i - last;
            last = i;
            ++k;
        }
        r.elems[0].steps = n - last;
        r.nodes.resize(k);
        r.elems.resize(k);
    }
}

void BoundaryTracer::execute() {
    w_ = map_.width();
    h_ = map_.height();
    nodeIndex_.assign(static_cast<size_t>(w_ + 1) * (h_ + 1), 0);
    visited_.assign(static_cast<size_t>(w_ + 1) * (h_ + 1), 0);

    const int numNodes = enumerateBoundaryNodes();
    image_.nodes.clear();
    image_.setNodesSize(numNodes);
    image_.initializeEdgeNodeCounts(w_, h_);
    for (int y = 0; y <= h_; ++y)
        for (int x = 0; x <= w_; ++x)
            if (nodeIndex(x, y) >= 0) image_.setInitialNodePositionAndType(nodeIndex(x, y), x, y);

    std::fill(visited_.begin(), visited_.end(), -1);
    for (Region& r : regions_) r.nodes.clear();

    // Trace each region once, from the first corner node of its first pixel in scan order.
    for (int y = 0; y < h_; ++y)
        for (int x = 0; x < w_; ++x) {
            const int seg = id(x, y);
            if (!regions_[seg].nodes.empty()) continue;
            for (int k = 0; k < 4; ++k) {
                const int nx = x + kCornerDx[k], ny = y + kCornerDy[k];
                if (nodeIndex(nx, ny) >= 0) {
                    traceBoundary(seg, nx, ny);
                    break;
                }
            }
        }

    // Junctions: corners where the four surrounding pixels (outside = -1) have at most one equal
    // pair become fixed.
    for (int y = 0; y <= h_; ++y)
        for (int x = 0; x <= w_; ++x) {
            const int ni = nodeIndex(x, y);
            if (ni < 0) continue;
            int p[4];
            for (int k = 0; k < 4; ++k) {
                const int px = x + kPixelDx[k], py = y + kPixelDy[k];
                p[k] = (px >= 0 && px < w_ && py >= 0 && py < h_) ? id(px, py) : -1;
            }
            const int equal = (p[0] == p[1]) + (p[0] == p[2]) + (p[0] == p[3]) + (p[1] == p[2]) +
                              (p[1] == p[3]) + (p[2] == p[3]);
            if (equal <= 1) image_.nodes[ni].type = 3;
        }

    removeRedundantNodes(antiAliased_);

    // Parent: the first neighbour whose own (outer) contour does not touch this region.
    const int nr = static_cast<int>(regions_.size());
    for (int i = 0; i < nr; ++i) {
        Region& r = regions_[i];
        r.parent = -1;
        for (int nb : r.neighbors) {
            const std::vector<int>& back = regions_[nb].neighbors;
            if (std::find(back.begin(), back.end(), i) == back.end()) {
                r.parent = nb;
                break;
            }
        }
    }

    nodeIndex_.clear();
    visited_.clear();
}

}  // namespace isvg
