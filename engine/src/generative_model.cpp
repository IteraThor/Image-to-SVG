#include "isvg/generative_model.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace isvg {

namespace {
constexpr double kPi = 3.141592653589793;
constexpr double kTwoPi = 6.283185307179586;
}  // namespace

void GenerativeModel::regionColorF(int r, float out[4]) const {
    const uint8_t* c = vi_.regions[r].color;
    for (int k = 0; k < 4; ++k) out[k] = static_cast<float>(c[k]) * kInv255f;
}

void GenerativeModel::observedF(int px, int py, float out[4]) const {
    const uint8_t* p = img_.at(py * img_.w + px);
    for (int k = 0; k < 4; ++k) out[k] = static_cast<float>(p[k]) * kInv255f;
}

void GenerativeModel::setup() {
    w_ = map_.width();
    h_ = map_.height();
    cellHead_.assign(static_cast<size_t>(w_) * h_, -1);
    cellTail_.assign(static_cast<size_t>(w_) * h_, -1);
    mask_.assign(static_cast<size_t>(w_) * h_, 0);
    model_.assign(static_cast<size_t>(w_) * h_ * 4, 0.0f);
    crumbs_.clear();
    const std::vector<int32_t>& ids = map_.rawIds();
    for (size_t i = 0; i < ids.size(); ++i) regionColorF(ids[i], &model_[i * 4]);
}

// --- crumbs ---------------------------------------------------------------------------------

void GenerativeModel::dropBreadCrumb(int ix, int iy, int region, int k, uint32_t flags, double x, double y) {
    if (ix < 0 || ix >= w_ || iy < 0 || iy >= h_) return;
    const int ci = newCrumb();
    BreadCrumb& c = crumbs_[ci];
    c.x = static_cast<float>(x);
    c.y = static_cast<float>(y);
    c.next = -1;
    c.region = region;
    c.flags = flags;
    c.k = k;
    float perim = -1.0f;
    if (flags & 1) {
        if (!(flags & 0x100)) {
            const double d = x - static_cast<double>(ix);
            perim = (y == static_cast<double>(iy)) ? static_cast<float>(d) : static_cast<float>(3.0 - d);
        } else {
            const double d = y - static_cast<double>(iy);
            perim = (x == static_cast<double>(ix)) ? static_cast<float>(4.0 - d) : static_cast<float>(d + 1.0);
        }
    }
    c.perim = perim;
    c.winding = static_cast<float>(winding_);
    int& head = cellHead(ix, iy);
    int& tail = cellTail_[iy * w_ + ix];
    if (head < 0) head = ci;
    else crumbs_[tail].next = ci;
    tail = ci;
}

void GenerativeModel::dropBreadCrumbs(int region, int k, int kNext) {
    const Region& reg = vi_.regions[region];
    const Node& a = vi_.nodes[reg.nodes[k]];
    const Node& b = vi_.nodes[reg.nodes[kNext]];
    const double x0 = a.x, y0 = a.y, x1 = b.x, y1 = b.y;
    const double dx = x1 - x0, dy = y1 - y0;

    // Start point
    {
        const double fx = std::floor(x0), fy = std::floor(y0);
        int ix = static_cast<int>(fx), iy = static_cast<int>(fy);
        uint32_t flags;
        if (x0 != fx && y0 != fy) {
            flags = 4;
        } else {
            const int lx = 0.0 > dx, ly = 0.0 > dy;
            if (x0 == fx && y0 == fy) {
                ix -= lx;
                iy -= ly;
                flags = 0x187;
            } else if (x0 == fx) {
                ix -= lx;
                flags = 0x105;
            } else {
                iy -= ly;
                flags = 0x85;
            }
        }
        dropBreadCrumb(ix, iy, region, k, flags, x0, y0);
    }

    // Grid-line crossings, in order along the edge (t from the start node).
    double cx = x0, cy = y0;
    for (;;) {
        double nx = 0, ny = 0, tX, tY;
        if (dx < 0.0) {
            nx = std::ceil(cx) - 1.0;
            tX = (nx - x0) / dx;
        } else if (dx > 0.0) {
            nx = std::floor(cx) + 1.0;
            tX = (nx - x0) / dx;
        } else {
            tX = 2.0;
        }
        bool vertical;  // crossing a vertical grid line (x = nx) next
        if (dy != 0.0 || dx == 0.0) {
            if (dy < 0.0) {
                ny = std::ceil(cy) - 1.0;
                tY = (ny - y0) / dy;
            } else if (dy > 0.0) {
                ny = std::floor(cy) + 1.0;
                tY = (ny - y0) / dy;
            } else {
                break;  // dx == dy == 0
            }
            if (tX < 1.0 || std::isnan(tX)) {
                vertical = tY > tX;
            } else {
                if (!(tY < 1.0)) break;
                vertical = false;
            }
        } else {  // dy == 0, dx != 0
            if (!(tX < 1.0) && !std::isnan(tX)) break;
            tY = 2.0;
            vertical = tY > tX;
        }
        uint32_t flags;
        int dix, diy;
        if (vertical) {
            cx = nx;
            cy = tX * dy + y0;
            flags = 0x100;
            dix = -1;
            diy = 0;
        } else {
            cx = tY * dx + x0;
            cy = ny;
            flags = 0x80;
            dix = 0;
            diy = -1;
        }
        const int ix = static_cast<int>(cx), iy = static_cast<int>(cy);
        if (static_cast<double>(ix) == cx && static_cast<double>(iy) == cy) {
            dropBreadCrumb(ix, iy, region, k, 0x183, cx, cy);
            dropBreadCrumb(ix - 1, iy, region, k, 0x183, cx, cy);
            dropBreadCrumb(ix, iy - 1, region, k, 0x183, cx, cy);
            dropBreadCrumb(ix - 1, iy - 1, region, k, 0x183, cx, cy);
        } else {
            dropBreadCrumb(ix, iy, region, k, flags | 1, cx, cy);
            dropBreadCrumb(ix + dix, iy + diy, region, k, flags | 1, cx, cy);
        }
    }

    // End point (only when it lies on a grid line)
    const double fx1 = std::floor(x1), fy1 = std::floor(y1);
    if (x1 != fx1 && y1 != fy1) return;
    int ix = static_cast<int>(x1), iy = static_cast<int>(y1);
    const int gx = dx > 0.0, gy = dy > 0.0;
    uint32_t flags;
    if (x1 == fx1 && y1 == fy1) {
        ix -= gx;
        iy -= gy;
        flags = 0x197;
    } else if (x1 == fx1) {
        ix -= gx;
        flags = 0x115;
    } else {
        iy -= gy;
        flags = 0x95;
    }
    dropBreadCrumb(ix, iy, region, k, flags, x1, y1);
}

// --- closing the per-pixel polygons -------------------------------------------------------------

int GenerativeModel::insertCornerNode(int c, int n, int px, int py, int corner) {
    const int k = (corner + 4) % 4;
    const int ci = newCrumb();
    BreadCrumb& nc = crumbs_[ci];
    nc.next = n;
    nc.winding = 0.0f;
    nc.region = crumbs_[c].region;
    nc.perim = static_cast<float>(static_cast<double>(k));
    nc.flags = 0x18b;
    nc.k = -2;
    nc.x = static_cast<float>(px + (static_cast<unsigned>(k - 1) <= 1u ? 1 : 0));
    nc.y = static_cast<float>(py + (static_cast<unsigned>(k - 2) <= 1u ? 1 : 0));
    crumbs_[c].next = ci;
    return ci;
}

int GenerativeModel::insertCornerNodes(int c, int n, int px, int py) {
    // Walk the cell boundary from crumb c to crumb n, inserting the cell corners passed, in the
    // direction given by the sign of the contour's turning between them.
    const BreadCrumb& C = crumbs_[c];
    const BreadCrumb& N = crumbs_[n];
    double d = static_cast<double>(N.winding - C.winding);
    if (C.k > N.k) d = d + kTwoPi;
    const float np = N.perim, cp = C.perim;
    if (d < 0.0) {
        int lo = static_cast<int>(std::floor(np));
        if (np > cp) lo -= 4;
        int i = static_cast<int>(static_cast<double>(cp) - 5.0) + 4;
        while (lo < i) {
            c = insertCornerNode(c, n, px, py, i);
            --i;
        }
    } else {
        int hi = static_cast<int>(std::ceil(np));
        if (cp > np) hi += 4;
        int i = static_cast<int>(static_cast<double>(cp) + 1.0);
        while (hi > i) {
            c = insertCornerNode(c, n, px, py, i);
            ++i;
        }
    }
    return c;
}

int GenerativeModel::checkInsertCornerNodes(int c, int px, int py) {
    const BreadCrumb& C = crumbs_[c];
    const int n = C.next;
    const BreadCrumb& N = crumbs_[n];
    const uint32_t both = C.flags & N.flags;
    if (!(both & 1)) return c;
    if (C.k == N.k) return c;
    const int len = static_cast<int>(vi_.regions[C.region].nodes.size());
    if (N.k == (C.k + 1) % len && C.perim == N.perim) return c;
    if (both & 2) {
        if (C.x == N.x) return c;
        if (C.y == N.y) return c;
    }
    return insertCornerNodes(c, n, px, py);
}

int GenerativeModel::checkCloseRegion(int last, int first, int px, int py) {
    const BreadCrumb& L = crumbs_[last];
    const BreadCrumb& F = crumbs_[first];
    if (!((L.flags & F.flags) & 1)) return last;
    const int len = static_cast<int>(vi_.regions[L.region].nodes.size());
    if ((L.k + 1) % len == F.k && L.perim == F.perim) return last;
    const int after = L.next;
    const int r = insertCornerNodes(last, first, px, py);
    crumbs_[r].next = after;
    return r;
}

void GenerativeModel::cleanOutTwoCrumbRegions(int px, int py) {
    // Drop runs of fewer than three crumbs. `prev` keeps its last value: after a
    // removal it points at the removed run, so a second consecutive removal is not unlinked.
    int& head = cellHead(px, py);
    int prev = -1;
    int r = head;
    while (r >= 0) {
        int next = crumbs_[r].next;
        int last;
        if (next < 0) {
            last = r;
        } else {
            const int region = crumbs_[r].region;
            last = r;
            int extra = 0;
            int a = next;
            for (;;) {
                if (crumbs_[a].region != region) {
                    if (extra > 1) {
                        prev = last;
                        r = crumbs_[prev].next;
                        goto next_run;
                    }
                    next = a;
                    break;
                }
                const int nn = crumbs_[a].next;
                ++extra;
                last = a;
                if (nn < 0) {
                    if (extra != 1) return;
                    next = -1;
                    break;
                }
                a = nn;
            }
        }
        if (head == r) head = next;
        else crumbs_[prev].next = next;
        prev = last;
        r = crumbs_[prev].next;
    next_run:;
    }
}

void GenerativeModel::closePixel(int px, int py) {
    int& head = cellHead(px, py);
    int d = head;
    while (d >= 0) {
        int b = d;
        int a = crumbs_[b].next;
        bool found = false;
        while (a >= 0) {
            if (crumbs_[b].region == crumbs_[a].region) {
                found = true;
                break;
            }
            if (d == head) {
                head = a;
                d = a;
            }
            b = a;
            a = crumbs_[b].next;
        }
        if (!found) {
            if (d == head) head = -1;
            break;
        }
        // b starts a run of >= 2 crumbs of one region; close it into a polygon.
        int c = b;
        do {
            c = checkInsertCornerNodes(c, px, py);
            c = crumbs_[c].next;
        } while (crumbs_[c].next >= 0 && crumbs_[crumbs_[c].next].region == crumbs_[c].region);
        const int r = checkCloseRegion(c, b, px, py);
        d = crumbs_[r].next;
    }
    cleanOutTwoCrumbRegions(px, py);
}

// --- colour and potential -----------------------------------------------------------------------

int GenerativeModel::regionToIsland(int region) const {
    for (int i = 0; i < slotCount_; ++i)
        if (slotRegion_[i] == region) return i;
    return -1;
}

void GenerativeModel::resetPixelToObserved(int px, int py) {
    // Unresolvable coverage: drop the pixel's crumbs and make it match the input exactly.
    cellHead(px, py) = -1;
    float o[4];
    observedF(px, py, o);
    float* m = modelAt(px, py);
    for (int k = 0; k < 4; ++k) m[k] = 0.0f + o[k];
}

double GenerativeModel::computePixelColor(int px, int py) {
    islandParent_ = -1;
    int c = cellHead(px, py);
    if (c < 0) return 1.0;
    slotCount_ = 0;
    float* m = modelAt(px, py);
    m[0] = m[1] = m[2] = m[3] = 0.0f;
    double total = 0.0;
    bool bad = false;
    while (c >= 0) {
        // Shoelace area of this region's polygon in the cell.
        const BreadCrumb& f = crumbs_[c];
        const int region = f.region;
        const float fx = f.x, fy = f.y;
        float lx = fx, ly = fy;
        double sum = 0.0;
        int last = c;
        for (int a = f.next; a >= 0 && crumbs_[a].region == region; a = crumbs_[a].next) {
            const float cxv = crumbs_[a].x, cyv = crumbs_[a].y;
            sum = sum + static_cast<double>(lx * cyv - ly * cxv);
            lx = cxv;
            ly = cyv;
            last = a;
        }
        const double area = (static_cast<double>(lx * fy - ly * fx) + sum) * 0.5;
        float rc[4];
        regionColorF(region, rc);
        const float fa = static_cast<float>(area);
        for (int k = 0; k < 4; ++k) m[k] = m[k] + fa * rc[k];
        total = total + area;
        if (slotCount_ > 9) std::abort();
        slotRegion_[slotCount_] = region;
        slotArea_[slotCount_] = area;
        if (0.0 > area) {
            mask_[py * w_ + px] |= 1;
            if (static_cast<int>(vi_.regions[crumbs_[last].region].nodes.size()) >= 5) bad = true;
        }
        ++slotCount_;
        c = crumbs_[last].next;
    }

    if (!(std::fabs(total - 1.0) > 0.001) && !bad) return total;
    if (!bad) {
        // Coverage does not add up: a region in the cell may be an island whose enclosing region
        // left no crumbs here. Credit the missing area to that parent.
        for (int i = 0; i < slotCount_; ++i) {
            const int parent = vi_.regions[slotRegion_[i]].parent;
            if (parent >= 0 && regionToIsland(parent) < 0) {
                if (islandParent_ < 0) islandParent_ = parent;
            }
        }
        const int n0 = slotCount_;
        if (islandParent_ >= 0) {
            float pc[4];
            regionColorF(islandParent_, pc);
            for (int k = 0; k < 4; ++k) m[k] = m[k] + pc[k];
            total = total + 1.0;
            if (slotCount_ > 9) std::abort();
            slotRegion_[slotCount_] = islandParent_;
            slotArea_[slotCount_] = 1.0;
            ++slotCount_;
        }
        for (int i = 0; i < n0; ++i) {
            const int parent = vi_.regions[slotRegion_[i]].parent;
            const int j = regionToIsland(parent);
            if (j >= 0) {
                const double ar = slotArea_[i];
                float pc[4];
                regionColorF(parent, pc);
                const float fa = static_cast<float>(ar);
                for (int k = 0; k < 4; ++k) m[k] = m[k] - fa * pc[k];
                slotArea_[j] = slotArea_[j] - ar;
                total = total - slotArea_[i];
            } else if (parent < 0 && islandParent_ >= 0) {  // top-level region beside the island
                const double ar = slotArea_[i];
                float pc[4];
                regionColorF(islandParent_, pc);
                const float fa = static_cast<float>(ar);
                for (int k = 0; k < 4; ++k) m[k] = m[k] - fa * pc[k];
                if (slotCount_ > 9) std::abort();
                slotArea_[slotCount_] = slotArea_[slotCount_] - ar;  // original indexes one past the parent
                total = total - slotArea_[i];
            }
        }
        bool negative = false;
        for (int i = 0; i < slotCount_; ++i)
            if (0.0 > slotArea_[i]) {
                negative = true;
                break;
            }
        if (!negative) {
            if (!(std::fabs(total - 1.0) > 0.001)) return total;
            mask_[py * w_ + px] |= 2;
            for (int i = 0; i < slotCount_; ++i) flagged_.insert(slotRegion_[i]);
            resetPixelToObserved(px, py);
            return -1.0;
        }
    }
    for (int i = 0; i < slotCount_; ++i) flagged_.insert(slotRegion_[i]);
    resetPixelToObserved(px, py);
    return -1.0;
}

double GenerativeModel::findPotential() {
    crumbs_.clear();
    flagged_.clear();
    std::fill(cellHead_.begin(), cellHead_.end(), -1);
    std::fill(cellTail_.begin(), cellTail_.end(), -1);

    const int nr = static_cast<int>(vi_.regions.size());
    for (int r = 0; r < nr; ++r) {
        const Region& reg = vi_.regions[r];
        const int n = static_cast<int>(reg.nodes.size());
        if (n < 2) std::abort();
        const Node& p0 = vi_.nodes[reg.nodes[0]];
        const Node& p1 = vi_.nodes[reg.nodes[1]];
        double prev = std::atan2(p1.y - p0.y, p1.x - p0.x);
        winding_ = 0.0;
        dropBreadCrumbs(r, 0, 1 % n);
        for (int i = 1; i <= n; ++i) {
            const Node& a = vi_.nodes[reg.nodes[i % n]];
            const Node& b = vi_.nodes[reg.nodes[(i + 1) % n]];
            const double ang = std::atan2(b.y - a.y, b.x - a.x);
            double d = ang - prev;
            if (d > kPi) winding_ = (d - kTwoPi) + winding_;
            else if (-kPi > d) winding_ = (d + kTwoPi) + winding_;
            else winding_ = d + winding_;
            if (i < n) dropBreadCrumbs(r, i, (i + 1) % n);
            prev = ang;
        }
        if (std::fabs(winding_ - kTwoPi) > 0.001) flagged_.insert(r);
    }

    double e = 0.0;
    for (int y = 0; y < h_; ++y)
        for (int x = 0; x < w_; ++x) {
            if (cellHead(x, y) >= 0) {
                closePixel(x, y);
                computePixelColor(x, y);
            }
            float o[4];
            observedF(x, y, o);
            const float* m = modelAt(x, y);
            float s = 0.0f;
            const float d0 = m[0] - o[0], d1 = m[1] - o[1], d2 = m[2] - o[2], d3 = m[3] - o[3];
            s = d0 * d0 + s;
            s = d1 * d1 + s;
            s = s + d2 * d2;
            s = d3 * d3 + s;
            e = e + static_cast<double>(s);
        }
    return 0.5 * e;
}

// --- gradient -----------------------------------------------------------------------------------

namespace {
// GenerativeModel::computePartialGradients: a crumb on edge a→b lying on a grid line
// has its free coordinate i interpolated between the nodes; d(pos_i)/d(a), d(pos_i)/d(b).
void partialGradients(const float pos[2], const Node& a, const Node& b, int i, double outA[2], double outB[2]) {
    const int j = (i == 0) ? 1 : 0;
    const double A[2] = {a.x, a.y}, B[2] = {b.x, b.y};
    const double bj = B[j];
    const double denom = A[j] - bj;
    const double t = (static_cast<double>(pos[j]) - bj) / denom;
    outA[i] = t;
    const double m = -((A[i] - B[i]) / denom);
    outA[j] = m * t;
    outB[i] = 1.0 - t;
    outB[j] = m * (1.0 - t);
}
}  // namespace

void GenerativeModel::findPixelGradient(std::vector<double>& grad, int px, int py) {
    computePixelColor(px, py);
    int c = cellHead(px, py);
    float o[4];
    observedF(px, py, o);
    const float* m = modelAt(px, py);
    const float res[4] = {m[0] - o[0], m[1] - o[1], m[2] - o[2], m[3] - o[3]};
    if (c < 0) return;

    int prevRegion = -1, runStart = c, len = 0;
    double coef = 0.0;
    const Region* reg = nullptr;
    auto addNode = [&](int node, double gx, double gy) {
        grad[2 * node] = grad[2 * node] + gx;
        grad[2 * node + 1] = grad[2 * node + 1] + gy;
    };
    while (c >= 0) {
        const BreadCrumb& C = crumbs_[c];
        if (C.region != prevRegion) {
            reg = &vi_.regions[C.region];
            float col[4];
            regionColorF(C.region, col);
            const int parent = reg->parent;
            if ((parent | islandParent_) >= 0) {
                float pc[4];
                regionColorF(parent, pc);
                for (int k = 0; k < 4; ++k) col[k] = col[k] - pc[k];
            }
            len = static_cast<int>(reg->nodes.size());
            float s = col[0] * res[0] + 0.0f;
            s = s + col[1] * res[1];
            s = s + col[2] * res[2];
            s = s + col[3] * res[3];
            coef = static_cast<double>(s) * 0.5;
            prevRegion = C.region;
            runStart = c;
        }
        const int nx = C.next;
        const int n2i = (nx >= 0 && crumbs_[nx].region == C.region) ? nx : runStart;
        const BreadCrumb& N = crumbs_[n2i];

        // Area term 0.5·(C.x·N.y − C.y·N.x) w.r.t. C's position.
        if (C.k >= 0) {
            int k = C.k;
            if (C.flags & 0x10) k = (k + 1) % len;
            if (C.flags & 1) {
                const int ia = reg->nodes[k], ib = reg->nodes[(k + 1) % len];
                const Node& a = vi_.nodes[ia];
                const Node& b = vi_.nodes[ib];
                int i = 1;
                float f;
                if ((!(C.flags & 0x80) && a.x != b.x) || a.y == b.y) {
                    f = -N.x;
                } else {
                    i = 0;
                    f = N.y;
                }
                double ga[2], gb[2];
                const float pos[2] = {C.x, C.y};
                partialGradients(pos, a, b, i, ga, gb);
                const double s = static_cast<double>(f) * coef;
                addNode(ia, ga[0] * s, ga[1] * s);
                addNode(ib, s * gb[0], s * gb[1]);
            } else {
                addNode(reg->nodes[k], static_cast<double>(N.y) * coef, static_cast<double>(-N.x) * coef);
            }
        }
        // ... and w.r.t. N's position.
        if (N.k >= 0) {
            int k = N.k;
            if (N.flags & 0x10) k = (k + 1) % len;
            if (!(N.flags & 1)) {
                addNode(reg->nodes[k], static_cast<double>(-C.y) * coef, static_cast<double>(C.x) * coef);
            } else {
                const int ia = reg->nodes[k], ib = reg->nodes[(k + 1) % len];
                const Node& a = vi_.nodes[ia];
                const Node& b = vi_.nodes[ib];
                int i = 1;
                float f;
                if ((!(N.flags & 0x80) && a.x != b.x) || a.y == b.y) {
                    f = C.x;
                } else {
                    i = 0;
                    f = -C.y;
                }
                double ga[2], gb[2];
                const float pos[2] = {N.x, N.y};
                partialGradients(pos, a, b, i, ga, gb);
                const double s = static_cast<double>(f) * coef;
                addNode(ia, ga[0] * s, ga[1] * s);
                addNode(ib, s * gb[0], s * gb[1]);
            }
        }
        c = nx;
    }
}

void GenerativeModel::findGradient(std::vector<double>& grad) {
    for (int y = 0; y < h_; ++y)
        for (int x = 0; x < w_; ++x)
            if (cellHead(x, y) >= 0) findPixelGradient(grad, x, y);
}

bool GenerativeModel::checkFlaggedRegionsForInversions() {
    int n = 0;
    for (int r : flagged_) n += vi_.flagInversionNodesForRegion(r);
    flagged_.clear();
    return n > 0;
}

}  // namespace isvg
