#include "isvg/contour_smoother.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <utility>
#include <vector>

#include "isvg/prng.h"

namespace isvg {

ContourSmoother::ContourSmoother(VectorImage& vi, PixelMap& map, const Image& img, const Palette& palette,
                                 bool antiAliased, const ContourSmootherParams& p)
    : vi_(vi), map_(map), img_(img), palette_(palette), antiAliased_(antiAliased), p_(p), gm_(vi, map, img) {}

// --- priors -------------------------------------------------------------------------------------

double ContourSmoother::quadraticPrior(int a, int b, int c) const {
    const Node& A = vi_.nodes[a];
    const Node& B = vi_.nodes[b];
    const Node& C = vi_.nodes[c];
    const double sn = stepsNext_, st = stepsThis_;
    const double dx = (C.x - B.x) / sn - (B.x - A.x) / st;
    const double dy = (C.y - B.y) / sn - (B.y - A.y) / st;
    return dy * dy + dx * dx + 0.0;
}

void ContourSmoother::quadraticPriorGradient(int a, int b, int c, double o1[2], double o2[2]) const {
    const Node& A = vi_.nodes[a];
    const Node& B = vi_.nodes[b];
    const Node& C = vi_.nodes[c];
    const double sn = stepsNext_, st = stepsThis_;
    const double dx = (C.x - B.x) / sn - (B.x - A.x) / st;
    const double dy = (C.y - B.y) / sn - (B.y - A.y) / st;
    const double tx = dx + dx, ty = dy + dy;
    o1[0] = tx / st;
    o1[1] = ty / st;
    o2[0] = tx / sn;
    o2[1] = ty / sn;
}

double ContourSmoother::angularPrior(int a, int b, int c) const {
    const Node& A = vi_.nodes[a];
    const Node& B = vi_.nodes[b];
    const Node& C = vi_.nodes[c];
    const double ux = B.x - A.x, uy = B.y - A.y;
    const double vx = C.x - B.x, vy = C.y - B.y;
    const double lv = std::sqrt(vy * vy + vx * vx + 0.0);
    const double lu = std::sqrt(uy * uy + ux * ux + 0.0);
    const double dot = (vy / lv) * (uy / lu) + (vx / lv) * (ux / lu) + 0.0;
    const double s = (2.0 - (dot + dot)) + 0.001;
    const double ang = std::sqrt(s);
    const double len = lv / static_cast<double>(stepsNext_) - lu / static_cast<double>(stepsThis_);
    return ang + 0.0 + len * len * p_.lengthPenaltyWeights[phase_];
}

void ContourSmoother::angularPriorGradient(int a, int b, int c, double o1[2], double o2[2]) const {
    const Node& A = vi_.nodes[a];
    const Node& B = vi_.nodes[b];
    const Node& C = vi_.nodes[c];
    const double ux = B.x - A.x, uy = B.y - A.y;
    const double vx = C.x - B.x, vy = C.y - B.y;
    const double lu = std::sqrt(uy * uy + ux * ux + 0.0);
    const double hux = ux / lu, huy = uy / lu;
    const double lv = std::sqrt(vy * vy + vx * vx + 0.0);
    const double hvx = vx / lv, hvy = vy / lv;
    double dot = huy * hvy + hux * hvx + 0.0;
    const double r = std::sqrt((2.0 - (dot + dot)) + 0.001);
    dot = -dot;
    const double iu = 1.0 / lu, iv = -1.0 / lv;
    o1[0] = ((hux * dot + hvx) * iu) / r;
    o1[1] = ((huy * dot + hvy) * iu) / r;
    o2[0] = ((dot * hvx + hux) * iv) / r;
    o2[1] = ((dot * hvy + huy) * iv) / r;
    const double sn = stepsNext_, st = stepsThis_;
    double len = lv / sn - lu / st;
    len = (len + len) * p_.lengthPenaltyWeights[phase_];
    o1[0] = (hux * len) / st + o1[0];
    o1[1] = (huy * len) / st + o1[1];
    o2[0] = (hvx * len) / sn + o2[0];
    o2[1] = (hvy * len) / sn + o2[1];
}

double ContourSmoother::prior(int a, int b, int c) {
    switch (p_.priorTypes[phase_]) {
        case 0: return quadraticPrior(a, b, c);
        case 3: return angularPrior(a, b, c);
        default: std::abort();  // sqrt prior (type 1) is not used by any parameter file
    }
}

void ContourSmoother::priorGradient(int a, int b, int c, double o1[2], double o2[2]) {
    switch (p_.priorTypes[phase_]) {
        case 0: quadraticPriorGradient(a, b, c, o1, o2); return;
        case 3: angularPriorGradient(a, b, c, o1, o2); return;
        default: std::abort();
    }
}

// --- potential and gradient ---------------------------------------------------------------------

double ContourSmoother::findPotential() {
    const double ps = p_.priorStrengths[phase_];
    const double flaggedScale = p_.antiInvPotPriorScale / ps;
    double e = 0.0;
    const int nr = static_cast<int>(vi_.regions.size());
    for (int r = 0; r < nr; ++r) {
        const Region& reg = vi_.regions[r];
        const int n = static_cast<int>(reg.nodes.size());
        for (int k = 0; k < n; ++k) {
            const int a = reg.nodes[(n - 1 + k) % n], b = reg.nodes[k], c = reg.nodes[(k + 1) % n];
            stepsThis_ = reg.elems[k].steps;
            stepsNext_ = reg.elems[(k + 1) % n].steps;
            if (vi_.nodes[b].flags & 0x10) e = e + quadraticPrior(a, b, c) * flaggedScale;
            if (reg.n < 8) {
                e = e + quadraticPrior(a, b, c);
            } else if (!(reg.elems[k].flag & 1)) {
                e = e + prior(a, b, c);
            }
        }
    }
    e = e * ps;
    const double air = ps * p_.airPressureWeight;
    for (int r = 0; r < nr; ++r) {
        if (vi_.regions[r].n >= 8) continue;
        vi_.updateArea(r);
        const double d = vi_.regions[r].area - static_cast<double>(vi_.regions[r].n);
        e = e + d * d * air;
    }
    switch (p_.measurementTypes[phase_]) {
        case 1: e = e + gm_.findPotential(); break;
        case 0:  // squared displacement from the traced position
            for (int i = 0; i < static_cast<int>(vi_.nodes.size()); ++i) {
                double mv[2];
                vi_.movedVector(i, mv);
                e = e + (mv[1] * mv[1] + (mv[0] * mv[0] + 0.0));
            }
            break;
        case 2:  // squared displacement beyond 0.7
            for (int i = 0; i < static_cast<int>(vi_.nodes.size()); ++i) {
                double mv[2];
                vi_.movedVector(i, mv);
                double d = std::sqrt(mv[1] * mv[1] + (mv[0] * mv[0] + 0.0)) - 0.7;
                d = 0.0 < d ? d * d : 0.0;
                e = e + d;
            }
            break;
        default: std::abort();
    }
    for (int i = 0; i < static_cast<int>(vi_.nodes.size()); ++i) {
        const uint8_t ai = vi_.nodes[i].antiInv;
        if (ai == 0) continue;
        const double coef = static_cast<double>(ai) * p_.antiInvPotMeasScale;
        double mv[2];
        vi_.movedVector(i, mv);
        e = e + (mv[1] * mv[1] + (mv[0] * mv[0] + 0.0)) * coef;
    }
    return e;
}

void ContourSmoother::findGradient() {
    std::fill(grad_.begin(), grad_.end(), 0.0);
    const double ps = p_.priorStrengths[phase_];
    const double flaggedScale = p_.antiInvPotPriorScale / ps;
    auto apply = [&](int a, int b, int c, const double o1[2], const double o2[2]) {
        grad_[2 * a] = grad_[2 * a] + o1[0];
        grad_[2 * a + 1] = grad_[2 * a + 1] + o1[1];
        grad_[2 * b] = grad_[2 * b] - o1[0];
        grad_[2 * b + 1] = grad_[2 * b + 1] - o1[1];
        grad_[2 * b] = grad_[2 * b] - o2[0];
        grad_[2 * b + 1] = grad_[2 * b + 1] - o2[1];
        grad_[2 * c] = grad_[2 * c] + o2[0];
        grad_[2 * c + 1] = grad_[2 * c + 1] + o2[1];
    };
    const int nr = static_cast<int>(vi_.regions.size());
    for (int r = 0; r < nr; ++r) {
        const Region& reg = vi_.regions[r];
        const int n = static_cast<int>(reg.nodes.size());
        for (int k = 0; k < n; ++k) {
            const int a = reg.nodes[(n - 1 + k) % n], b = reg.nodes[k], c = reg.nodes[(k + 1) % n];
            stepsThis_ = reg.elems[k].steps;
            stepsNext_ = reg.elems[(k + 1) % n].steps;
            double o1[2], o2[2];
            if (vi_.nodes[b].flags & 0x10) {
                quadraticPriorGradient(a, b, c, o1, o2);
                for (double* v : {o1, o2}) {
                    v[0] = v[0] * flaggedScale;
                    v[1] = v[1] * flaggedScale;
                }
                apply(a, b, c, o1, o2);
            }
            if (reg.n <= 7) {
                quadraticPriorGradient(a, b, c, o1, o2);
                apply(a, b, c, o1, o2);
            } else if (!(reg.elems[k].flag & 1)) {
                priorGradient(a, b, c, o1, o2);
                apply(a, b, c, o1, o2);
            }
        }
    }
    for (double& g : grad_) g = g * ps;

    const double air = p_.airPressureWeight * ps;
    for (int r = 0; r < nr; ++r) {
        const Region& reg = vi_.regions[r];
        if (reg.n > 7) continue;
        const double coef = (reg.area - static_cast<double>(reg.n)) * (air + air);
        const int n = static_cast<int>(reg.nodes.size());
        if (n <= 0) continue;
        for (int i = 0; i < n; ++i) {
            const int p = reg.nodes[i], q = reg.nodes[(i + 1) % n];
            const Node& Q = vi_.nodes[q];
            grad_[2 * p] = (Q.y * 0.5) * coef + grad_[2 * p];
            grad_[2 * p + 1] = (Q.x * -0.5) * coef + grad_[2 * p + 1];
            const Node& P = vi_.nodes[p];
            grad_[2 * q] = (P.y * -0.5) * coef + grad_[2 * q];
            grad_[2 * q + 1] = (P.x * 0.5) * coef + grad_[2 * q + 1];
        }
    }

    switch (p_.measurementTypes[phase_]) {
        case 1: gm_.findGradient(grad_); break;
        case 0:
            for (int i = 0; i < static_cast<int>(vi_.nodes.size()); ++i) {
                double mv[2];
                vi_.movedVector(i, mv);
                grad_[2 * i] = grad_[2 * i] + (mv[0] + mv[0]);
                grad_[2 * i + 1] = grad_[2 * i + 1] + (mv[1] + mv[1]);
            }
            break;
        default: std::abort();  // type 2's gradient is not used by any parameter file
    }

    const double meas2 = p_.antiInvPotMeasScale + p_.antiInvPotMeasScale;
    for (int i = 0; i < static_cast<int>(vi_.nodes.size()); ++i) {
        const uint8_t ai = vi_.nodes[i].antiInv;
        if (ai == 0) continue;
        double mv[2];
        vi_.movedVector(i, mv);
        const double t = static_cast<double>(ai) * meas2;
        grad_[2 * i] = grad_[2 * i] + t * mv[0];
        grad_[2 * i + 1] = grad_[2 * i + 1] + t * mv[1];
    }

    // Border nodes only slide along their image side.
    const int* ec = vi_.edgeNodeCounts;
    grad_[1] = 0.0;
    int last = -1, idx = 0;
    for (int i = 0; i < ec[0]; ++i) {
        grad_[2 * idx] = 0.0;
        last = idx++;
    }
    grad_[2 * last + 1] = 0.0;
    grad_[2 * idx + 1] = 0.0;
    int j = idx;
    for (int i = 0; i < ec[1]; ++i) grad_[2 * (j++)] = 0.0;
    grad_[2 * (j - 1) + 1] = 0.0;
    for (int i = 0; i < ec[2]; ++i) grad_[2 * (j++) + 1] = 0.0;
    for (int i = 0; i < ec[3]; ++i) grad_[2 * (j++) + 1] = 0.0;
}

std::vector<double*> ContourSmoother::statePointers() {
    std::vector<double*> v;
    for (Node& nd : vi_.nodes) {
        v.push_back(&nd.x);
        v.push_back(&nd.y);
    }
    return v;
}

std::vector<double*> ContourSmoother::gradientPointers() {
    std::vector<double*> v;
    for (double& g : grad_) v.push_back(&g);
    return v;
}

void ContourSmoother::afterMove() {
    vi_.clampEdgeNodes();
    vi_.clampMarginNodes();
}

void ContourSmoother::afterCGReset(bool keepBest) {
    if (antiAliased_ && !keepBest) opt_.resetBestYet();
}

void ContourSmoother::endOfIteration() {
    if (!antiAliased_) return;
    if (p_.optimizerType[phase_] == 0) {
        bool inv;
        if (opt_.iteration() % p_.inversionCheckPeriod == p_.inversionCheckPeriod - 1)
            inv = vi_.flagInversionNodes() != 0;
        else
            inv = gm_.checkFlaggedRegionsForInversions();
        if (inv) opt_.forceCGReset(false);
    } else {
        const int dec = vi_.decrementNodeCounts();
        const int chk = gm_.checkFlaggedRegionsForInversions() ? 1 : 0;
        if (chk + dec != 0) opt_.resetBestYet();
    }
}

// --- initialisation -------------------------------------------------------------------------------

void ContourSmoother::initializeNodePosition(int i) {
    // Nudge the node towards where the anti-aliased pixel colours say the boundary is: the
    // mixture weights of the four surrounding pixels (palette-index WeightsFinder) give an offset
    // along the two diagonals.
    const Node& nd = vi_.nodes[i];
    const int i0 = static_cast<int>(nd.fx) - 1, j0 = static_cast<int>(nd.fy) - 1;
    const std::vector<int32_t>& ids = map_.rawIds();
    const int w = map_.width(), h = map_.height();
    int P[4];
    P[0] = vi_.regions[ids[j0 * w + i0]].paletteIndex;
    P[1] = vi_.regions[ids[j0 * w + i0 + 1]].paletteIndex;
    P[2] = vi_.regions[ids[(j0 + 1) * w + i0]].paletteIndex;
    P[3] = vi_.regions[ids[(j0 + 1) * w + i0 + 1]].paletteIndex;
    double S[4] = {0, 0, 0, 0};
    for (int dy = 0; dy < 2; ++dy)
        for (int dx = 0; dx < 2; ++dx) {
            const int x = i0 + dx, y = j0 + dy;
            if (x < 0 || x >= w || y < 0 || y >= h) continue;
            wf_.computePixelF(x, y, -1, -1, nullptr, false);
            const int k = wf_.numNeighbors();
            if (k <= 0) continue;
            double& s = S[dy * 2 + dx];
            for (int q = 0; q < 4; ++q)
                for (int a = 0; a < k; ++a)
                    if (wf_.neighborIds()[a] == P[q]) {
                        if (a != 0) s = s + wf_.weight(a);
                        break;
                    }
        }
    const double d0 = S[0] - S[3], d1 = S[2] - S[1];
    double X = (M_[2] * d1 + M_[0] * d0) * 0.5;
    double Y = (d0 * M_[1] + d1 * M_[3]) * 0.5;
    X = -0.333 > X ? -0.333 : (0.333 < X ? 0.333 : X);
    Y = -0.333 > Y ? -0.333 : (0.333 < Y ? 0.333 : Y);
    Node& n = vi_.nodes[i];
    n.x = X + n.x;
    n.y = Y + n.y;
}

void ContourSmoother::initializeNodePositions() {
    if (!antiAliased_) return;
    wf_.clear();
    wf_.attach(map_, img_, [this](int idx) { return palette_.color(idx); });
    wf_.setIdMap([this](int id) { return vi_.regions[id].paletteIndex; });
    const double r = p_.perturbationRange;
    for (int i = vi_.numEdgeNodes(); i < static_cast<int>(vi_.nodes.size()); ++i) {
        initializeNodePosition(i);
        const double u1 = prngNext();
        const double u2 = prngNext();
        Node& n = vi_.nodes[i];
        n.x = n.x + (u1 * (r + r) - r);
        n.y = n.y + (u2 * (r + r) - r);
    }
}

void ContourSmoother::initialize(bool initPositions) {
    for (int m : p_.measurementTypes)
        if (m == 1) {
            gm_.setup();
            break;
        }
    grad_.assign(2 * vi_.nodes.size(), 0.0);
    opt_.setup(*this);
    vi_.initializeMarginNodes(map_.width(), map_.height());
    if (antiAliased_) vi_.buildIricNeighbors();
    if (initPositions) initializeNodePositions();
}

void ContourSmoother::executePhase(int phase) {
    phase_ = phase;
    if (!p_.isEnabled[phase]) return;
    if (p_.measurementTypes[phase] >= 3) std::abort();
    opt_.setSettings(p_.phase[phase]);
    opt_.doCG();
    if (p_.doPunctureCorners[phase]) punctureCorners();
}

void ContourSmoother::finish() {
    if (antiAliased_) vi_.flagInversionNodes();
}

void ContourSmoother::execute(bool initPositions) {
    initialize(initPositions);
    for (int p = 0; p < 3; ++p) executePhase(p);
    finish();
}

// --- straight edges -------------------------------------------------------------------------------

namespace {

constexpr double kLineMaxDev = 0.45;  // largest distance of a node from the line (pixels)
// A run grows until a node strays too far or the best parabola through it bows too much. A gentle
// arc always stops on the bow, so it ends up bowing close to kLineMaxBow; a noisy straight edge
// stops at a corner long before. Only runs bowing at most kLineAcceptBow are straightened.
constexpr double kLineMaxBow = 0.3;  // pixels
constexpr double kLineAcceptBow = 0.15;  // pixels
// Run ends that drift off the line into a smooth turn (the start of a rounded corner) are dropped.
constexpr double kLineEndDev = 0.1;  // pixels
constexpr double kLineSharpTurn = 0.7;  // sine of ~45 degrees: beyond, the end is a sharp corner
constexpr double kLineMinLength = 6.0;  // pixels
constexpr double kLineMinSmoothLength = 20.0;  // pixels, for a run without a sharp corner at either end
constexpr double kAxisSnapSin = 0.035;  // lines within ~2 degrees of horizontal/vertical become exact
constexpr int kLineMinNodes = 4;
constexpr double kCornerMinSin = 0.26;  // two runs meet in a corner above ~15 degrees
constexpr double kCornerMaxMove = 1.0;  // pixels
// Fairing of the curved stretches between straight runs and corners (fairCurves).
constexpr int kFairHalfWindow = 10;  // nodes; irons out pixel staircases up to ~20 px long
constexpr double kFairMaxMove = 0.5;  // pixels
constexpr double kFairCornerCos = 0.87;  // turns above ~30 degrees are corners

struct LineFit {
    double cx = 0, cy = 0, dx = 1, dy = 0;  // centroid and unit direction
    double maxDev = 0, bow = 0, length = 0;
};

// Total least-squares line through the nodes seq[i..j]. The bow tells a noisy straight edge from a
// gentle arc: a line fits both within a fraction of a pixel, but only the arc has a parabola
// h = a + b·s + c·s² (s along the line, h across) that bends by a real amount over the run.
LineFit fitLine(const std::vector<Node>& nodes, const std::vector<int>& seq, int i, int j) {
    LineFit f;
    const int m = j - i + 1;
    for (int k = i; k <= j; ++k) {
        f.cx += nodes[seq[k]].x;
        f.cy += nodes[seq[k]].y;
    }
    f.cx /= m;
    f.cy /= m;
    double sxx = 0, sxy = 0, syy = 0;
    for (int k = i; k <= j; ++k) {
        const double x = nodes[seq[k]].x - f.cx, y = nodes[seq[k]].y - f.cy;
        sxx += x * x;
        sxy += x * y;
        syy += y * y;
    }
    const double ang = 0.5 * std::atan2(2.0 * sxy, sxx - syy);
    f.dx = std::cos(ang);
    f.dy = std::sin(ang);
    double smin = 0, smax = 0, A[3][3] = {}, B[3] = {};
    for (int k = i; k <= j; ++k) {
        const double x = nodes[seq[k]].x - f.cx, y = nodes[seq[k]].y - f.cy;
        const double s = x * f.dx + y * f.dy, h = y * f.dx - x * f.dy;
        f.maxDev = std::max(f.maxDev, std::fabs(h));
        smin = std::min(smin, s);
        smax = std::max(smax, s);
        const double p[3] = {1.0, s, s * s};
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) A[r][c] += p[r] * p[c];
            B[r] += p[r] * h;
        }
    }
    f.length = smax - smin;
    // Quadratic coefficient by Cramer's rule.
    auto det3 = [](const double M[3][3]) {
        return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
               M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
    };
    const double d = det3(A);
    if (m >= 3 && std::fabs(d) > 1e-12) {
        double C[3][3];
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) C[r][c] = c == 2 ? B[r] : A[r][c];
        f.bow = std::fabs(det3(C) / d) * f.length * f.length * 0.25;
    }
    return f;
}

}  // namespace

void ContourSmoother::straightenEdges() {
    const int numEdge = vi_.numEdgeNodes();
    auto movable = [&](int node) { return node >= numEdge && !(vi_.nodes[node].flags & 0x80); };
    std::vector<char> onLine(vi_.nodes.size(), 0);
    for (const Region& reg : vi_.regions) {
        const int n = static_cast<int>(reg.nodes.size());
        if (n < kLineMinNodes) continue;
        // Walk the closed contour from a node a run cannot pass: an image-border node if there is one,
        // else the sharpest turn (most likely a corner). The start node is repeated at the end.
        int start = -1;
        double sharpest = -2.0;
        for (int k = 0; k < n && start < 0; ++k)
            if (!movable(reg.nodes[k])) start = k;
        if (start < 0)
            for (int k = 0; k < n; ++k) {
                const Node& a = vi_.nodes[reg.nodes[(k + n - 1) % n]];
                const Node& b = vi_.nodes[reg.nodes[k]];
                const Node& c = vi_.nodes[reg.nodes[(k + 1) % n]];
                const double ux = b.x - a.x, uy = b.y - a.y, vx = c.x - b.x, vy = c.y - b.y;
                const double lu = std::hypot(ux, uy), lv = std::hypot(vx, vy);
                if (lu == 0.0 || lv == 0.0) continue;
                const double turn = -(ux * vx + uy * vy) / (lu * lv);
                if (turn > sharpest) {
                    sharpest = turn;
                    start = k;
                }
            }
        if (start < 0) continue;
        std::vector<int> seq(n + 1);
        for (int k = 0; k <= n; ++k) seq[k] = reg.nodes[(start + k) % n];
        const int m = n + 1;

        auto devOf = [&](const LineFit& f, int node) {
            return std::fabs((vi_.nodes[node].y - f.cy) * f.dx - (vi_.nodes[node].x - f.cx) * f.dy);
        };
        // Whether the contour leaves the run end at position k (dir -1 start, +1 end) without a sharp
        // turn, judged from the node 3 steps further on.
        auto smoothEnd = [&](const LineFit& f, const std::vector<int>& seq, int k, int dir) {
            const Node& p = vi_.nodes[seq[k]];
            const Node& q = vi_.nodes[reg.nodes[((start + k + 3 * dir) % n + n) % n]];
            const double ux = q.x - p.x, uy = q.y - p.y, l = std::hypot(ux, uy);
            return l == 0.0 || std::fabs(ux * f.dy - uy * f.dx) / l < kLineSharpTurn;
        };

        // Greedy: from each position take the longest straight run, continue from its end.
        std::vector<std::pair<int, int>> runs;
        for (int i = 0; i + kLineMinNodes - 1 < m;) {
            int best = -1;
            for (int j = i + kLineMinNodes - 1; j < m; ++j) {
                const LineFit f = fitLine(vi_.nodes, seq, i, j);
                if (f.maxDev > kLineMaxDev || f.bow > kLineMaxBow) break;
                if (f.length >= kLineMinLength) best = j;
                if (!movable(seq[j])) break;
            }
            if (best < 0) {
                ++i;
                continue;
            }
            // An edge between two sharp corners is no piece of an arc, whatever its bow. A short run
            // that blends into curves at both ends is a flat step of a pixel staircase along a gentle
            // curve; straightened, it would lock the staircase in.
            const LineFit whole = fitLine(vi_.nodes, seq, i, best);
            const bool smoothA = smoothEnd(whole, seq, i, -1), smoothB = smoothEnd(whole, seq, best, 1);
            if ((whole.bow > kLineAcceptBow && (smoothA || smoothB)) ||
                (smoothA && smoothB && whole.length < kLineMinSmoothLength)) {
                ++i;
                continue;
            }
            int a = i, b = best;
            for (bool trimmed = true; trimmed && b - a + 1 > kLineMinNodes;) {
                const LineFit f = fitLine(vi_.nodes, seq, a, b);
                trimmed = false;
                if (smoothEnd(f, seq, a, -1) && devOf(f, seq[a]) > kLineEndDev) {
                    ++a;
                    trimmed = true;
                } else if (smoothEnd(f, seq, b, 1) && devOf(f, seq[b]) > kLineEndDev) {
                    --b;
                    trimmed = true;
                }
            }
            if (fitLine(vi_.nodes, seq, a, b).length >= kLineMinLength) runs.emplace_back(a, b);
            i = best;
        }
        if (runs.empty()) continue;

        // Up to two lines per contour position (a position shared by two runs is their corner).
        std::vector<LineFit> fits;
        std::vector<std::pair<int, int>> lineAt(m, {-1, -1});
        auto addLine = [&](int pos, int line) {
            if (lineAt[pos].first < 0) lineAt[pos].first = line;
            else lineAt[pos].second = line;
        };
        for (const auto& [a, b] : runs) {
            fits.push_back(fitLine(vi_.nodes, seq, a, b));
            LineFit& f = fits.back();
            if (std::fabs(f.dy) < kAxisSnapSin) {
                f.dx = f.dx < 0 ? -1.0 : 1.0;
                f.dy = 0.0;
            } else if (std::fabs(f.dx) < kAxisSnapSin) {
                f.dy = f.dy < 0 ? -1.0 : 1.0;
                f.dx = 0.0;
            }
            for (int k = a; k <= b; ++k) addLine(k, static_cast<int>(fits.size()) - 1);
        }
        // The start node appears at both ends of seq.
        if (lineAt[m - 1].first >= 0) addLine(0, lineAt[m - 1].first);
        lineAt[m - 1] = {-1, -1};

        for (int k = 0; k < m - 1; ++k) {
            const int node = seq[k];
            if (lineAt[k].first < 0 || !movable(node)) continue;
            onLine[node] = 1;
            Node& nd = vi_.nodes[node];
            const LineFit& f = fits[lineAt[k].first];
            if (lineAt[k].second < 0) {
                const double s = (nd.x - f.cx) * f.dx + (nd.y - f.cy) * f.dy;
                nd.x = f.cx + s * f.dx;
                nd.y = f.cy + s * f.dy;
                continue;
            }
            const LineFit& g = fits[lineAt[k].second];
            const double cross = f.dx * g.dy - f.dy * g.dx;
            if (std::fabs(cross) < kCornerMinSin) continue;
            const double t = ((g.cx - f.cx) * g.dy - (g.cy - f.cy) * g.dx) / cross;
            const double x = f.cx + t * f.dx, y = f.cy + t * f.dy;
            if (std::hypot(x - nd.x, y - nd.y) > kCornerMaxMove) continue;
            nd.x = x;
            nd.y = y;
        }
    }
    fairCurves(onLine);
}

void ContourSmoother::fairCurves(const std::vector<char>& onLine) {
    // Barriers: junctions (nodes with more than two contour neighbours; a node on the boundary of
    // two regions is in both contours with the same neighbours), sharp corners and border nodes.
    // They stay put and no smoothing window reaches past them. Nodes on straight lines stay put too,
    // but anchor the curves next to them.
    const int numNodes = static_cast<int>(vi_.nodes.size());
    std::vector<char> barrier(numNodes, 0);
    std::vector<int> nb(2 * numNodes, -1);
    auto addNeighbor = [&](int node, int other) {
        int* s = &nb[2 * node];
        if (s[0] == other || s[1] == other) return;
        if (s[0] < 0) s[0] = other;
        else if (s[1] < 0) s[1] = other;
        else barrier[node] = 1;
    };
    for (const Region& reg : vi_.regions) {
        const int n = static_cast<int>(reg.nodes.size());
        for (int k = 0; k < n; ++k) {
            const int node = reg.nodes[k];
            addNeighbor(node, reg.nodes[(k + n - 1) % n]);
            addNeighbor(node, reg.nodes[(k + 1) % n]);
            if (n < 5) {
                barrier[node] = 1;
                continue;
            }
            // A sharp turn, judged over two nodes either side, is a corner to keep.
            const Node& a = vi_.nodes[reg.nodes[(k + n - 2) % n]];
            const Node& b = vi_.nodes[node];
            const Node& c = vi_.nodes[reg.nodes[(k + 2) % n]];
            const double ux = b.x - a.x, uy = b.y - a.y, vx = c.x - b.x, vy = c.y - b.y;
            const double lu = std::hypot(ux, uy), lv = std::hypot(vx, vy);
            if (lu == 0.0 || lv == 0.0 || (ux * vx + uy * vy) / (lu * lv) < kFairCornerCos) barrier[node] = 1;
        }
    }
    const int numEdge = vi_.numEdgeNodes();
    for (int i = 0; i < numNodes; ++i)
        if (i < numEdge || (vi_.nodes[i].flags & 0x80) || vi_.nodes[i].type == 3 || nb[2 * i + 1] < 0) barrier[i] = 1;

    // Every free node moves onto the parabola fitted (least squares) through the kFairHalfWindow
    // nodes either side of it; a window ends at a barrier. A parabola follows a real curve, so the
    // curve does not shrink, but not a pixel staircase along it. Positions are read before any move.
    std::vector<double> nx(numNodes), ny(numNodes);
    std::vector<char> done(numNodes, 0);
    for (int i = 0; i < numNodes; ++i) {
        nx[i] = vi_.nodes[i].x;
        ny[i] = vi_.nodes[i].y;
    }
    for (const Region& reg : vi_.regions) {
        const int n = static_cast<int>(reg.nodes.size());
        for (int k = 0; k < n; ++k) {
            const int node = reg.nodes[k];
            if (barrier[node] || onLine[node] || done[node]) continue;
            done[node] = 1;
            int lo = 0, hi = 0;
            while (lo < kFairHalfWindow && !barrier[reg.nodes[(k - lo - 1 + n) % n]]) ++lo;
            if (lo < kFairHalfWindow) ++lo;  // the barrier itself still lies on the curve
            while (hi < kFairHalfWindow && !barrier[reg.nodes[(k + hi + 1) % n]]) ++hi;
            if (hi < kFairHalfWindow) ++hi;
            if (lo + hi + 1 >= n || lo < 2 || hi < 2) continue;
            const Node& first = vi_.nodes[reg.nodes[(k - lo + n) % n]];
            const Node& last = vi_.nodes[reg.nodes[(k + hi) % n]];
            double dx = last.x - first.x, dy = last.y - first.y;
            const double len = std::hypot(dx, dy);
            if (len < 1.0) continue;
            dx /= len;
            dy /= len;
            const Node& p = vi_.nodes[node];
            double A[3][3] = {}, B[3] = {};
            for (int j = -lo; j <= hi; ++j) {
                const Node& q = vi_.nodes[reg.nodes[(k + j + n) % n]];
                const double s = (q.x - p.x) * dx + (q.y - p.y) * dy, h = (q.y - p.y) * dx - (q.x - p.x) * dy;
                const double v[3] = {1.0, s, s * s};
                for (int r = 0; r < 3; ++r) {
                    for (int c = 0; c < 3; ++c) A[r][c] += v[r] * v[c];
                    B[r] += v[r] * h;
                }
            }
            // The fitted h at s = 0 (the node itself) by Cramer's rule.
            auto det3 = [](const double M[3][3]) {
                return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) -
                       M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
                       M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
            };
            const double d = det3(A);
            if (std::fabs(d) < 1e-12) continue;
            double C[3][3];
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c) C[r][c] = c == 0 ? B[r] : A[r][c];
            const double h0 = std::clamp(det3(C) / d, -kFairMaxMove, kFairMaxMove);
            nx[node] = p.x - h0 * dy;
            ny[node] = p.y + h0 * dx;
        }
    }
    for (int i = 0; i < numNodes; ++i) {
        vi_.nodes[i].x = nx[i];
        vi_.nodes[i].y = ny[i];
    }
}

namespace {

// isProbablyPuncturedNode: learned decision tree over the corner features (only
// features 0, 2, 6, 8, 15 and 16 are consulted).
bool isProbablyPuncturedNode(const double f[21]) {
    const double c = f[6];
    if (0.611194 <= c) {
        if (f[15] < 0.5) return 1.72002 <= c;
        if (!(f[16] < 0.097537)) return false;
        if (!(c < 1.42496)) return true;
        if (!(0.122841 <= f[2])) return true;
        if (!(0.445642 <= f[0])) return true;
        return f[8] < 0.0046405;
    }
    if (!(0.282594 <= c)) return false;
    if (!(f[0] < 0.322933)) return false;
    if (!(0.0286571 <= f[0])) return true;
    if (!(c < 0.483023)) return true;
    return f[2] < 0.00209524;
}

}  // namespace

void ContourSmoother::punctureCorners() {
    // Slide a 7-node window along every contour; the "turn" 2 − 2·cos at its five inner nodes
    // decides (via the learned tree) whether the centre node is a sharp corner. Corners become
    // fixed nodes (type 3) and lose the angular prior.
    double feat[21] = {};
    const int nr = static_cast<int>(vi_.regions.size());
    for (int r = 0; r < nr; ++r) {
        Region& reg = vi_.regions[r];
        const int n = static_cast<int>(reg.nodes.size());
        if (n <= 6) continue;
        int W[7];
        for (int j = 0; j < 7; ++j) W[j] = reg.nodes[(n - 3 + j) % n];
        double L[6], U[6][2], T[5];
        auto segment = [&](int from, int to, double& len, double u[2]) {
            const Node& a = vi_.nodes[from];
            const Node& b = vi_.nodes[to];
            double dx = b.x - a.x, dy = b.y - a.y;
            len = std::sqrt(dx * dx + 0.0 + dy * dy);
            if (len > 0.0) {
                dx = dx / len;
                dy = dy / len;
            }
            u[0] = dx;
            u[1] = dy;
        };
        auto turn = [](const double u[2], const double v[2]) {
            const double d = u[1] * v[1] + (u[0] * v[0] + 0.0);
            return 2.0 - (d + d);
        };
        for (int e = 0; e < 6; ++e) segment(W[e], W[e + 1], L[e], U[e]);
        for (int e = 0; e < 5; ++e) T[e] = turn(U[e], U[e + 1]);
        for (int k = 0;;) {
            for (int j = 2; j < 21; ++j) feat[j] = 0.0;
            feat[0] = ((T[0] + T[4]) + 0.0) + (T[1] + T[3]);
            feat[1] = (T[1] + T[3]) + 0.0;
            // Mirror so that the side with the larger inner turn comes first.
            const bool mirror = !(T[1] > T[3]);
            const double far = mirror ? T[0] : T[4], near = mirror ? T[1] : T[3], opp = mirror ? T[4] : T[0];
            feat[4] = far;
            feat[5] = near;
            feat[6] = T[2];
            feat[7] = mirror ? T[3] : T[1];
            feat[8] = opp;
            for (int e = 0; e < 6; ++e) feat[9 + e] = L[e];
            double mx = T[0];
            for (int e = 1; e < 5; ++e) mx = T[e] > mx ? T[e] : mx;
            feat[15] = (T[2] == mx) ? 1.0 : 0.0;
            if (T[2] == mx) reg.elems[k].flag |= 8;
            feat[3] = ((far + 0.0) + near) + 0.0;
            feat[2] = ((far + 0.0) + near) + opp;
            double mn = T[1] < T[0] ? T[1] : T[0];
            mn = T[3] < mn ? T[3] : mn;
            mn = T[4] < mn ? T[4] : mn;
            feat[16] = mn;
            if (isProbablyPuncturedNode(feat)) {
                vi_.nodes[W[3]].type = 3;
                reg.elems[k].flag |= 1;
            }
            // Slide the window by one node.
            const int next = reg.nodes[(k + 4) % n];
            for (int j = 0; j < 6; ++j) W[j] = W[j + 1];
            W[6] = next;
            double len, u[2];
            segment(W[5], W[6], len, u);
            for (int e = 0; e < 5; ++e) {
                L[e] = L[e + 1];
                U[e][0] = U[e + 1][0];
                U[e][1] = U[e + 1][1];
            }
            L[5] = len;
            U[5][0] = u[0];
            U[5][1] = u[1];
            for (int e = 0; e < 4; ++e) T[e] = T[e + 1];
            T[4] = 2.0 - ((u[1] * U[4][1] + (u[0] * U[4][0] + 0.0)) * 2.0);
            ++k;
            if (k == n) break;
        }
    }
}

}  // namespace isvg
