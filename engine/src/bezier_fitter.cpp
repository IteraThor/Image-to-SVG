#include "isvg/bezier_fitter.h"

#include <cmath>
#include <cstdlib>

#include "isvg/blas_ref.h"
#include "isvg/math_util.h"

namespace isvg {

namespace {

// sla::BezierCurve<double,2,3>::fitDataFixedEnds: least-squares inner control points
// for fixed ends (ridge 1e-5), via blas dgemm + LAPACK dgesv; returns the squared residual sum.
double fitCubicFixedEnds(BezierCurve& c, const std::vector<double>& t, const std::vector<double>& data, int n) {
    std::vector<double> b(2 * n), A(2 * n * 4, 0.0);
    double M[16] = {}, r[4] = {};
    const int lda = 2 * n;
    for (int i = 0; i < n; ++i) {
        const double ti = t[i], s = 1.0 - ti;
        const double t2 = ti * ti, s2 = s * s;
        const double t3 = ti * t2, s3 = s * s2;
        b[2 * i] = data[2 * i] - (c.p[0][0] * s3 + c.p[3][0] * t3);
        b[2 * i + 1] = data[2 * i + 1] - (s3 * c.p[0][1] + t3 * c.p[3][1]);
        const double b1 = (s2 * 3.0) * ti, b2 = (s * 3.0) * t2;
        A[2 * i + 0 * lda] = b1;
        A[2 * i + 2 * lda] = b2;
        A[2 * i + 1 + 1 * lda] = b1;
        A[2 * i + 1 + 3 * lda] = b2;
    }
    for (int i = 0; i < 4; ++i) M[i * 5] = 1e-05;
    blas::dgemmTN(4, 4, lda, 1.0, A.data(), lda, A.data(), lda, 1.0, M, 4);
    blas::dgemmTN(4, 1, lda, 1.0, A.data(), lda, b.data(), lda, 0.0, r, 4);
    int ipiv[4];
    blas::dgesv(4, 1, M, 4, ipiv, r, 4);
    c.p[1][0] = r[0];
    c.p[1][1] = r[1];
    c.p[2][0] = r[2];
    c.p[2][1] = r[3];
    double err = 0.0;
    for (int i = 0; i < n; ++i) {
        const double ti = t[i], s = 1.0 - ti;
        const double tt = ti * ti, ss = s * s;
        const double t3 = ti * tt;
        const double b1 = (ti * 3.0) * ss, b2 = (tt * 3.0) * s, b0 = s * ss;
        const double x = ((c.p[0][0] * b0 + b1 * c.p[1][0]) + b2 * c.p[2][0]) + c.p[3][0] * t3;
        const double y = ((b0 * c.p[0][1] + b1 * c.p[1][1]) + b2 * c.p[2][1]) + t3 * c.p[3][1];
        const double dx = data[2 * i] - x, dy = data[2 * i + 1] - y;
        err = err + (dy * dy + (dx * dx + 0.0));
    }
    return err;
}

// sla::BezierCurve<double,2,2>::fitDataFixedEnds (q: 3 control points, q[1] solved).
double fitQuadFixedEnds(double q[3][2], const std::vector<double>& t, const std::vector<double>& data, int n) {
    std::vector<double> b(2 * n), A(2 * n * 2, 0.0);
    double M[4] = {}, r[2] = {};
    const int lda = 2 * n;
    for (int i = 0; i < n; ++i) {
        const double ti = t[i], s = 1.0 - ti;
        const double s2 = s * s, b1 = (s + s) * ti, t2 = ti * ti;
        b[2 * i] = data[2 * i] - (q[0][0] * s2 + q[2][0] * t2);
        b[2 * i + 1] = data[2 * i + 1] - (s2 * q[0][1] + t2 * q[2][1]);
        A[2 * i] = b1;
        A[2 * i + 1 + lda] = b1;
    }
    M[0] = 1e-05;
    M[3] = 1e-05;
    blas::dgemmTN(2, 2, lda, 1.0, A.data(), lda, A.data(), lda, 1.0, M, 2);
    blas::dgemmTN(2, 1, lda, 1.0, A.data(), lda, b.data(), lda, 0.0, r, 2);
    int ipiv[2];
    blas::dgesv(2, 1, M, 2, ipiv, r, 2);
    q[1][0] = r[0];
    q[1][1] = r[1];
    double err = 0.0;
    for (int i = 0; i < n; ++i) {
        const double ti = t[i], s = 1.0 - ti;
        const double tt = ti * ti, b1 = (ti + ti) * s, ss = s * s;
        const double x = (q[0][0] * ss + b1 * q[1][0]) + q[2][0] * tt;
        const double y = (ss * q[0][1] + b1 * q[1][1]) + tt * q[2][1];
        const double dx = data[2 * i] - x, dy = data[2 * i + 1] - y;
        err = err + (dy * dy + (dx * dx + 0.0));
    }
    return err;
}

// sla::BezierCurve<double,2,D>::getCoefficients: power-basis coefficients of a Bézier of degree D.
void bezierToPower(const double (*p)[2], int degree, double (*coef)[2]) {
    int prevBinom = 1;
    for (int i = 0; i <= degree; ++i) {
        double binom;
        if (i == 0) {
            binom = 1.0;
            prevBinom = 1;
        } else {
            prevBinom = ((degree + 1 - i) * prevBinom) / i;
            binom = static_cast<double>(prevBinom);
        }
        double s = (i & 1) ? -1.0 : 1.0;
        double ax = 0.0, ay = 0.0;
        int c = 1;
        for (int j = 0;;) {
            ++j;
            ax = ax + p[j - 1][0] * s;
            ay = ay + p[j - 1][1] * s;
            if (i < j) break;
            c = (c * ((i + 1) - j)) / j;
            s = static_cast<double>(c);
            if ((j + i) & 1) s = -s;
        }
        coef[i][0] = binom * ax;
        coef[i][1] = binom * ay;
    }
}

// sla::BezierCurve<double,2,3>::setCoefficients: power basis → cubic control points.
void powerToCubic(const double coef[4][2], BezierCurve& c) {
    c.p[0][0] = coef[0][0];
    c.p[0][1] = coef[0][1];
    int binom = 1;
    for (int k = 1; k <= 3; ++k) {
        binom = ((4 - k) * binom) / k;
        double px = coef[k][0] * (1.0 / static_cast<double>(binom));
        double py = coef[k][1] * (1.0 / static_cast<double>(binom));
        double d = 1.0;
        int ck = 1;
        unsigned parity = static_cast<unsigned>(k);
        for (int j = 0; j < k; ++j) {
            if (parity & 1) d = -d;
            py = py - c.p[j][1] * d;
            px = px - d * c.p[j][0];
            if (j + 1 == k) break;
            ck = (ck * ((k + 1) - (j + 1))) / (j + 1);
            parity = static_cast<unsigned>(j + 1 + k);
            d = static_cast<double>(ck);
        }
        c.p[k][0] = px;
        c.p[k][1] = py;
    }
}

}  // namespace

// --- contour segmentation ---------------------------------------------------------------------

void BezierFitter::initialize() {
    int maxLen = -1;
    for (const Region& r : vi_.regions)
        if (maxLen < static_cast<int>(r.nodes.size())) maxLen = static_cast<int>(r.nodes.size());
    frags_.assign(vi_.regions.empty() ? 1 : maxLen + 2, Fragment{});
    for (Node& nd : vi_.nodes) {
        if (nd.type == 4) nd.type = 3;
        else if (nd.type == 1 || nd.type == 2 || nd.type == 5) nd.type = 0;
    }
}

int BezierFitter::followContour(int r, int k, int dir, uint8_t type) const {
    const Region& reg = vi_.regions[r];
    const int n = static_cast<int>(reg.nodes.size());
    const int self = reg.nodes[k % n];
    for (int i = k + dir + n;; i += dir) {
        const int j = i % n;
        const int node = reg.nodes[j];
        if (vi_.nodes[node].type != type || self == node) return j;
    }
}

void BezierFitter::initializeFragments(int r, int a, int b) {
    curRegion_ = r;
    curStart_ = a;
    runLength_ = b - a + 1;
    const int count = b - a + 1;
    if (count > static_cast<int>(frags_.size())) std::abort();
    for (int i = 0; i < count; ++i) {
        Fragment& f = frags_[i];
        f.start = i;
        f.len = 1;
        f.cost = 0.0;
        f.dirty = true;
        f.next = i + 1;
    }
    frags_[b - a].next = -1;
}

double BezierFitter::selfCost(int f) {
    Fragment& F = frags_[f];
    if (F.dirty) {
        const int s = curStart_ + F.start;
        F.cost = fitBezierCurve(scratch_, curRegion_, s, s + F.len, false);
        F.dirty = false;
    }
    return F.cost;
}

double BezierFitter::computeDislikeOf(int f, int ds, int de) {
    const int s = curStart_ + frags_[f].start;
    const int len = frags_[f].len;
    const double own = selfCost(f);
    double d = own - fitBezierCurve(scratch_, curRegion_, s + ds, len + s + de, false);
    if (ds < de) d = -d;
    return d;
}

double BezierFitter::computeMergeMetric(int f) {
    const int n = frags_[f].next;
    const int s = curStart_ + frags_[f].start;
    const int lf = frags_[f].len, ln = frags_[n].len;
    const double cf = selfCost(f);
    const double cn = selfCost(n);
    const double merged = fitBezierCurve(scratch_, curRegion_, s, ln + lf + s, false);
    return merged - (cf + cn);
}

void BezierFitter::considerMergeAndSwap(int f) {
    Fragment& F = frags_[f];
    const int n = F.next;
    Fragment& N = frags_[n];
    const double m = computeMergeMetric(f);
    if (thresh_ > m) {
        F.len += N.len;
        F.dirty = true;
        F.next = N.next;
        return;
    }
    const double d1 = computeDislikeOf(n, -1, 0);
    const double d2 = computeDislikeOf(f, 0, -1);
    if (0.0 > d1 - d2) {
        N.len += 1;
        F.len -= 1;
        N.dirty = true;
        N.start -= 1;
        F.dirty = true;
        if (F.start < N.start) return;
        std::abort();  // selectKeepers invariant violated
    }
    const double e1 = computeDislikeOf(f, 0, 1);
    const double e2 = computeDislikeOf(n, 1, 0);
    if (0.0 > e1 - e2) {
        F.len += 1;
        N.start += 1;
        N.len -= 1;
        F.dirty = true;
        N.dirty = true;
        if (N.next >= 0 && frags_[N.next].start > N.start) return;
        if (runLength_ > N.start) return;
        std::abort();
    }
}

void BezierFitter::mergeAndSwapLoop(int r, int a, int b) {
    initializeFragments(r, a, b);
    thresh_ = p_.statThreshInitial;
    const double l = math::log(p_.statThreshFinal / p_.statThreshInitial);
    const double rate = static_cast<double>(static_cast<float>(
        math::exp(static_cast<double>(static_cast<float>(l / (static_cast<double>(p_.maxIterations) * p_.rampFraction))))));
    for (int it = 0; it < p_.maxIterations; ++it) {
        for (int f = 0; f >= 0 && frags_[f].next >= 0 && frags_[frags_[f].next].next >= 0; f = frags_[f].next)
            considerMergeAndSwap(f);
        if (static_cast<double>(p_.maxIterations) * p_.rampFraction > static_cast<double>(it))
            thresh_ = thresh_ * rate;
        else
            thresh_ = p_.statThreshFinal;
    }
}

void BezierFitter::selectKeepers(int r, int a, int b) {
    mergeAndSwapLoop(r, a, b);
    const Region& reg = vi_.regions[r];
    const int n = static_cast<int>(reg.nodes.size());
    int idx = a;
    for (int f = 0; f >= 0;) {
        const Fragment& F = frags_[f];
        const int nextStart = F.next >= 0 ? frags_[F.next].start : runLength_;
        vi_.nodes[reg.nodes[idx % n]].type = 5;
        if (nextStart <= F.start + 1) {
            ++numCurves_;
            ++idx;
        } else {
            const int end = idx - F.start + nextStart;
            for (int i = idx + 1; i != end; ++i) {
                Node& nd = vi_.nodes[reg.nodes[i % n]];
                nd.i18 = numCurves_;
                nd.type = 1;
            }
            ++numCurves_;
            idx = end;
        }
        f = F.next;
    }
    vi_.nodes[reg.nodes[a % n]].type = 4;
    vi_.nodes[reg.nodes[b % n]].type = 4;
}

void BezierFitter::segmentContours() {
    numCurves_ = 0;
    for (int r = 0; r < static_cast<int>(vi_.regions.size()); ++r) {
        const Region& reg = vi_.regions[r];
        const int n = static_cast<int>(reg.nodes.size());
        for (int k = 0; k < n; ++k) {
            if (vi_.nodes[reg.nodes[k]].type != 0) continue;
            const int a = followContour(r, k, -1, 0);
            int b = followContour(r, k, 1, 0);
            while (b <= a) b += n;
            selectKeepers(r, a, b);
        }
    }
    curves_.assign(numCurves_, BezierCurve{});
}

// --- fitting -------------------------------------------------------------------------------------

int BezierFitter::computeDataAndTime(const BezierCurve& c, int r, int a, int b, std::vector<double>& t,
                                     std::vector<double>& data, bool final) {
    const Region& reg = vi_.regions[r];
    const int n = static_cast<int>(reg.nodes.size());
    while (b < a) b += n;
    const int count = b - a - 1;
    if (count < 1) return 0;
    t.resize(count);
    data.resize(2 * count);
    for (int i = 0; i < count; ++i) {
        Node& nd = vi_.nodes[reg.nodes[(a + 1 + i) % n]];
        const double x = nd.x, y = nd.y;
        data[2 * i] = x;
        data[2 * i + 1] = y;
        if (i == 0) {
            const double dx = x - c.p[0][0], dy = y - c.p[0][1];
            t[0] = std::sqrt(dy * dy + dx * dx + 0.0);
        } else {
            const double dx = x - data[2 * i - 2], dy = y - data[2 * i - 1];
            t[i] = std::sqrt(dy * dy + dx * dx + 0.0) + t[i - 1];
        }
        if (final) nd.type = 2;
    }
    const double ex = c.p[3][0] - data[2 * count - 2], ey = c.p[3][1] - data[2 * count - 1];
    double total = t[count - 1] + std::sqrt(ey * ey + ex * ex + 0.0);
    if (total == 0.0) total = 1.0;
    for (int i = 0; i < count; ++i) t[i] = t[i] / total;
    return count;
}

double BezierFitter::fitBezierCurve(BezierCurve& c, int r, int a, int b, bool final) {
    const Region& reg = vi_.regions[r];
    const int n = static_cast<int>(reg.nodes.size());
    while (b < a) b += n;
    auto pos = [&](int k, double out[2]) {
        const Node& nd = vi_.nodes[reg.nodes[k % n]];
        out[0] = nd.x;
        out[1] = nd.y;
    };
    pos(a, c.p[0]);
    pos(a + 1, c.p[1]);
    pos(n - 1 + b, c.p[2]);
    pos(b, c.p[3]);
    std::vector<double> t, data;
    const int count = computeDataAndTime(c, r, a, b, t, data, final);
    if (count < 1) return 0.0;
    if (count >= 4) return fitCubicFixedEnds(c, t, data, count);
    if (count == 3) {
        double q[3][2] = {{c.p[0][0], c.p[0][1]}, {0, 0}, {c.p[3][0], c.p[3][1]}};
        const double err = fitQuadFixedEnds(q, t, data, 3);
        double coef[4][2] = {};
        bezierToPower(q, 2, coef);
        powerToCubic(coef, c);
        return err;
    }
    // One or two inner points: a straight line with control points at quarters; the cost only
    // looks at the first inner point.
    const double p0x = c.p[0][0], p0y = c.p[0][1], p3x = c.p[3][0], p3y = c.p[3][1];
    c.p[1][0] = p0x * 0.75 + p3x * 0.25;
    c.p[1][1] = p0y * 0.75 + p3y * 0.25;
    c.p[2][0] = p3x * 0.75 + p0x * 0.25;
    c.p[2][1] = p3y * 0.75 + p0y * 0.25;
    const double t0 = t[0];
    const double ex = (p0x * (1.0 - t0) + p3x * t0) - data[0];
    const double ey = ((1.0 - t0) * p0y + t0 * p3y) - data[1];
    return ey * ey + ex * ex + 0.0;
}

void BezierFitter::fitBezierCurves() {
    for (int r = 0; r < static_cast<int>(vi_.regions.size()); ++r) {
        const Region& reg = vi_.regions[r];
        const int n = static_cast<int>(reg.nodes.size());
        for (int k = 0; k < n; ++k) {
            const Node& nd = vi_.nodes[reg.nodes[k]];
            if (nd.type != 1) continue;
            const int a = followContour(r, k, -1, 1);
            int b = followContour(r, k, 1, 1);
            while (b <= a) b += n;
            fitBezierCurve(curves_[nd.i18], r, a, b, true);
        }
    }
}

// --- curve references ----------------------------------------------------------------------------

void BezierFitter::buildCurveReferenceArrays() {
    const int nr = static_cast<int>(vi_.regions.size());
    refs_.assign(nr, {});
    for (int r = 0; r < nr; ++r) {
        const Region& reg = vi_.regions[r];
        const int n = static_cast<int>(reg.nodes.size());
        if (n <= 0) continue;
        // Start at a node that does not lie inside a curve.
        int k = 0, idx = 0;
        while (vi_.nodes[reg.nodes[idx]].type == 2) {
            --k;
            idx = (k + n) % n;
        }
        if (k != 0)
            do k += n;
            while (k < 0);
        const int end = n + k;
        std::vector<CurveReference>& out = refs_[r];
        while (true) {
            const int iPrev = (n + k - 1) % n;
            const int prevNode = reg.nodes[iPrev];
            int iCur = k % n;
            const int curNode = reg.nodes[iCur];
            const int iNext = (k + 1) % n;
            const int nextNode = reg.nodes[iNext];
            const int nb = reg.elems[iCur].neighbor;
            if (vi_.nodes[curNode].type == 2) {
                // A curve: walk to its end node.
                int e = iCur;
                while (vi_.nodes[reg.nodes[e]].type == 2) {
                    ++k;
                    e = k % n;
                }
                const int endNode = reg.nodes[e];
                const int ci = vi_.nodes[curNode].i18;
                const BezierCurve& c = curves_[ci];
                const Node& s = vi_.nodes[prevNode];
                const double eps = 9.999999974752427e-07;
                const bool reversed = std::fabs(c.p[0][0] - s.x) > eps || std::fabs(c.p[0][1] - s.y) > eps;
                out.push_back({ci, reversed ? 2 : 1, nb, prevNode, endNode, iPrev, e});
                if (end <= k) break;
                continue;
            }
            if (vi_.nodes[nextNode].type != 2) out.push_back({iCur, 0, nb, curNode, nextNode, iCur, iNext});
            k = k + 1;
            if (end <= k) break;
        }
    }
}

void BezierFitter::execute() {
    initialize();
    segmentContours();
    fitBezierCurves();
    buildCurveReferenceArrays();
}

}  // namespace isvg
