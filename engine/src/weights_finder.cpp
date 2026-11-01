#include "isvg/weights_finder.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "isvg/blas_ref.h"

namespace isvg {
namespace {

constexpr int kLd = 9;
constexpr int kN8i[8] = {1, 0, -1, 0, -1, 1, -1, 1};
constexpr int kN8j[8] = {0, 1, 0, -1, 1, -1, -1, 1};
constexpr int kN24i[24] = {1, 0, -1, 0, -1, 1, -1, 1, -2, -2, -2, -2, -2, 2, 2, 2, 2, 2, -1, 0, 1, -1, 0, 1};
constexpr int kN24j[24] = {0, 1, 0, -1, 1, -1, -1, 1, -2, -1, 0, 1, 2, -2, -1, 0, 1, 2, -2, -2, -2, 2, 2, 2};

inline double& at(double* M, int r, int c) { return M[r + c * kLd]; }

[[noreturn]] void fatal(const char* what) {
    std::fprintf(stderr, "WeightsFinder: %s\n", what);
    std::abort();
}

// In-place lower Cholesky, row by row, as in setupQP / solveQP.
void cholesky(double* M, int n) {
    double d = at(M, 0, 0);
    for (int j = 0;; ++j) {
        if (!(d > 0.0)) fatal("matrix not positive definite");
        at(M, j, j) = std::sqrt(d);
        if (j + 1 == n) return;
        const int i = j + 1;
        for (int k = 0; k <= j; ++k) {
            double s = at(M, i, k);
            for (int m = 0; m < k; ++m) s = s - at(M, i, m) * at(M, k, m);
            at(M, i, k) = s / at(M, k, k);
        }
        double s = at(M, i, i);
        for (int m = 0; m <= j; ++m) s = s - at(M, i, m) * at(M, i, m);
        at(M, i, i) = s;
        d = s;
    }
}

// Solve L Lᵀ x = b in place. Forward terms are subtracted from the highest index down; backward
// terms from the lowest index up (fixed loop order, for reproducible results).
void cholSolve(double* L, int n, double* x) {
    for (int i = 0; i < n; ++i) {
        double s = x[i];
        for (int m = i - 1; m >= 0; --m) s = s - at(L, i, m) * x[m];
        x[i] = s / at(L, i, i);
    }
    for (int i = n - 1; i >= 0; --i) {
        double s = x[i];
        for (int r = i + 1; r < n; ++r) s = s - at(L, r, i) * x[r];
        x[i] = s / at(L, i, i);
    }
}

}  // namespace

void WeightsFinder::clear() {
    k_ = 0;
    ids_.fill(0);
    for (double& v : w_) v = 0.0;
    for (double& v : A_) v = 0.0;
    for (double& v : H_) v = 0.0;
    for (double& v : f_) v = 0.0;
    for (double& v : E_) v = 0.0;
    for (int i = 0; i < 9; ++i) E_[i] = 1.0;  // constraint 0: Σw = 1
    for (double& v : Y_) v = 0.0;
    for (double& v : S_) v = 0.0;
    for (double& v : lambda_) v = 0.0;
    for (double& v : r_) v = 0.0;
    for (bool& b : pinned_) b = false;
    for (double& v : p_) v = 0.0;
    numCandidates_ = -1;
}

int WeightsFinder::idWithReplace(int x, int y, int replaceTo, int replaceFrom) {
    const int id = map_->idOfReplace(x, y, replaceTo, replaceFrom);
    return (id >= 0 && idMap_) ? idMap_(id) : id;
}

int WeightsFinder::updateNeighbors(int x, int y, int replaceTo, int replaceFrom, bool use24) {
    const int own = idWithReplace(x, y, replaceTo, replaceFrom);
    ids_[0] = own;
    int k = own >= 0 ? 1 : 0;
    const int n = use24 ? 24 : 8;
    for (int i = 0; i < n; ++i) {
        const int id = use24 ? idWithReplace(x + kN24i[i], y + kN24j[i], replaceTo, replaceFrom)
                             : idWithReplace(x + kN8i[i], y + kN8j[i], replaceTo, replaceFrom);
        if (id < 0) continue;
        bool seen = false;
        for (int j = 0; j < k; ++j)
            if (ids_[j] == id) { seen = true; break; }
        if (seen) continue;
        ids_[k++] = id;
        if (k > 8) return k;
    }
    return k;
}

void WeightsFinder::setupQP() {
    for (int i = 0; i < k_; ++i) {
        const float* c = colorOf_(ids_[i]);
        A_[3 * i + 0] = static_cast<double>(c[0]);
        A_[3 * i + 1] = static_cast<double>(c[1]);
        A_[3 * i + 2] = static_cast<double>(c[2]);
    }
    blas::dgemmTN(k_, k_, 3, 1.0, A_, 3, A_, 3, 0.0, H_, kLd);
    for (int i = 0; i < k_; ++i) at(H_, i, i) = at(H_, i, i) + 0.01;
    cholesky(H_, k_);
}

void WeightsFinder::solveQP() {
    blas::dgemv('T', 3, k_, 1.0, A_, 3, p_, 0.0, f_);
    for (bool& b : pinned_) b = false;
    int t = 0;
    for (;;) {
        const int m = t + 1;
        for (int c = 0; c <= t; ++c) {
            double* y = &Y_[c * kLd];
            for (int i = 0; i < k_; ++i) y[i] = E_[c * kLd + i];
            cholSolve(H_, k_, y);
        }
        blas::dgemmTN(m, m, k_, 1.0, E_, kLd, Y_, kLd, 0.0, S_, kLd);
        blas::dgemv('T', k_, m, 1.0, Y_, kLd, f_, 0.0, r_);
        r_[0] = r_[0] - 1.0;
        cholesky(S_, m);
        for (int i = 0; i < m; ++i) lambda_[i] = r_[i];
        cholSolve(S_, m, lambda_);
        for (int i = 0; i < k_; ++i) w_[i] = f_[i];
        blas::dgemv('N', k_, m, -1.0, E_, kLd, lambda_, 1.0, w_);
        cholSolve(H_, k_, w_);

        int worst = -1;
        double minW = 0.0;
        for (int i = 0; i < k_; ++i) {
            if (w_[i] < minW && !pinned_[i]) {
                worst = i;
                minW = w_[i];
            }
        }
        if (worst < 0 || !(w_[worst] < -0.0001)) break;
        ++t;
        if (k_ <= t) fatal("too many active constraints");
        for (int i = 0; i < k_; ++i) E_[t * kLd + i] = 0.0;
        E_[t * kLd + worst] = 1.0;
        pinned_[worst] = true;
    }
    blas::dgemv('N', 3, k_, 1.0, A_, 3, w_, 0.0, p_);
}

Color4f WeightsFinder::computePixelF(int x, int y, int replaceTo, int replaceFrom, double* ownWeight,
                                     bool use24) {
    const int k = updateNeighbors(x, y, replaceTo, replaceFrom, use24);
    k_ = k;
    Color4f out;
    if (k >= 10) fatal("too many neighbours");
    if (k < 2) {
        if (ownWeight) *ownWeight = 1.0;
        const float* c = colorOf_(ids_[0]);
        for (int i = 0; i < 4; ++i) out.c[i] = c[i];
        return out;
    }
    if (x < 0 || y < 0 || x >= img_->w || y >= img_->h) fatal("pixel outside image");
    const uint8_t* px = img_->at(y * img_->w + x);
    float p[4];
    for (int i = 0; i < 4; ++i) p[i] = static_cast<float>(px[i]) * kInv255f;
    if (k == 2) {
        const float* c0 = colorOf_(ids_[0]);
        const float* c1 = colorOf_(ids_[1]);
        float d[4], e[4], dd = 0.0f, dot = 0.0f;
        for (int i = 0; i < 4; ++i) d[i] = c1[i] - c0[i];
        dd = d[0] * d[0] + dd;
        dd = d[1] * d[1] + dd;
        dd = dd + d[2] * d[2];
        dd = dd + d[3] * d[3];
        if (!(dd > 0.0f)) {
            if (ownWeight) *ownWeight = 0.5;
            w_[0] = w_[1] = 0.5;
            for (int i = 0; i < 4; ++i) out.c[i] = c0[i];
            return out;
        }
        for (int i = 0; i < 4; ++i) e[i] = (p[i] - c0[i]) * d[i];
        dot = 0.0f + e[0];
        dot = e[1] + dot;
        dot = dot + e[2];
        dot = dot + e[3];
        const double ddd = static_cast<double>(dd);
        const double t0 = static_cast<double>(dot) / ddd;
        const double shrink = 0.01 / (ddd + 0.02);
        double t = (1.0 - (t0 + t0)) * shrink + t0;
        double w0;
        float f;
        if (0.0 > t) {
            w0 = 1.0;
            f = 0.0f;
            t = 0.0;
        } else if (t > 1.0) {
            f = 1.0f;
            t = 1.0;
            w0 = 0.0;
        } else {
            w0 = 1.0 - t;
            f = static_cast<float>(t);
        }
        w_[0] = w0;
        w_[1] = t;
        if (ownWeight) *ownWeight = w0;
        for (int i = 0; i < 4; ++i) out.c[i] = f * d[i] + c0[i];
        return out;
    }
    setupQP();
    for (int i = 0; i < 4; ++i) p_[i] = static_cast<double>(p[i]);
    solveQP();
    if (ownWeight) *ownWeight = w_[0];
    for (int i = 0; i < 4; ++i) out.c[i] = static_cast<float>(p_[i]);
    return out;
}

int WeightsFinder::heaviestSegmentId() const {
    if (numCandidates_ < 2) return ids_[0];
    int best = 0;
    for (int i = 1; i < numCandidates_; ++i)
        if (w_[best] <= w_[i] && w_[i] != w_[best]) best = i;
    return ids_[best];
}

}  // namespace isvg
