// Port of the parts of ANN 1.1.1 (kd_tree.cpp, kd_split.cpp, kd_util.cpp, kd_search.cpp,
// pr_queue_k.h) used by the engine. Copyright (c) 1997-2005 University of Maryland and
// Sunil Arya and David Mount. Distributed under the GNU Lesser General Public License 2.1.
#include "isvg/ann_kdtree.h"

#include <cstdlib>

namespace isvg::ann {

namespace {

constexpr double kSplitErr = 0.001;  // kd_split.cpp ERR

// ANNmin_k: the k smallest keys seen so far, ascending; equal keys keep insertion order and a
// key equal to the current maximum of a full set is dropped.
struct MinK {
    struct Item {
        double key;
        int info;
    };
    int k, n = 0;
    std::vector<Item> mk;
    explicit MinK(int kk) : k(kk), mk(kk + 1) {}
    double maxKey() const { return n == k ? mk[k - 1].key : DBL_MAX; }
    void insert(double kv, int inf) {
        int i;
        for (i = n; i > 0; --i) {
            if (mk[i - 1].key > kv) mk[i] = mk[i - 1];
            else break;
        }
        mk[i].key = kv;
        mk[i].info = inf;
        if (n < k) ++n;
    }
};

struct SearchCtx {
    const double* q;
    double maxErr;
    const std::vector<const double*>* pts;
    MinK* mk;
};

}  // namespace

struct KdTree::Node {
    bool leaf = false;
    // leaf
    int nPts = 0;
    const int* bkt = nullptr;
    // split
    int cutDim = 0;
    double cutVal = 0, cdLo = 0, cdHi = 0;
    Node* lo = nullptr;
    Node* hi = nullptr;

    void search(double boxDist, SearchCtx& s) const {
        if (leaf) {
            double minDist = s.mk->maxKey();
            for (int i = 0; i < nPts; ++i) {
                const double* pp = (*s.pts)[bkt[i]];
                double dist = 0;
                int d;
                for (d = 0; d < 2; ++d) {
                    const double t = s.q[d] - pp[d];
                    if ((dist = dist + t * t) > minDist) break;
                }
                if (d >= 2) {  // ANN_ALLOW_SELF_MATCH
                    s.mk->insert(dist, bkt[i]);
                    minDist = s.mk->maxKey();
                }
            }
            return;
        }
        const double cutDiff = s.q[cutDim] - cutVal;
        if (cutDiff < 0) {
            lo->search(boxDist, s);
            double boxDiff = cdLo - s.q[cutDim];
            if (boxDiff < 0) boxDiff = 0;
            boxDist = boxDist + (cutDiff * cutDiff - boxDiff * boxDiff);  // ANN_DIFF(x,y) = y - x
            if (boxDist * s.maxErr < s.mk->maxKey()) hi->search(boxDist, s);
        } else {
            hi->search(boxDist, s);
            double boxDiff = s.q[cutDim] - cdHi;
            if (boxDiff < 0) boxDiff = 0;
            boxDist = boxDist + (cutDiff * cutDiff - boxDiff * boxDiff);  // ANN_DIFF(x,y) = y - x
            if (boxDist * s.maxErr < s.mk->maxKey()) lo->search(boxDist, s);
        }
    }
};

KdTree::KdTree(const std::vector<const double*>& pts) : pts_(pts) {
    const int n = static_cast<int>(pts.size());
    pidx_.resize(n);
    for (int i = 0; i < n; ++i) pidx_[i] = i;
    if (n == 0) return;
    // annEnclRect
    for (int d = 0; d < 2; ++d) {
        double lo = pts_[pidx_[0]][d], hi = pts_[pidx_[0]][d];
        for (int i = 0; i < n; ++i) {
            const double v = pts_[pidx_[i]][d];
            if (v < lo) lo = v;
            else if (v > hi) hi = v;
        }
        bndLo_[d] = lo;
        bndHi_[d] = hi;
    }
    double lo[2] = {bndLo_[0], bndLo_[1]}, hi[2] = {bndHi_[0], bndHi_[1]};
    root_ = build(pidx_.data(), n, lo, hi);
}

KdTree::~KdTree() = default;

KdTree::Node* KdTree::build(int* pidx, int n, double lo[2], double hi[2]) {
    nodes_.push_back(std::make_unique<Node>());
    Node* node = nodes_.back().get();
    if (n <= 1) {  // bucket size 1 (n == 0 → the shared trivial leaf; an empty leaf here)
        node->leaf = true;
        node->nPts = n;
        node->bkt = pidx;
        return node;
    }
    auto PA = [&](int i, int d) { return pts_[pidx[i]][d]; };
    // sl_midpt_split
    double maxLength = hi[0] - lo[0];
    for (int d = 1; d < 2; ++d) {
        const double length = hi[d] - lo[d];
        if (length > maxLength) maxLength = length;
    }
    double maxSpread = -1;
    int cd = 0;
    for (int d = 0; d < 2; ++d) {
        if ((hi[d] - lo[d]) >= (1 - kSplitErr) * maxLength) {
            double mn = PA(0, d), mx = PA(0, d);  // annSpread
            for (int i = 1; i < n; ++i) {
                const double c = PA(i, d);
                if (c < mn) mn = c;
                else if (c > mx) mx = c;
            }
            const double spr = mx - mn;
            if (spr > maxSpread) {
                maxSpread = spr;
                cd = d;
            }
        }
    }
    const double ideal = (lo[cd] + hi[cd]) / 2;
    double mn = PA(0, cd), mx = PA(0, cd);  // annMinMax
    for (int i = 1; i < n; ++i) {
        const double c = PA(i, cd);
        if (c < mn) mn = c;
        else if (c > mx) mx = c;
    }
    double cv;
    if (ideal < mn) cv = mn;
    else if (ideal > mx) cv = mx;
    else cv = ideal;
    // annPlaneSplit
    auto swap = [&](int a, int b) {
        const int t = pidx[a];
        pidx[a] = pidx[b];
        pidx[b] = t;
    };
    int l = 0, r = n - 1;
    for (;;) {
        while (l < n && PA(l, cd) < cv) ++l;
        while (r >= 0 && PA(r, cd) >= cv) --r;
        if (l > r) break;
        swap(l, r);
        ++l;
        --r;
    }
    const int br1 = l;
    r = n - 1;
    for (;;) {
        while (l < n && PA(l, cd) <= cv) ++l;
        while (r >= br1 && PA(r, cd) > cv) --r;
        if (l > r) break;
        swap(l, r);
        ++l;
        --r;
    }
    const int br2 = l;
    int nLo;
    if (ideal < mn) nLo = 1;
    else if (ideal > mx) nLo = n - 1;
    else if (br1 > n / 2) nLo = br1;
    else if (br2 < n / 2) nLo = br2;
    else nLo = n / 2;

    const double lv = lo[cd], hv = hi[cd];
    hi[cd] = cv;
    Node* loChild = build(pidx, nLo, lo, hi);
    hi[cd] = hv;
    lo[cd] = cv;
    Node* hiChild = build(pidx + nLo, n - nLo, lo, hi);
    lo[cd] = lv;
    node->cutDim = cd;
    node->cutVal = cv;
    node->cdLo = lv;
    node->cdHi = hv;
    node->lo = loChild;
    node->hi = hiChild;
    return node;
}

void KdTree::kSearch(const double* q, int k, int* idx, double* dist, double eps) const {
    if (k > static_cast<int>(pts_.size())) std::abort();
    MinK mk(k);
    SearchCtx s{q, (1.0 + eps) * (1.0 + eps), &pts_, &mk};
    // annBoxDistance
    double bd = 0.0;
    for (int d = 0; d < 2; ++d) {
        if (q[d] < bndLo_[d]) {
            const double t = bndLo_[d] - q[d];
            bd = bd + t * t;
        } else if (q[d] > bndHi_[d]) {
            const double t = q[d] - bndHi_[d];
            bd = bd + t * t;
        }
    }
    if (root_) root_->search(bd, s);
    for (int i = 0; i < k; ++i) {
        dist[i] = i < mk.n ? mk.mk[i].key : DBL_MAX;
        idx[i] = i < mk.n ? mk.mk[i].info : -1;
    }
}

}  // namespace isvg::ann
