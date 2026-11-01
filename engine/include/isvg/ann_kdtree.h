// k-nearest-neighbour search as performed by the ANN library (ANN 1.1.1, Arya & Mount, LGPL-2.1), which
// the engine uses for VectorImage::buildIricNeighbors. Only what that call needs is
// reproduced: an ANNkd_tree with bucket size 1 and the sliding-midpoint split rule
// (ANN_KD_SUGGEST), and annkSearch. Tie-breaking (which of several equidistant points survive the
// k cut-off) depends on the exact tree layout and visiting order, so both follow ANN's code.
//
// This file is a port of parts of ANN and is distributed under the GNU LGPL 2.1.
#pragma once

#include <cfloat>
#include <memory>
#include <vector>

namespace isvg::ann {

class KdTree {
public:
    // pts: n 2-D points (x, y); kept by reference (ANN stores point pointers).
    KdTree(const std::vector<const double*>& pts);
    ~KdTree();
    // k nearest points to q (dist² ascending, ties in ANN order); eps as in annkSearch.
    void kSearch(const double* q, int k, int* idx, double* dist, double eps) const;

private:
    struct Node;
    Node* build(int* pidx, int n, double lo[2], double hi[2]);

    const std::vector<const double*>& pts_;
    std::vector<int> pidx_;
    double bndLo_[2] = {0, 0}, bndHi_[2] = {0, 0};
    std::vector<std::unique_ptr<Node>> nodes_;
    Node* root_ = nullptr;
};

}  // namespace isvg::ann
