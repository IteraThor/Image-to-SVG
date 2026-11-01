// BezierFitter: splits every contour into runs between
// fixed nodes, merges the pixel steps of each run into curve fragments (annealed merge/swap on the
// least-squares fit cost) and fits one cubic Bézier per fragment.
#pragma once

#include <vector>

#include "isvg/vector_image.h"

namespace isvg {

struct BezierCurve {  // sla::BezierCurve<double,2,3> (0x40 bytes)
    double p[4][2] = {};
};

// CurveReference (0x1c bytes, Region: one piece of a region's outline.
struct CurveReference {
    int index = -1;  // curve index (types 1/2) or contour index of the step start (type 0)
    int type = 0;  // 0 straight step, 1 Bézier forward, 2 Bézier reversed
    int neighbor = -1;  // region on the other side
    int startNode = -1, endNode = -1;
    int startIdx = -1, endIdx = -1;  // contour indices
};

struct BezierFitterParams {  // EP_LOGO_HIG_AA
    double statThreshInitial = 1e-8;
    double statThreshFinal = 0.2;
    double rampFraction = 0.5;
    int maxIterations = 20;
};

class BezierFitter {
public:
    BezierFitter(VectorImage& vi, const BezierFitterParams& p) : vi_(vi), p_(p) {}
    void execute();

    std::vector<BezierCurve>& curves() { return curves_; }
    std::vector<std::vector<CurveReference>>& references() { return refs_; }
    int numCurves() const { return numCurves_; }

private:
    struct Fragment {  // 0x20 bytes
        int start = 0, len = 1;
        int next = -1;
        double cost = 0;
        bool dirty = true;
    };

    void initialize();
    void segmentContours();
    int followContour(int r, int k, int dir, uint8_t type) const;
    void selectKeepers(int r, int a, int b);
    void mergeAndSwapLoop(int r, int a, int b);
    void initializeFragments(int r, int a, int b);
    void considerMergeAndSwap(int f);
    double selfCost(int f);
    double computeDislikeOf(int f, int ds, int de);
    double computeMergeMetric(int f);
    void fitBezierCurves();
    double fitBezierCurve(BezierCurve& c, int r, int a, int b, bool final);
    int computeDataAndTime(const BezierCurve& c, int r, int a, int b, std::vector<double>& t,
                           std::vector<double>& data, bool final);
    void buildCurveReferenceArrays();

    VectorImage& vi_;
    BezierFitterParams p_;
    std::vector<Fragment> frags_;
    int numCurves_ = 0;
    int runLength_ = 0;
    double thresh_ = 0;
    BezierCurve scratch_;
    int curRegion_ = 0, curStart_ = 0;
    std::vector<BezierCurve> curves_;  // VectorImage
    std::vector<std::vector<CurveReference>> refs_;  // per region
};

}  // namespace isvg
