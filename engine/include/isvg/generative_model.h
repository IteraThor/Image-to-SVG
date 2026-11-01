// GenerativeModel: renders the region contours as an
// anti-aliased image by exact area coverage and measures it against the input image.
//
// Every contour edge is walked across the pixel grid ("bread crumbs" at the points where it enters,
// crosses and leaves pixel cells). Per pixel, the crumbs of each region form a polygon once the
// cell corners between them are inserted; the polygon areas weight the region colours.
#pragma once

#include <cstdint>
#include <set>
#include <vector>

#include "isvg/image.h"
#include "isvg/pixel_map.h"
#include "isvg/vector_image.h"

namespace isvg {

struct BreadCrumb {          
    int next = -1;  // (pointer)
    uint32_t flags = 0;  //: 1 on the cell boundary, 2 cell corner, 4 contour node,
                             // 0x10 end node, 0x80 horizontal / 0x100 vertical grid line
    int k = 0;  // contour index of the edge (-2 = inserted cell corner)
    int region = 0;
    float x = 0, y = 0;
    float perim = -1;  // position along the cell boundary, 0..4 from the top-left
    float winding = 0;  // turning angle of the contour so far
};

class GenerativeModel {
public:
    GenerativeModel(VectorImage& vi, PixelMap& map, const Image& img) : vi_(vi), map_(map), img_(img) {}

    void setup();
    double findPotential();
    void findGradient(std::vector<double>& grad);  // (Vec2 per node, added)
    bool checkFlaggedRegionsForInversions();

    const std::vector<float>& model() const { return model_; }
    const std::vector<uint8_t>& mask() const { return mask_; }

private:
    void dropBreadCrumbs(int region, int k, int kNext);
    void dropBreadCrumb(int ix, int iy, int region, int k, uint32_t flags, double x, double y);
    void closePixel(int px, int py);
    int checkInsertCornerNodes(int c, int px, int py);
    int checkCloseRegion(int last, int first, int px, int py);
    int insertCornerNodes(int c, int n, int px, int py);
    int insertCornerNode(int c, int n, int px, int py, int corner);
    void cleanOutTwoCrumbRegions(int px, int py);
    double computePixelColor(int px, int py);
    void findPixelGradient(std::vector<double>& grad, int px, int py);
    int regionToIsland(int region) const;
    void regionColorF(int r, float out[4]) const;
    void observedF(int px, int py, float out[4]) const;
    void resetPixelToObserved(int px, int py);

    int& cellHead(int px, int py) { return cellHead_[py * w_ + px]; }
    float* modelAt(int px, int py) { return &model_[(static_cast<size_t>(py) * w_ + px) * 4]; }
    int newCrumb() {
        crumbs_.emplace_back();
        return static_cast<int>(crumbs_.size()) - 1;
    }

    VectorImage& vi_;  // (regions, nodes
    PixelMap& map_;
    const Image& img_;
    int w_ = 0, h_ = 0;
    std::vector<BreadCrumb> crumbs_;  // BreadCrumbContainer
    std::vector<int> cellHead_, cellTail_;  // (head, tail) per pixel
    std::vector<float> model_;  // model colour per pixel
    std::vector<uint8_t> mask_;  // 1 = negative area, 2 = unresolved coverage
    double winding_ = 0;
    // Per pixel scratch regions, areas, count, island parent). The
    // area array has room for 10; index `count` can be touched once (original quirk), hence 11.
    int slotRegion_[10] = {};
    double slotArea_[11] = {};
    int slotCount_ = 0;
    int islandParent_ = -1;
    std::set<int> flagged_;  // regions to check for inversions
};

}  // namespace isvg
