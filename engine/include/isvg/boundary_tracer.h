// BoundaryTracer: turns the final segment map into a
// graph of pixel-corner nodes and one closed outer contour per region.
#pragma once

#include <vector>

#include "isvg/pixel_map.h"
#include "isvg/pixel_segmentation.h"
#include "isvg/vector_image.h"

namespace isvg {

// Segmenter::buildRegions: one Region per segment id 0..maxId with its pixel count,
// palette index and byte colour (PixelSegmentation::getColorUc: (int)(c*255f + 0.5f)).
void buildRegions(PixelSegmentation& seg, std::vector<Region>& regions);

class BoundaryTracer {
public:
    BoundaryTracer(PixelMap& map, VectorImage& image, bool antiAliased)
        : map_(map), regions_(image.regions), image_(image), antiAliased_(antiAliased) {}

    void execute();

private:
    int enumerateBoundaryNodes();
    void traceBoundary(int seg, int x, int y);
    void removeRedundantNodes(bool antiAliased);

    int id(int x, int y) const { return map_.rawIds()[y * w_ + x]; }
    int& nodeIndex(int x, int y) { return nodeIndex_[y * (w_ + 1) + x]; }
    int& visited(int x, int y) { return visited_[y * (w_ + 1) + x]; }

    PixelMap& map_;
    std::vector<Region>& regions_;
    VectorImage& image_;  // (nodes:
    bool antiAliased_;  // config [4]
    int w_ = 0, h_ = 0;
    std::vector<int> nodeIndex_;  //: (w+1)×(h+1) node index per pixel corner, -1 = none
    std::vector<int> visited_;  //: (w+1)×(h+1) region that last passed through a corner
};

}  // namespace isvg
