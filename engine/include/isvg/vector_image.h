// Geometry produced by the boundary tracer and refined by the later stages
// VectorImage at: regions at, nodes at.
#pragma once

#include <climits>
#include <cstdint>
#include <utility>
#include <vector>

namespace isvg {

// Node (0x30 bytes): a pixel-corner position on a region boundary.
struct Node {
    double x = 0, y = 0;  // current position (optimised)
    float fx = 0, fy = 0;  // initial (pixel-corner) position
    int i18 = -1;
    uint8_t type = 0;  //: 3 = fixed (image corner or junction of >= 3 regions)
    //: 0x01 IRIC neighbours built, 0x02 inverted, 0x04 checked, 0x10 use quadratic prior,
    // 0x40 redundant (removed), 0x80 on the image border
    uint8_t flags = 0;
    uint8_t antiInv = 0;  // anti-inversion penalty count-down (b1e)
    //: nearby contour elements (region, contour index) from buildIricNeighbors
    std::vector<std::pair<int, int>> iric;
};

// ContourElement (12 bytes): one step of a region contour, from the previous node to this one.
struct ContourElement {
    uint8_t flag = 0;
    int neighbor = -1;  //: region on the other side of the step (-1 = image border)
    int steps = 1;  //: pixel edges covered since the previous kept node
};

// Region (0x90 bytes): one final segment with its outer contour (holes are the contours of the
// regions it contains).
struct Region {
    int n = 0;  // pixel count
    double area = 0;  // signed contour area (VectorImage::updateArea)
    int paletteIndex = 0;
    uint8_t color[4] = {0, 0, 0, 0};  // (image channel order)
    std::vector<int> nodes;  // contour node indices
    std::vector<ContourElement> elems;  // alongside `nodes`
    std::vector<int> neighbors;  // sorted, without -1
    int parent = -1;  // enclosing region
    double insideNormal[2] = {0, 0};
    int baseSignOffset = INT_MAX;
};

struct VectorImage {
    std::vector<Region> regions;
    std::vector<Node> nodes;
    int edgeNodeCounts[4] = {};  //: nodes on the left, right, top, bottom image side
    int numNodesInitial = 0;
    int w = 0, h = 0;
    std::vector<int> marginNodes;

    int numEdgeNodes() const {
        return edgeNodeCounts[0] + edgeNodeCounts[1] + edgeNodeCounts[2] + edgeNodeCounts[3];
    }

    // VectorImage::setNodesSize (fresh default nodes).
    void setNodesSize(int n) {
        nodes.assign(n, Node{});
        numNodesInitial = n;
    }
    //: the left and right columns hold h+1 nodes each (corners included), the top
    // and bottom rows the w-1 nodes between the corners.
    void initializeEdgeNodeCounts(int width, int height) {
        w = width;
        h = height;
        edgeNodeCounts[0] = height + 1;
        edgeNodeCounts[1] = height + 1;
        edgeNodeCounts[2] = width - 1;
        edgeNodeCounts[3] = width - 1;
    }
    //: image corners are fixed (type 3), every other node starts as type 0.
    void setInitialNodePositionAndType(int i, int x, int y) {
        Node& nd = nodes[i];
        nd.x = x;
        nd.y = y;
        nd.type = ((x == 0 || x == w) && (y == 0 || y == h)) ? 3 : 0;
        nd.fx = static_cast<float>(nd.x);
        nd.fy = static_cast<float>(nd.y);
    }

    // Offset of the node from its initial position.
    void movedVector(int i, double out[2]) const {
        const Node& nd = nodes[i];
        out[0] = nd.x - static_cast<double>(nd.fx);
        out[1] = nd.y - static_cast<double>(nd.fy);
    }

    void initializeMarginNodes(int width, int height);
    void clampEdgeNodes();
    void clampMarginNodes();
    void updateArea(int r);
    int decrementNodeCounts();

    // Inversion (self-/mutual-intersection) detection, anti-aliased mode only
    // (src/vector_image_inversion.cpp).
    void buildIricNeighbors();  // (ANN kd-tree, k=50, r=5)
    int flagInversionNodes();
    int flagInversionNodesForRegion(int r);
    int computeBaseSignOffset(int r, int start);
    void computeInsideNormal(int r, int k, double out[2]) const;
};

}  // namespace isvg
