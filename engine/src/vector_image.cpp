#include "isvg/vector_image.h"

namespace isvg {

void VectorImage::initializeMarginNodes(int width, int height) {
    // Non-border nodes within 1.1 px of the image edge are clamped after every optimiser move.
    w = width;
    h = height;
    const double lo = 1.100000023841858;  // double constant
    const float m = 1.1f;
    const double hiX = static_cast<double>(static_cast<float>(w) - m);
    const double hiY = static_cast<double>(static_cast<float>(h) - m);
    marginNodes.clear();
    for (int i = numEdgeNodes(); i < static_cast<int>(nodes.size()); ++i) {
        const Node& nd = nodes[i];
        if (nd.x < lo || hiX < nd.x || nd.y < lo || hiY < nd.y) marginNodes.push_back(i);
    }
}

void VectorImage::clampEdgeNodes() {
    // Border nodes stay on their image side; the column ends sit on the image corners.
    const double W = w, H = h;
    nodes[0].y = 0.0;
    int idx = 0, last = -1;
    for (int i = 0; i < edgeNodeCounts[0]; ++i) {
        last = idx;
        nodes[idx++].x = 0.0;
    }
    nodes[last].y = H;
    nodes[idx].y = 0.0;
    for (int i = 0; i < edgeNodeCounts[1]; ++i) {
        last = idx;
        nodes[idx++].x = W;
    }
    nodes[last].y = H;
    for (int i = 0; i < edgeNodeCounts[2]; ++i) nodes[idx++].y = 0.0;
    for (int i = 0; i < edgeNodeCounts[3]; ++i) nodes[idx + i].y = H;
}

void VectorImage::clampMarginNodes() {
    const double lo = 0.01;
    const double hiX = static_cast<double>(w) + -0.01;
    const double hiY = static_cast<double>(h) + -0.01;
    for (int i : marginNodes) {
        Node& nd = nodes[i];
        nd.x = nd.x < lo ? lo : (hiX < nd.x ? hiX : nd.x);
        nd.y = nd.y < lo ? lo : (hiY < nd.y ? hiY : nd.y);
    }
}

void VectorImage::updateArea(int r) {
    Region& reg = regions[r];
    const int n = static_cast<int>(reg.nodes.size());
    double a = 0.0;
    if (n >= 1) {
        for (int i = 0; i < n; ++i) {
            const Node& p = nodes[reg.nodes[i]];
            const Node& q = nodes[reg.nodes[(i + 1) % n]];
            a = (p.x * q.y - p.y * q.x) + a;
        }
        a = a * 0.5;
    }
    reg.area = a;
}

int VectorImage::decrementNodeCounts() {
    int n = 0;
    for (Node& nd : nodes)
        if (nd.antiInv > 1) {
            ++n;
            --nd.antiInv;
        }
    return n;
}

}  // namespace isvg
