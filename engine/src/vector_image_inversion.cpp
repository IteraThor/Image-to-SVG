// VectorImage inversion detection: contour steps that cross other steps (or fold back) get their
// nodes flagged; the smoother then pulls them back with an anti-inversion penalty.
#include <cmath>
#include <initializer_list>

#include "isvg/ann_kdtree.h"
#include "isvg/vector_image.h"

namespace isvg {

namespace {

// VectorImage::isPointOnLineSegment (dir = unit direction of a→b).
bool isPointOnLineSegment(const double p[2], const double a[2], const double b[2], const double dir[2]) {
    const double dx = p[0] - a[0], dy = p[1] - a[1];
    if (!(std::fabs(dir[0] * dy - dir[1] * dx) < 0.0001)) return false;
    const double proj = dir[1] * dy + dir[0] * dx + 0.0;
    if (!(proj < 0.0)) {
        const double ex = b[0] - a[0], ey = b[1] - a[1];
        if (!(std::sqrt(ey * ey + ex * ex + 0.0) < proj)) return true;
    }
    if (!(1e-08 <= dy * dy + dx * dx + 0.0)) return true;
    const double fx = p[0] - b[0], fy = p[1] - b[1];
    return fy * fy + fx * fx + 0.0 < 1e-08;
}

//: ±1 if the segments touch at an endpoint, ±2 if they cross properly, 0 otherwise;
// the sign is that of the cross product of their directions.
int intersectionCrossProductSign(const double p1[2], const double p2[2], const double q1[2], const double q2[2]) {
    const double qdx = q2[0] - q1[0], qdy = q2[1] - q1[1];
    const double pdx = p2[0] - p1[0], pdy = p2[1] - p1[1];
    double cross = qdy * pdx - qdx * pdy;
    int sign;
    if (std::isnan(cross)) {
        sign = 1;
    } else {
        if (cross == 0.0) return 0;
        sign = cross < 0.0 ? -1 : 1;
    }
    double pu[2] = {pdx, pdy}, qu[2] = {qdx, qdy};
    const double pl = std::sqrt(pdy * pdy + pdx * pdx + 0.0);
    if (0.0 < pl) {
        pu[0] = pdx / pl;
        pu[1] = pdy / pl;
    }
    const double ql = std::sqrt(qdx * qdx + 0.0 + qdy * qdy);
    if (0.0 < ql) {
        qu[0] = qdx / ql;
        qu[1] = qdy / ql;
    }
    if (isPointOnLineSegment(p1, q1, q2, qu)) return sign;
    if (isPointOnLineSegment(p2, q1, q2, qu)) return sign;
    if (isPointOnLineSegment(q1, p1, p2, pu)) return sign;
    if (isPointOnLineSegment(q2, p1, p2, pu)) return sign;
    const double sx = p1[0] - q1[0], sy = p1[1] - q1[1];
    double u = qdx * sy - qdy * sx;
    double v = pdx * sy - sx * pdy;
    if (cross < 0.0) {
        cross = -cross;
        u = -u;
        v = -v;
    }
    if (u > 0.0 && u < cross && v > 0.0 && v < cross) return sign * 2;
    return 0;
}

//: the same test on the directions from `c`, normalised to the unit circle.
int intersectionCrossProductSignOnUnitCircle(double p1[2], double p2[2], double q1[2], double q2[2],
                                              const double c[2]) {
    for (double* v : {p1, p2, q1, q2}) {
        const double dx = v[0] - c[0], dy = v[1] - c[1];
        const double l = std::sqrt(dy * dy + dx * dx + 0.0);
        if (0.0 < l) {
            v[0] = dx / l;
            v[1] = dy / l;
        } else {
            v[0] = dx;
            v[1] = dy;
        }
    }
    return intersectionCrossProductSign(p1, p2, q1, q2);
}

//: ray p + t·r against the segment q1q2.
int rayIntersectionCrossProductSign(const double p[2], const double r[2], const double q1[2], const double q2[2]) {
    const double dy = q2[1] - q1[1], dx = q2[0] - q1[0];
    const double c = dy * r[0] - dx * r[1];
    if (c == 0.0) return 0;
    if (std::fabs((q2[1] - p[1]) * r[0] - (q2[0] - p[0]) * r[1]) < 0.0001 ||
        std::fabs((q1[1] - p[1]) * r[0] - (q1[0] - p[0]) * r[1]) < 0.0001)
        return 0.0 < c ? 1 : -1;
    const double u = dx * (p[1] - q1[1]) - dy * (p[0] - q1[0]);
    const double v = r[0] * (p[1] - q1[1]) - r[1] * (p[0] - q1[0]);
    if (c <= 0.0) {
        if (u <= 0.0 && v < 0.0) return c < v ? -2 : 0;
    } else if (0.0 <= u && 0.0 < v) {
        return v < c ? 2 : 0;
    }
    return 0;
}

}  // namespace

void VectorImage::computeInsideNormal(int r, int k, double out[2]) const {
    const Region& reg = regions[r];
    const int n = static_cast<int>(reg.nodes.size());
    const Node& p = nodes[reg.nodes[k]];
    const Node& pa = nodes[reg.nodes[(n - 1 + k) % n]];
    double ux = pa.x - p.x, uy = pa.y - p.y;
    const double ul = std::sqrt(uy * uy + ux * ux + 0.0);
    if (0.0 < ul) {
        ux = ux / ul;
        uy = uy / ul;
    }
    const Node& pb = nodes[reg.nodes[(k + 1) % n]];
    double vx = pb.x - p.x, vy = pb.y - p.y;
    const double vl = std::sqrt(vy * vy + vx * vx + 0.0);
    if (0.0 < vl) {
        vx = vx / vl;
        vy = vy / vl;
    }
    double sx = vx + ux, sy = vy + uy;
    const double d3 = uy - vy, d2 = vx - ux;
    const double m1 = d2 * d2 + d3 * d3 + 0.0;
    double m2 = sy * sy + sx * sx + 0.0;
    if (m1 < m2) {
        if (0.0 < vy * ux - uy * vx) {
            sx = -sx;
            sy = -sy;
        }
    } else {
        sx = d3;
        sy = d2;
        m2 = m1;
    }
    out[0] = sx;
    out[1] = sy;
    const double l = std::sqrt(m2);
    if (0.0 < l) {
        out[0] = sx / l;
        out[1] = sy / l;
    }
}

int VectorImage::computeBaseSignOffset(int r, int start) {
    double nrm[2];
    computeInsideNormal(r, start, nrm);
    Region& reg = regions[r];
    const int n = static_cast<int>(reg.nodes.size());
    const Node& p = nodes[reg.nodes[start]];
    reg.insideNormal[0] = nrm[0];
    reg.insideNormal[1] = nrm[1];
    if (n < 3) return 2;
    const double pp[2] = {p.x, p.y};
    unsigned sum = 0;
    int i = start + 1;
    int i1;
    do {
        i1 = i + 1;
        const Node& a = nodes[reg.nodes[i % n]];
        const Node& b = nodes[reg.nodes[i1 % n]];
        const double qa[2] = {a.x, a.y}, qb[2] = {b.x, b.y};
        sum += static_cast<unsigned>(rayIntersectionCrossProductSign(pp, nrm, qa, qb));
        i = i1;
    } while (i1 != n - 1 + start);
    if ((sum & 1) == 0) return static_cast<int>(2 - sum);
    return 2;  // odd parity: invalid intersection sign sum
}

int VectorImage::flagInversionNodesForRegion(int r) {
    Region& reg = regions[r];
    const int n = static_cast<int>(reg.nodes.size());
    reg.baseSignOffset = computeBaseSignOffset(r, 0);
    if (n < 1) return 0;
    for (int i = 0; i < n; ++i) nodes[reg.nodes[i]].flags &= 0xf0;
    int run = 0, acc = 0, count = 0;
    for (int i = 0; i < n; ++i) {
        const int ia = reg.nodes[i], ib = reg.nodes[(i + 1) % n];
        Node& a = nodes[ia];
        Node& b = nodes[ib];
        bool checkB = true;
        if (!a.iric.empty()) {
            int total = 0, signedSum = 0;
            const double pa[2] = {a.x, a.y}, pb[2] = {b.x, b.y};
            for (const auto& [r2, idx2] : a.iric) {
                const Region& reg2 = regions[r2];
                const int n2 = static_cast<int>(reg2.nodes.size());
                const int ic = reg2.nodes[idx2], id = reg2.nodes[(idx2 + 1) % n2];
                const Node& c = nodes[ic];
                const Node& d = nodes[id];
                if (r != r2 || ib != id) {
                    if (ia != ic && ia != id && ib != ic && ib != id) {
                        const double pc[2] = {c.x, c.y}, pd[2] = {d.x, d.y};
                        const int s = intersectionCrossProductSign(pa, pb, pc, pd);
                        total += s < 0 ? -s : s;
                        if (r == r2) signedSum += s;
                    }
                } else {
                    const Node& an = nodes[reg.nodes[(i + 2) % n]];
                    const Node& dn = nodes[reg2.nodes[(idx2 + 2) % n2]];
                    double v0[2] = {a.x, a.y}, v1[2] = {an.x, an.y}, v2[2] = {c.x, c.y}, v3[2] = {dn.x, dn.y};
                    const int s = intersectionCrossProductSignOnUnitCircle(v0, v1, v2, v3, pb);
                    total += s < 0 ? -s : s;
                    signedSum += s;
                }
            }
            acc += signedSum;
            if (total != 0) {
                const uint8_t fa = a.flags;
                a.antiInv = 3;
                a.flags = fa | 0x16;
                b.antiInv = 3;
                const uint8_t fb = b.flags;
                b.flags = fb | 0x16;
                count += ((fb >> 1) ^ 1) & 1;
                count += ((fa >> 1) ^ 1) & 1;
                checkB = false;
            }
        }
        if (checkB) {
            const uint8_t fb = b.flags;
            if (acc + reg.baseSignOffset == 0) {
                run = 0;
                b.flags = fb | 4;
            } else {
                ++run;
                count += ((fb >> 1) ^ 1) & 1;
                b.flags = fb | 0x16;
                b.antiInv = 3;
                if (run > 20) {
                    reg.baseSignOffset = computeBaseSignOffset(r, i);
                    run = -0x7fffffff;
                    acc = 0;
                }
            }
        }
    }
    return count;
}

int VectorImage::flagInversionNodes() {
    int count = decrementNodeCounts();
    for (int r = 0; r < static_cast<int>(regions.size()); ++r) count += flagInversionNodesForRegion(r);
    return count;
}

void VectorImage::buildIricNeighbors() {
    // Every contour element (region, index) is a kd-tree point at its node's position; each node
    // records the elements (other than itself) within distance 5 among its 50 nearest.
    std::vector<std::pair<int, int>> elems;
    std::vector<const double*> pts;
    for (int r = 0; r < static_cast<int>(regions.size()); ++r)
        for (int i = 0; i < static_cast<int>(regions[r].nodes.size()); ++i) {
            elems.emplace_back(r, i);
            pts.push_back(&nodes[regions[r].nodes[i]].x);
        }
    ann::KdTree tree(pts);
    const int numNodes = static_cast<int>(nodes.size());
    const int k = numNodes < 0x33 ? numNodes : 0x32;
    std::vector<int> idx(k > 0 ? k : 1);
    std::vector<double> dist(k > 0 ? k : 1);
    for (int r = 0; r < static_cast<int>(regions.size()); ++r) {
        const Region& reg = regions[r];
        for (int i = 0; i < static_cast<int>(reg.nodes.size()); ++i) {
            const int ni = reg.nodes[i];
            Node& nd = nodes[ni];
            if (nd.flags & 1) continue;
            tree.kSearch(&nd.x, k, idx.data(), dist.data(), 1e-12);
            nd.iric.clear();
            for (int j = 0; j < k; ++j) {
                if (!(dist[j] < 25.0)) continue;
                const auto& e = elems[idx[j]];
                if (regions[e.first].nodes[e.second] != ni) nd.iric.push_back(e);
            }
            nd.flags |= 1;
        }
    }
}

}  // namespace isvg
