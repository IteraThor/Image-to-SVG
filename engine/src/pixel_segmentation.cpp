#include "isvg/pixel_segmentation.h"

#include <list>

namespace isvg {

void PixelSegmentation::initialize() {
    for (int i = 0; i < image_.numPixels(); ++i) map_.raw(i) = i;
    map_.resetRedirects();
    initialized_ = false;
    segments_.clear();
}

// Relabel segments by 4-connected id components in scan order of their first pixel; disconnected
// pieces of one id become separate segments. Pending redirects are resolved first.
void PixelSegmentation::renumber() {
    bool carry = initialized_ && hasSegmentRecords() && paletteIndexOf(0) >= 0;
    std::vector<int> palette;
    if (carry) {
        palette.resize(map_.maxId() + 1);
        for (int i = 0; i <= map_.maxId(); ++i) palette[i] = paletteIndexOf(i);
    }
    map_.flushRedirectsIfDirty();
    std::list<int> carried;
    int k = 0;
    const int n = image_.numPixels();
    for (int idx = 0; idx < n; ++idx) {
        const int v = map_.raw(idx);
        if (v >= 0) {
            if (carry) carried.push_back(palette[v]);
            map_.floodFill(v, ~k, idx);
            ++k;
        }
    }
    updateMaxId(k - 1);
    if (carry) {
        auto it = carried.begin();
        for (int i = 0; i < k; ++i, ++it) paletteIndexOf(i) = *it;
    }
    for (int idx = 0; idx < n; ++idx) map_.raw(idx) = ~map_.raw(idx);
}

void PixelSegmentation::updateMaxId(int maxId) {
    map_.setMaxId(maxId);
    const int count = static_cast<int>(segments_.size());
    if (count <= maxId) {
        const int want = maxId + 1;
        if (initialized_ && hasSegmentRecords() && paletteIndexOf(0) >= 0) {
            segments_.resize(want);  // keep existing segments (Array1::resize(n, copy=true))
        } else if (want != count) {
            segments_.assign(want, PixelSegment{});
        }
    }
    initialized_ = true;
}

void PixelSegmentation::updateStatistics() {
    renumber();
    for (int i = 0; i <= map_.maxId(); ++i) segments_[i].clearStats();
    const int n = image_.numPixels();
    for (int idx = 0; idx < n; ++idx) {
        PixelSegment& s = segments_[map_.idOf(idx)];
        const uint8_t* p = image_.at(idx);
        const float c[4] = {p[0] * kInv255f, p[1] * kInv255f, p[2] * kInv255f, p[3] * kInv255f};
        s.addPixelF(c);
    }
}

namespace {

// PixelMap::wouldChangeTopology
bool wouldChangeTopology(PixelMap& map, int x, int y, int dx, int dy) {
    const int a = map.idOf(x, y);
    const int b = map.idOf(x, y + dy);
    const int c = map.idOf(x + dx, y);
    const int d = map.idOf(x + dx, y + dy);
    return a == b && a == c && a != d;
}

// PixelSegmentation::computeColorDistance: Σ((p(x1,y1) − p(x2,y2)) / 255)²
float colorDistance(const Image& img, int x1, int y1, int x2, int y2) {
    const uint8_t* p2 = img.at(y2 * img.w + x2);
    const uint8_t* p1 = img.at(y1 * img.w + x1);
    float s = 0.0f;
    float d[4];
    for (int k = 0; k < 4; ++k) {
        const float t = (static_cast<float>(p1[k]) - static_cast<float>(p2[k])) / 255.0f;
        d[k] = t * t;
    }
    s = d[0] + s;
    s = s + d[1];
    s = d[2] + s;
    s = s + d[3];
    return s;
}

// PixelSegmentation::computeDiagonalCurvature: Σ(((p1 − 2·p2) + p3) / 255)², 0 if an
// endpoint is outside the image.
float diagonalCurvature(const Image& img, int x1, int y1, int x2, int y2, int x3, int y3) {
    if (!(x1 >= 0 && x1 < img.w && y1 >= 0 && y1 < img.h && x3 >= 0 && x3 < img.w && y3 >= 0 && y3 < img.h))
        return 0.0f;
    const uint8_t* p3 = img.at(y3 * img.w + x3);
    const uint8_t* p2 = img.at(y2 * img.w + x2);
    const uint8_t* p1 = img.at(y1 * img.w + x1);
    float d[4];
    for (int k = 0; k < 4; ++k) {
        const float m = static_cast<float>(p2[k]);
        const float t = ((static_cast<float>(p1[k]) - (m + m)) + static_cast<float>(p3[k])) * kInv255f;
        d[k] = t * t;
    }
    float s = 0.0f;
    s = d[0] + s;
    s = s + d[1];
    s = d[2] + s;
    s = s + d[3];
    return s;
}

}  // namespace

void PixelSegmentation::performDiagonalExtrusions() {
    const int w = image_.w, h = image_.h;
    if (h <= 1) return;
    for (int y = 0; y + 1 < h; ++y) {
        if (w <= 1) continue;
        for (int x = 0; x + 1 < w; ++x) {
            const int a = map_.idOf(x, y), b = map_.idOf(x + 1, y);
            const int c = map_.idOf(x, y + 1), d = map_.idOf(x + 1, y + 1);
            if (a == c || a == b || c == d || b == d) continue;
            const bool ad = a == d, bc = b == c;
            if (!ad && !bc) continue;
            bool mainDiagonal;  // extend a along TL–BR (true) or c along TR–BL (false)
            if (ad && bc) {
                float sMain = 0.0f, sAnti = 0.0f;
                for (int j = y - 1; j <= y; ++j) {
                    for (int i = x - 1; i <= x; ++i) {
                        sMain = sMain + diagonalCurvature(image_, i, j, i + 1, j + 1, i + 2, j + 2);
                        sAnti = sAnti + diagonalCurvature(image_, i, j + 2, i + 1, j + 1, i + 2, j);
                    }
                }
                if (sMain > sAnti + 0.0f) {
                    mainDiagonal = false;
                } else if (sAnti > sMain + 0.0f) {
                    mainDiagonal = true;
                } else {
                    continue;
                }
            } else {
                mainDiagonal = ad;
            }
            // Candidate pixels P (first) and Q (second) that could receive the diagonal's id.
            int px, py, qx, qy, id, t1dx, t1dy, t2dx, t2dy;
            if (mainDiagonal) {
                id = a;
                px = x + 1; py = y; qx = x; qy = y + 1;
                t1dx = 1; t1dy = -1; t2dx = -1; t2dy = 1;
            } else {
                id = c;
                px = x; py = y; qx = x + 1; qy = y + 1;
                t1dx = -1; t1dy = -1; t2dx = 1; t2dy = 1;
            }
            const bool t1 = wouldChangeTopology(map_, px, py, t1dx, t1dy);
            const bool t2 = wouldChangeTopology(map_, qx, qy, t2dx, t2dy);
            bool fillP;
            if (t1 && !t2) {
                fillP = false;
            } else if (!t1 && t2) {
                fillP = true;
            } else if (t1 && t2) {
                continue;
            } else {
                // The same four pixel pairs are compared.
                float dP, dQ;
                if (mainDiagonal) {
                    dP = colorDistance(image_, x, y, x + 1, y);
                    dP = dP + colorDistance(image_, x + 1, y + 1, x + 1, y);
                    float d3 = colorDistance(image_, x, y, x, y + 1);
                    float d4 = colorDistance(image_, x + 1, y + 1, x, y + 1);
                    dQ = d4 + d3;
                } else {
                    dP = colorDistance(image_, x, y + 1, x, y);
                    dP = dP + colorDistance(image_, x + 1, y, x, y);
                    float d3 = colorDistance(image_, x, y + 1, x + 1, y + 1);
                    float d4 = colorDistance(image_, x + 1, y, x + 1, y + 1);
                    dQ = d4 + d3;
                }
                fillP = dQ > dP;
            }
            if (fillP)
                map_.raw(py * w + px) = id;
            else
                map_.raw(qy * w + qx) = id;
        }
    }
}

void PixelSegmentation::getColorF(int id, float out[4]) {
    const int pi = paletteIndexOf(id);
    if (pi < 0) {
        ownColorF(id, out);
        return;
    }
    const float* c = palette_->color(pi);
    for (int k = 0; k < 4; ++k) out[k] = c[k];
}

void PixelSegmentation::generatePaletteFromSegments() {
    if (initialized_ && hasSegmentRecords() && paletteIndexOf(0) >= 0) return;  // "already using palette"
    palette_->clearMap();
    for (int i = 0; i <= map_.maxId(); ++i) {
        float c[4];
        getColorF(i, c);
        uint8_t b[4];
        for (int k = 0; k < 4; ++k) b[k] = static_cast<uint8_t>(static_cast<int>(c[k] * 255.0f + 0.5f));
        const int idx = palette_->tryAddColorToMap(b);
        paletteIndexOf(i) = idx;
    }
    palette_->setNewPaletteFromMap();
}

void PixelSegmentation::compact() {
    segments_.resize(map_.maxId() + 1);
    map_.compact();
}

}  // namespace isvg
