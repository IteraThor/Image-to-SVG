#include "isvg/shadows.h"

#include "isvg/math_util.h"

namespace isvg {
namespace {

inline float seqSum4(const float v[4]) {
    float s = 0.0f;
    s = v[0] + s;
    s = s + v[1];
    s = s + v[2];
    s = s + v[3];
    return s;
}

inline void loadF(const Image& img, int x, int y, float out[4]) {
    const uint8_t* p = img.at(y * img.w + x);
    for (int k = 0; k < 4; ++k) out[k] = static_cast<float>(p[k]) * kInv255f;
}

inline float sqDiff(const float a[4], const float b[4]) {
    float d[4];
    for (int k = 0; k < 4; ++k) {
        const float t = a[k] - b[k];
        d[k] = t * t;
    }
    return seqSum4(d);
}

}  // namespace

Image blurPixels(const Image& src) {
    double K[25];
    {
        int i = 0;
        for (int dy = -2; dy <= 2; ++dy) {
            const double ry = static_cast<double>(dy * dy) / 1.125;
            for (int dx = -2; dx <= 2; ++dx) {
                const double e = static_cast<double>(dx * dx) / 1.125 + ry;
                K[i++] = math::exp(-e);
            }
        }
        double sum = 0.0 + K[0];
        for (int k = 1; k < 25; ++k) sum = sum + K[k];
        const double norm = 1.0 / sum;
        for (double& v : K) v = v * norm;
    }
    Image out;
    out.w = src.w;
    out.h = src.h;
    out.px.resize(src.px.size());
    for (int y = 0; y < src.h; ++y) {
        for (int x = 0; x < src.w; ++x) {
            float acc[4] = {0, 0, 0, 0};
            double wsum = 0.0;
            int k = 0;
            for (int j = y - 2; j <= y + 2; ++j) {
                for (int i = x - 2; i <= x + 2; ++i, ++k) {
                    if (i < 0 || j < 0 || i >= src.w || j >= src.h) continue;
                    wsum = wsum + K[k];
                    const float w = static_cast<float>(K[k]);
                    float c[4];
                    loadF(src, i, j, c);
                    for (int ch = 0; ch < 4; ++ch) acc[ch] = w * c[ch] + acc[ch];
                }
            }
            const float scale = static_cast<float>(1.0 / wsum);
            uint8_t* o = &out.px[(static_cast<size_t>(y) * src.w + x) * 4];
            for (int ch = 0; ch < 4; ++ch) {
                const float v = scale * acc[ch];
                o[ch] = static_cast<uint8_t>(static_cast<int>(v * 255.0f + 0.5f));
            }
        }
    }
    return out;
}

SignalImages computeSignals(const Image& src) {
    const Image blur = blurPixels(src);
    SignalImages s;
    s.w = src.w;
    s.h = src.h;
    s.s1.assign(static_cast<size_t>(s.w) * s.h, 0.0f);
    s.s2.assign(static_cast<size_t>(s.w) * s.h, 0.0f);
    const int w = s.w, h = s.h;
    for (int y = 0; y < h; ++y) {
        const int down = y < h - 1 ? y + 1 : y - 1;
        const int up = y != 0 ? y - 1 : y + 1;
        for (int x = 0; x < w; ++x) {
            int right, left;
            if (x == 0) {
                right = w > 1 ? 1 : -1;
                left = 1;
            } else {
                right = x < w - 1 ? x + 1 : x - 1;
                left = x - 1;
            }
            float c[4], n[4];
            loadF(src, x, y, c);
            float& s1 = s.s1[static_cast<size_t>(y) * w + x];
            loadF(src, right, y, n);
            s1 = sqDiff(c, n);
            loadF(src, left, y, n);
            s1 = sqDiff(c, n) + s1;
            loadF(src, x, down, n);
            s1 = s1 + sqDiff(c, n);
            loadF(src, x, up, n);
            s1 = sqDiff(c, n) + s1;

            float bu[4], bd[4], bl[4], br[4], bc[4], avg[4];
            loadF(blur, x, up, bu);
            loadF(blur, x, down, bd);
            loadF(blur, left, y, bl);
            loadF(blur, right, y, br);
            for (int k = 0; k < 4; ++k) avg[k] = (((br[k] + bl[k]) + bd[k]) + bu[k]) * 0.25f;
            loadF(blur, x, y, bc);
            s.s2[static_cast<size_t>(y) * w + x] = sqDiff(bc, avg);
        }
    }
    return s;
}

std::vector<uint8_t> boundaryMask8(const PixelMap& map, int level) {
    const int w = map.width(), h = map.height();
    std::vector<uint8_t> m(static_cast<size_t>(w) * h, 0);
    const auto& ids = map.rawIds();
    for (int y = 0; y + 1 < h; ++y) {
        for (int x = 0; x + 1 < w; ++x) {
            const int a = ids[y * w + x];
            if (a == ids[y * w + x + 1] && a == ids[(y + 1) * w + x] && a == ids[(y + 1) * w + x + 1]) continue;
            m[y * w + x] = m[y * w + x + 1] = m[(y + 1) * w + x] = m[(y + 1) * w + x + 1] =
                static_cast<uint8_t>(level);
        }
    }
    for (int k = level - 1; k > 0; --k) {
        for (int y = 0; y + 1 < h; ++y) {
            for (int x = 0; x + 1 < w; ++x) {
                const int cells[4] = {y * w + x, y * w + x + 1, (y + 1) * w + x, (y + 1) * w + x + 1};
                bool any = false;
                for (int c : cells)
                    if (static_cast<int8_t>(m[c]) > k) any = true;
                if (!any) continue;
                for (int c : cells)
                    if (m[c] == 0) m[c] = static_cast<uint8_t>(k);
            }
        }
    }
    return m;
}

bool isProbablyShadowingSegment(const FeatureVector& f) {
    // Each comparison mirrors COMISD + Jcc: a jump taken on "below or equal" is also taken for NaN.
    if (!(0.00402429 > f[23])) {
        if (!(0.0241458 > f[23])) return false;
        if (!(0.0473801 > f[5])) return false;
        if (!(0.0138127 > f[22])) return true;
        return !(0.0102382 > f[39]);
    }
    if (!(1.5 > f[18])) {
        if (0.5 > f[20]) return false;
        if (0.5 > f[19]) return 0.0112894 > f[1];
        if (!(0.449972 > f[41])) return false;
        if (0.118876 > f[21]) return 41.535 > f[0];
        return 5.5 > f[11];
    }
    if (1.5 > f[14]) return false;
    if (!(0.0491163 > f[5])) return !(0.196468 > f[22]);
    return 1.5 > f[20];
}

bool isProbablyNotFalsePositive(const FeatureVector& f, int configT) {
    if (!(3.57554 > f[37])) return 0.00258463 > f[35];
    if (!(0.0442212 > f[28])) {
        if (configT == 0) return true;
        if (!(0.694472 > f[34])) return false;
        return !(3.09081 > f[0]);
    }
    if (!(0.00627504 > f[39])) return true;
    if (static_cast<unsigned>(configT - 1) <= 1u) return true;
    return !(0.934367 > f[4]);
}

}  // namespace isvg
