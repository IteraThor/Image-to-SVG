#include "isvg/palette_finder.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <utility>

#include "isvg/prng.h"
#include "isvg/shadows.h"
#include "isvg/weights_finder.h"

namespace isvg {

Image downsamplePixels(const Image& src) {
    const int W = src.w, H = src.h;
    if (W * H < 200001) return src;
    const double s = std::sqrt(150000.0 / static_cast<double>(W * H));
    Image out;
    out.w = static_cast<int>(static_cast<double>(W) * s + 0.5);
    out.h = static_cast<int>(s * static_cast<double>(H) + 0.5);
    out.px.resize(static_cast<size_t>(out.w) * out.h * 4);
    // downsampleAliased: nearest neighbour at the cell centres.
    const double rx = static_cast<double>(W) / static_cast<double>(out.w);
    const double ry = static_cast<double>(H) / static_cast<double>(out.h);
    for (int y = 0; y < out.h; ++y) {
        int sy = static_cast<int>((static_cast<double>(y) + 0.5) * ry);
        if (H < sy) sy = H;
        for (int x = 0; x < out.w; ++x) {
            int sx = static_cast<int>((static_cast<double>(x) + 0.5) * rx);
            if (W < sx) sx = W;
            std::memcpy(&out.px[(static_cast<size_t>(y) * out.w + x) * 4], src.at(sy * W + sx), 4);
        }
    }
    return out;
}

namespace {

inline float seqSum4(const float v[4]) {
    float s = 0.0f;
    for (int k = 0; k < 4; ++k) s = v[k] + s;
    return s;
}

inline void loadF(const Image& img, int i, float c[4]) {
    const uint8_t* p = img.at(i);
    for (int k = 0; k < 4; ++k) c[k] = static_cast<float>(p[k]) * kInv255f;
}

void randomizeColor(std::array<float, 4>& c) {
    c[2] = static_cast<float>(prngNext());
    c[1] = static_cast<float>(prngNext());
    c[3] = 1.0f;
    c[0] = static_cast<float>(prngNext());
}

}  // namespace

void PaletteFinder::doHistogram() {
    ColorHistogram hist;
    hist.setNumBins16();  // use_16_hist_bins
    const int n = img_.numPixels();
    for (int i = 0; i < n; ++i) hist.addColor(img_.at(i), aa_ ? edgeWeightLinear(i) : 1.0);
    hist.findPeaks(peaks_, peakWeights_);
}

void PaletteFinder::initializeK(int k) {
    cur_.colors.assign(k, {0, 0, 0, 0});  // Palette::setNumColors
    means_.assign(k, {0, 0, 0, 0});
    meanWeights_.assign(k, 0.0f);
    map_.setMaxId(k);
    for (int i = 0; i < k; ++i) {
        if (i < static_cast<int>(peaks_.size())) {
            for (int c = 0; c < 4; ++c) cur_.colors[i][c] = static_cast<float>(peaks_[i][c]) * kInv255f;
        } else {
            randomizeColor(cur_.colors[i]);
        }
    }
}

void PaletteFinder::assignK(bool final) {
    const int k = static_cast<int>(cur_.colors.size());
    if (final) cur_.weights.assign(k, 0.0);
    const int n = img_.numPixels();
    for (int i = 0; i < n; ++i) {
        float c[4];
        loadF(img_, i, c);
        auto dist = [&](int p) {
            float sq[4];
            for (int t = 0; t < 4; ++t) {
                const float d = c[t] - cur_.colors[p][t];
                sq[t] = d * d;
            }
            return seqSum4(sq);
        };
        float best = dist(0);
        int idx = 0;
        for (int p = 1; p < k; ++p) {
            const float d = dist(p);
            if (best > d) idx = p;
            best = d < best ? d : best;  // MINSS d, best
        }
        map_.raw(i) = idx;
        if (final) cur_.weights[idx] = (aa_ ? edgeIntensity(i) : 1.0) + cur_.weights[idx];
    }
}

void PaletteFinder::computeMeanColors(bool weighted) {
    const int k = static_cast<int>(means_.size());
    std::fill(meanWeights_.begin(), meanWeights_.end(), 0.0f);
    for (auto& m : means_) m = {0, 0, 0, 0};
    const int w = img_.w, h = img_.h;
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            const int i = y * w + x;
            const int id = map_.raw(i);
            float c[4];
            if (weighted) {
                const float wt = static_cast<float>(edgeWeightXY(x, y));
                loadF(img_, i, c);
                for (int t = 0; t < 4; ++t) means_[id][t] = means_[id][t] + c[t] * wt;
                meanWeights_[id] = wt + meanWeights_[id];
            } else {
                loadF(img_, i, c);
                for (int t = 0; t < 4; ++t) means_[id][t] = means_[id][t] + c[t];
                meanWeights_[id] = meanWeights_[id] + 1.0f;
            }
        }
    }
    for (int i = 0; i < k; ++i) {
        if (meanWeights_[i] > 0.0f) {
            for (int t = 0; t < 4; ++t) means_[i][t] = means_[i][t] / meanWeights_[i];
        } else {
            randomizeColor(means_[i]);
        }
    }
}

void PaletteFinder::doKMeans(int k) {
    initializeK(k);
    bool converged = false;
    for (int it = 0; it < 10; ++it) {
        assignK(it == 9 || converged);
        if (converged) break;
        computeMeanColors(aa_);
        float maxDiff = 0.0f;
        for (int i = 0; i < k; ++i) {
            float d[4];
            for (int t = 0; t < 4; ++t) d[t] = std::fabs(cur_.colors[i][t] - means_[i][t]);
            cur_.colors[i] = means_[i];
            float m = d[1] > d[0] ? d[1] : d[0];  // MAXSS chain
            m = d[2] > m ? d[2] : m;
            m = d[3] > m ? d[3] : m;
            maxDiff = m > maxDiff ? m : maxDiff;
        }
        converged = 0.01f > maxDiff;
    }
    for (int i = 0; i < k; ++i)
        if (meanWeights_[i] == 0.0f) cur_.colors[i] = {0, 0, 0, 0};
}

float PaletteFinder::computePixelCost(int x, int y, const std::vector<uint8_t>& mask, WeightsFinder& wf) {
    const int i = y * img_.w + x;
    float model[4];
    if (aa_ && static_cast<int8_t>(mask[i]) >= 1) {
        const Color4f m = wf.computePixelF(x, y, -1, -1, nullptr, false);
        std::memcpy(model, m.c, sizeof model);
    } else {
        std::memcpy(model, cur_.colors[map_.raw(i)].data(), sizeof model);
    }
    float c[4], sq[4];
    loadF(img_, i, c);
    for (int t = 0; t < 4; ++t) {
        const float d = model[t] - c[t];
        sq[t] = d * d;
    }
    return seqSum4(sq);
}

void PaletteFinder::computeCost() {
    const std::vector<uint8_t> mask = boundaryMask8(map_, p_.maskLevel);
    // The engine's WeightsFinder  in palette-colour mode.
    WeightsFinder wf;
    wf.clear();
    wf.attach(map_, img_, [this](int id) { return cur_.colors[id].data(); });
    float cost = 0.0f;
    for (int y = 0; y < img_.h; ++y)
        for (int x = 0; x < img_.w; ++x) cost = cost + computePixelCost(x, y, mask, wf);
    const float logCost = static_cast<float>(std::log10(static_cast<double>(cost) + 1.0));
    const float colorCost = static_cast<float>(static_cast<int>(cur_.colors.size()) - 2) * p_.costPerColor;
    cur_.logCost = logCost;
    cur_.colorCost = colorCost;
    cur_.cost = logCost + colorCost;
}

void PaletteFinder::findPalette(int k) {
    doKMeans(k);
    computeCost();
    for (auto& c : cur_.colors) {
        uint8_t b[4];
        for (int t = 0; t < 4; ++t) b[t] = static_cast<uint8_t>(static_cast<int>(c[t] * 255.0f + 0.5f));
        for (int t = 0; t < 4; ++t) c[t] = static_cast<float>(b[t]) * kInv255f;
    }
    palettes_[k] = cur_;
}

void PaletteFinder::computeFeatures() {
    double* f = features_;
    {
        std::vector<uint32_t> px(static_cast<size_t>(img_.numPixels()));
        std::memcpy(px.data(), img_.px.data(), px.size() * 4);
        std::sort(px.begin(), px.end());
        f[0] = static_cast<double>(std::unique(px.begin(), px.end()) - px.begin());
    }
    f[1] = static_cast<double>(img_.numPixels());
    using Entry = std::pair<double, int>;
    std::vector<Entry> a(13, {1000000.0, -1}), b(13, {1000000.0, -1});
    double norm[13];
    for (double& v : norm) v = 1.0;
    double mn = static_cast<double>(palettes_[2].logCost), mx = mn;
    for (int k = 2; k <= 12; ++k) {
        a[k] = {static_cast<double>(palettes_[k].cost), k};
        const double v = static_cast<double>(palettes_[k].logCost);
        norm[k] = v;
        b[k] = {v, k};
        mn = mn <= v ? mn : v;
        mx = v <= mx ? mx : v;
    }
    for (int k = 2; k <= 12; ++k) norm[k] = (norm[k] - mn) / (mx - mn);
    // sla::compareAscending<pair<double,int>>: by value, then by k (also when unordered).
    auto less = [](const Entry& x, const Entry& y) {
        if (x.first < y.first) return true;
        if (x.first > y.first) return false;
        return x.second < y.second;
    };
    std::sort(a.begin(), a.end(), less);
    std::sort(b.begin(), b.end(), less);
    int rankA[13], rankB[13];
    for (int i = 0; i < 13; ++i) rankA[i] = rankB[i] = -1;
    for (int i = 0; i < 13; ++i) {
        if (a[i].second > 1) rankA[a[i].second] = i;
        if (b[i].second > 1) rankB[b[i].second] = i;
    }
    for (int i = 0; i < 4; ++i) {
        f[2 + i] = static_cast<double>(a[i].second);
        f[6 + i] = static_cast<double>(b[i].second);
        f[10 + i] = b[i].first;
    }
    for (int k = 2; k <= 12; ++k) {
        f[14 + k - 2] = static_cast<double>(rankA[k]);
        f[25 + k - 2] = static_cast<double>(rankB[k]);
    }
    f[36] = static_cast<double>(peaks_.size());
    double pm = peakWeights_.empty() ? 0.0 : peakWeights_[0];
    for (double v : peakWeights_) pm = pm <= v ? pm : v;
    f[37] = pm;
    for (int k = 2; k <= 12; ++k) f[38 + k - 2] = norm[k];
    f[49] = mn;
    f[50] = mx;
    for (int k = 2; k <= 12; ++k) {
        const auto& w = palettes_[k].weights;
        double m = w[0];
        for (double v : w) m = m <= v ? m : v;
        f[51 + k - 2] = m;
    }
}

void PaletteFinder::findPalettesRange(int kMin, int kMax) {
    doHistogram();
    for (int k = kMin; k <= kMax; ++k) findPalette(k);
}

std::string PaletteFinder::findPalettes(const std::string& ext, const std::string& typeString) {
    (void)ext;
    doHistogram();
    for (int k = 2; k <= 12; ++k) findPalette(k);
    computeFeatures();
    return classifyPaletteMode(features_, typeString);
}

}  // namespace isvg
