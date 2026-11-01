#include "isvg/parameter_classifier.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iterator>
#include <list>

#include "isvg/palette_finder.h"
#include "isvg/pixel_segmenter.h"
#include "isvg/shadows.h"

namespace isvg {

// ---------------------------------------------------------------------------------------------
// ColorHistogram (16 bins per channel, shift 4)

void ColorHistogram::setNumBins16() {
    bins_.assign(static_cast<size_t>(kBins) * kBins * kBins * kBins, 0.0);
}

void ColorHistogram::addColor(const uint8_t c[4], double w) {
    double& b = bins_[index(c[0] >> kShift, c[1] >> kShift, c[2] >> kShift, c[3] >> kShift)];
    b = w + b;
}

double ColorHistogram::totalWeight() const {
    double s = 0.0;
    for (double v : bins_) s = s + v;
    return s;
}

void ColorHistogram::zeroOutColor(const uint8_t c[4], bool neighborhood) {
    const int i0 = c[0] >> kShift, i1 = c[1] >> kShift, i2 = c[2] >> kShift, i3 = c[3] >> kShift;
    if (!neighborhood) {
        bins_[index(i0, i1, i2, i3)] = 0.0;
        return;
    }
    for (int a = i0 - 1; a <= i0 + 1; ++a)
        for (int b = i1 - 1; b <= i1 + 1; ++b)
            for (int d = i2 - 1; d <= i2 + 1; ++d)
                for (int e = i3 - 1; e <= i3 + 1; ++e)
                    if (a >= 0 && a < kBins && b >= 0 && b < kBins && d >= 0 && d < kBins && e >= 0 && e < kBins)
                        bins_[index(a, b, d, e)] = 0.0;
}

namespace {

struct PeakNode {
    double weight;
    uint8_t color[4];
};

// The peak list's ordering (helper's comparator): heavier first; on equal weight, `a` goes
// first if any byte of `b` is smaller than the same byte of `a` (not a strict weak order, so the
// libstdc++ list::sort algorithm below has to be reproduced exactly).
bool peakBefore(const PeakNode& a, const PeakNode& b) {
    if (b.weight < a.weight) return true;
    if (!(b.weight <= a.weight)) return false;
    for (int k = 0; k < 4; ++k)
        if (b.color[k] < a.color[k]) return true;
    return false;
}

// libstdc++ std::list::merge.
void listMerge(std::list<PeakNode>& self, std::list<PeakNode>& x) {
    auto first1 = self.begin();
    auto first2 = x.begin();
    while (first1 != self.end() && first2 != x.end()) {
        if (peakBefore(*first2, *first1)) {
            auto next = std::next(first2);
            self.splice(first1, x, first2);
            first2 = next;
        } else {
            ++first1;
        }
    }
    if (first2 != x.end()) self.splice(self.end(), x, first2, x.end());
}

// libstdc++ std::list::sort (bottom-up merge with 64 bins).
void listSort(std::list<PeakNode>& l) {
    if (l.empty() || std::next(l.begin()) == l.end()) return;
    std::list<PeakNode> carry;
    std::list<PeakNode> tmp[64];
    std::list<PeakNode>* fill = tmp;
    std::list<PeakNode>* counter;
    do {
        carry.splice(carry.begin(), l, l.begin());
        for (counter = tmp; counter != fill && !counter->empty(); ++counter) {
            listMerge(*counter, carry);
            carry.swap(*counter);
        }
        carry.swap(*counter);
        if (counter == fill) ++fill;
    } while (!l.empty());
    for (counter = tmp + 1; counter != fill; ++counter) listMerge(*counter, *(counter - 1));
    l.swap(*(fill - 1));
}

}  // namespace

void ColorHistogram::findPeaks(std::vector<std::array<uint8_t, 4>>& colors, std::vector<double>& weights) const {
    std::list<PeakNode> peaks;
    for (int i3 = 0; i3 < kBins; ++i3)
        for (int i2 = 0; i2 < kBins; ++i2)
            for (int i1 = 0; i1 < kBins; ++i1)
                for (int i0 = 0; i0 < kBins; ++i0) {
                    const double v = bins_[index(i0, i1, i2, i3)];
                    if (!(2.0 <= v)) continue;
                    bool peak = true;
                    for (int a = i3 - 1; a <= i3 + 1 && peak; ++a)
                        for (int b = i2 - 1; b <= i2 + 1 && peak; ++b)
                            for (int d = i1 - 1; d <= i1 + 1 && peak; ++d)
                                for (int e = i0 - 1; e <= i0 + 1; ++e) {
                                    if (e < 0 || e >= kBins || d < 0 || d >= kBins || b < 0 || b >= kBins ||
                                        a < 0 || a >= kBins)
                                        continue;
                                    if (v < bins_[index(e, d, b, a)]) {
                                        peak = false;
                                        break;
                                    }
                                }
                    if (!peak) continue;
                    constexpr int half = 1 << (kShift - 1);
                    PeakNode n{v, {static_cast<uint8_t>((i0 << kShift) + half), static_cast<uint8_t>((i1 << kShift) + half),
                                   static_cast<uint8_t>((i2 << kShift) + half), static_cast<uint8_t>((i3 << kShift) + half)}};
                    peaks.push_back(n);
                }
    listSort(peaks);
    colors.clear();
    weights.clear();
    for (const PeakNode& n : peaks) {
        colors.push_back({n.color[0], n.color[1], n.color[2], n.color[3]});
        weights.push_back(n.weight);
    }
}

// ---------------------------------------------------------------------------------------------
// EdgeImage

std::vector<uint8_t> edgeImage(const Image& img) {
    const int w = img.w, h = img.h;
    std::vector<uint8_t> e(static_cast<size_t>(w) * h);
    auto load = [&](int x, int y, float out[4]) {
        const uint8_t* p = img.at(y * w + x);
        for (int k = 0; k < 4; ++k) out[k] = static_cast<float>(p[k]) * kInv255f;
    };
    auto sq = [](const float a[4], const float b[4]) {
        float s = 0.0f;
        for (int k = 0; k < 4; ++k) {
            const float d = b[k] - a[k];
            s = d * d + s;
        }
        return s;
    };
    for (int y = 0; y < h; ++y) {
        const int ym = y != 0 ? y - 1 : 1;
        const int yp = y < h - 1 ? y + 1 : y - 1;
        for (int x = 0; x < w; ++x) {
            const int xm = x != 0 ? x - 1 : 1;
            const int xp = x < w - 1 ? x + 1 : x - 1;
            float c[4], n[4];
            load(x, y, c);
            load(x, ym, n);
            double s = static_cast<double>(sq(c, n)) + 0.0;
            load(x, yp, n);
            s = s + static_cast<double>(sq(c, n));
            load(xm, y, n);
            s = static_cast<double>(sq(c, n)) + s;
            load(xp, y, n);
            s = static_cast<double>(sq(c, n)) + s;
            const double v = (std::sqrt(s) * 255.0) / 3.4641016151377544;
            uint8_t b;
            if (0.0 > v) b = 0;
            else if (v > 255.0) b = 0xff;
            else b = static_cast<uint8_t>(static_cast<int>(v));
            e[static_cast<size_t>(y) * w + x] = b;
        }
    }
    return e;
}

double edgeWeight(uint8_t e, double scale, double floor) {
    const double v = 1.0 - (static_cast<double>(e) / 255.0) * scale;
    return floor > v ? floor : v;
}

// ---------------------------------------------------------------------------------------------
// Edge profile (ParameterClassifier::isStrongEdge / computeEdgeProfileFeature)

namespace {

bool isStrongEdge(const Image& img, int x, int y, bool horizontal, float threshold, float p[5]) {
    float c[6][4];
    for (int i = 0; i < 6; ++i) {
        const int px = horizontal ? x - 2 + i : x;
        const int py = horizontal ? y : y - 2 + i;
        const uint8_t* q = img.at(py * img.w + px);
        for (int k = 0; k < 4; ++k) c[i][k] = static_cast<float>(q[k]) * kInv255f;
    }
    for (int i = 0; i < 5; ++i) {
        const float* a = c[i];
        const float* b = c[i + 1];
        const double next = (static_cast<double>(b[2]) * 0.3 + static_cast<double>(b[1]) * 0.59) +
                            static_cast<double>(b[0]) * 0.11;
        const double prev = (static_cast<double>(a[2]) * 0.3 + static_cast<double>(a[1]) * 0.59) +
                            static_cast<double>(a[0]) * 0.11;
        p[i] = static_cast<float>(next - prev);
    }
    const float p2 = p[2] < 0.0f ? -p[2] : p[2];
    if (p[2] < 0.0f)
        for (int i = 0; i < 5; ++i) p[i] = -p[i];
    float run;
    if (0.0f > p[3]) {
        p[3] = 0.0f;
        p[4] = 0.0f;
        if (0.0f > p[1]) {
            p[1] = 0.0f;
            p[0] = 0.0f;
            run = 0.0f;
        } else {
            goto check0;
        }
    } else if (0.0f > p[4]) {
        p[4] = 0.0f;
        if (0.0f > p[1]) {
            p[1] = 0.0f;
            p[0] = 0.0f;
            run = 0.0f;
        } else {
            goto check0;
        }
    } else if (0.0f > p[1]) {
        p[1] = 0.0f;
        p[0] = 0.0f;
        run = 0.0f;
    } else {
    check0:
        run = p[0];
        if (0.0f > run) {
            p[0] = 0.0f;
            run = 0.0f;
        }
    }
    for (int i = 1; i < 5; ++i) run = p[i] > run ? p[i] : run;  // MAXSS p[i], run
    const bool strong = (p2 > threshold) && (run == p2);
    float s = 0.0f;
    for (int i = 0; i < 5; ++i) s = s + p[i];
    for (int i = 0; i < 5; ++i) p[i] = p[i] / s;
    return strong;
}

void computeEdgeProfileFeature(const Image& img, float out[5]) {
    for (int k = 0; k < 5; ++k) out[k] = 0.0f;
    const int w = img.w, h = img.h;
    if (!(h > 5 && w > 5)) return;
    int count = 0;
    float p[5];
    for (int y = 2; y != h - 3; ++y) {
        for (int x = 2; x != w - 3; ++x) {
            if (isStrongEdge(img, x, y, true, 0.1f, p)) {
                for (int k = 0; k < 5; ++k) out[k] = out[k] + p[k];
                ++count;
            }
            if (isStrongEdge(img, x, y, false, 0.1f, p)) {
                for (int k = 0; k < 5; ++k) out[k] = out[k] + p[k];
                ++count;
            }
        }
    }
    if (count != 0)
        for (int k = 0; k < 5; ++k) out[k] = out[k] / static_cast<float>(count);
}

inline float seqSum4(const float v[4]) {
    float s = 0.0f;
    s = v[0] + s;
    s = v[1] + s;
    s = v[2] + s;
    s = v[3] + s;
    return s;
}

}  // namespace

// ---------------------------------------------------------------------------------------------
// ParameterClassifier::classifyParameters

void computeImageTypeFeatures(const Image& img, PixelMap& map, PixelSegmentation& seg,
                              const std::vector<uint8_t>& edges, const EdgeParams& ep, float f[20]) {
    for (int k = 0; k < 20; ++k) f[k] = 0.0f;
    const int w = img.w, h = img.h, n = w * h;
    const SignalImages sig = computeSignals(img);
    const std::vector<uint8_t> mask = boundaryMask8(map, 2);
    auto weightAt = [&](int i) { return edgeWeight(edges[i], ep.scale, ep.floor); };

    ColorHistogram hist;
    hist.setNumBins16();
    for (int i = 0; i < n; ++i) hist.addColor(img.at(i), weightAt(i));
    std::vector<std::array<uint8_t, 4>> peaks;
    std::vector<double> peakW;
    hist.findPeaks(peaks, peakW);
    int numPeaks = static_cast<int>(peaks.size());
    const double total0 = hist.totalWeight();
    for (const auto& c : peaks) hist.zeroOutColor(c.data(), false);
    const double totalPeaks = hist.totalWeight();
    const ColorHistogram backup = hist;
    double totalNeighborhoods, totalLines;
    if (numPeaks < 1) {
        totalNeighborhoods = hist.totalWeight();
        hist = backup;
        totalLines = hist.totalWeight();
        peaks.assign(1, {0, 0, 0, 0});
        peakW.assign(1, static_cast<double>(n));
        numPeaks = 1;
    } else {
        for (const auto& c : peaks) hist.zeroOutColor(c.data(), true);
        totalNeighborhoods = hist.totalWeight();
        hist = backup;
        for (int i = 0; i + 1 < numPeaks; ++i) {
            for (int j = i + 1; j < numPeaks; ++j) {
                float a[4], d[4];
                for (int k = 0; k < 4; ++k) a[k] = static_cast<float>(peaks[i][k]) * kInv255f;
                for (int k = 0; k < 4; ++k) d[k] = static_cast<float>(peaks[j][k]) * kInv255f - a[k];
                float sq[4];
                for (int k = 0; k < 4; ++k) sq[k] = d[k] * d[k];
                const float len = std::sqrt(seqSum4(sq));
                const int steps = static_cast<int>(len * 16.0f);
                const float denom = static_cast<float>(steps + 2);
                for (int s = 0; s < steps + 3; ++s) {
                    const float t = static_cast<float>(s) / denom;
                    uint8_t c[4];
                    for (int k = 0; k < 4; ++k) {
                        const float v = t * d[k] + a[k];
                        c[k] = static_cast<uint8_t>(static_cast<int>(v * 255.0f + 0.5f));
                    }
                    hist.zeroOutColor(c, false);
                }
            }
        }
        totalLines = hist.totalWeight();
    }
    const int cap = std::min(numPeaks, 12);

    std::vector<uint8_t> label(static_cast<size_t>(n));
    std::vector<std::array<float, 4>> clusterSum(cap, {0, 0, 0, 0}), clusterMean(cap, {0, 0, 0, 0});
    std::vector<float> clusterW(cap, 0.0f);
    float totalW = 0.0f, resAll = 0.0f, resBoundary = 0.0f, resInterior = 0.0f;
    int cntBoundary = 0, cntInterior = 0;
    std::vector<uint32_t> distinct(static_cast<size_t>(n));
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            const int i = y * w + x;
            const uint8_t* px = img.at(i);
            float c[4];
            for (int k = 0; k < 4; ++k) c[k] = static_cast<float>(px[k]) * kInv255f;
            int best = 0;
            if (numPeaks > 1) {
                // Distance to the peak colour in byte units (the peak is not scaled).
                auto dist = [&](int p) {
                    float sq[4];
                    for (int k = 0; k < 4; ++k) {
                        const float t = c[k] - static_cast<float>(peaks[p][k]);
                        sq[k] = t * t;
                    }
                    return static_cast<double>(seqSum4(sq));
                };
                double bd = dist(0);
                for (int p = 1; p < cap; ++p) {
                    const double d = dist(p);
                    if (bd > d) best = p;
                    bd = d < bd ? d : bd;
                }
            }
            label[i] = static_cast<uint8_t>(best);
            const float wt = static_cast<float>(weightAt(i));
            totalW = totalW + wt;
            for (int k = 0; k < 4; ++k) clusterSum[best][k] = clusterSum[best][k] + wt * c[k];
            clusterW[best] = clusterW[best] + wt;
            const float s1 = sig.s1[i];
            f[0] = s1 * s1 + f[0];
            const float s2 = sig.s2[i];
            f[1] = s2 * s2 + f[1];
            PixelSegment& sg = seg.segment(map.idOf(x, y));
            const float cnt = static_cast<float>(sg.n);
            float r[4];
            for (int k = 0; k < 4; ++k) {
                const float d = c[k] - sg.sum[k] / cnt;
                r[k] = d * d;
            }
            const float res = seqSum4(r);
            std::memcpy(&distinct[i], px, 4);
            if (mask[i] == 0) {
                ++cntInterior;
                resInterior = resInterior + res;
            } else {
                ++cntBoundary;
                resBoundary = resBoundary + res;
            }
            resAll = resAll + res;
        }
    }
    for (int p = 0; p < cap; ++p)
        for (int k = 0; k < 4; ++k) clusterMean[p][k] = clusterSum[p][k] / clusterW[p];

    float m2 = 1e-6f, m4 = 0.0f;
    for (int i = 0; i < n; ++i) {
        const uint8_t* px = img.at(i);
        const float wt = static_cast<float>(weightAt(i));
        const auto& m = clusterMean[label[i]];
        float d = 0.0f;
        for (int k = 0; k < 4; ++k) d = (static_cast<float>(px[k]) * kInv255f - m[k]) + d;
        const float d2 = d * d;
        m2 = m2 + wt * d2;
        m4 = m4 + (d2 * d2) * wt;
    }
    m2 = m2 / totalW;
    m4 = m4 / totalW;
    const float kurt = m4 / (m2 * m2) - 3.0f;

    if (cntInterior != 0) {
        resInterior = resInterior / static_cast<float>(cntInterior);
        if (cntBoundary != 0) {
            resBoundary = resBoundary / static_cast<float>(cntBoundary);
            resAll = resAll / static_cast<float>(cntBoundary + cntInterior);
        } else {
            resAll = resAll / static_cast<float>(cntInterior);
        }
    } else if (cntBoundary != 0) {
        resBoundary = resBoundary / static_cast<float>(cntBoundary);
        resAll = resAll / static_cast<float>(cntBoundary);
    }
    const float sum = resInterior + resBoundary;
    f[5] = sum > 0.0f ? resBoundary / sum : 0.0f;
    f[6] = resAll;
    f[7] = resBoundary;
    f[8] = resInterior;
    const float nf = static_cast<float>(n);
    f[0] = std::sqrt(f[0] / nf);
    f[1] = std::sqrt(f[1] / nf);
    const float segs = static_cast<float>(map.maxId() + 1);
    f[9] = nf;
    f[10] = segs;
    f[11] = segs / nf;
    f[12] = segs / std::sqrt(nf);
    float prof[5];
    computeEdgeProfileFeature(img, prof);
    f[2] = (prof[4] + prof[0]) * 0.5f;
    f[3] = (prof[3] + prof[1]) * 0.5f;
    f[4] = prof[2];
    std::sort(distinct.begin(), distinct.end());
    f[13] = static_cast<float>(std::unique(distinct.begin(), distinct.end()) - distinct.begin());
    f[14] = kurt;
    f[15] = m2;
    f[16] = m4;
    f[17] = static_cast<float>(totalPeaks / total0);
    f[18] = static_cast<float>(totalNeighborhoods / total0);
    f[19] = static_cast<float>(totalLines / total0);
}

// ---------------------------------------------------------------------------------------------

Classification classifyImage(const Image& img, const std::string& ext, PixelMap& map) {
    // Image classification: segmentation with the classify parameter set (not anti-aliased,
    // flexible palette): overSegmentation(false) + generatePaletteFromSegments.
    PixelSegmenterParams p;
    p.max_iterations = 5;
    p.lambda_pre = 1e-5f;
    p.lambda_initial = 1e-4f;
    p.lambda_final = 0.1f;
    p.min_num_pixels = 0;
    p.num_iter_btwn_renum = 1;
    p.rise_fraction = 0.4f;
    Palette palette;
    map.setSize(img.w, img.h);
    PixelSegmentation seg(map, img);
    seg.setPalette(&palette);
    PixelSegmenter segmenter(seg, p);
    segmenter.initialize();
    segmenter.overSegmentation(false);
    seg.generatePaletteFromSegments();
    segmenter.enforceMinNumPixels();

    // ParameterClassifier::classifyImage
    Classification r;
    const EdgeParams ep;
    const std::vector<uint8_t> edges = edgeImage(img);
    computeImageTypeFeatures(img, map, seg, edges, ep, r.typeFeatures);
    r.typeString = classifyImageType(r.typeFeatures, ext);
    const bool antiAliased = r.typeString.rfind("logo ", 0) != 0;  // Shared::is_anti_aliased = type != logo
    const Image small = downsamplePixels(img);
    map.setSize(small.w, small.h);  // resize
    PaletteFinder finder(small, map, edges, img.w, antiAliased, PaletteFinderParams{}, ep);
    r.paletteString = finder.findPalettes(ext, r.typeString);
    std::memcpy(r.paletteFeatures, finder.features(), sizeof r.paletteFeatures);
    for (int k = 2; k <= 12; ++k) r.palettes[k] = finder.palette(k);

    // ImageClassification::setImageTypeAndQualityFromString / setPaletteModeFromString (the maps
    // image_classification_type_map etc. as filled at run time).
    const std::string type = r.typeString.substr(0, r.typeString.find(' '));
    r.imageType = type == "logo" ? 2 : type == "logo_aa" ? 1 : 0;
    r.quality = 2 - (r.typeString.back() - '0');
    const int colors = std::stoi(r.paletteString);
    r.paletteMode = colors == 2 ? 0 : colors == 13 ? 2 : 1;
    r.maxColors = colors == 13 ? 12 : colors;
    return r;
}

void findPalettesStage(const Image& img, PixelMap& map, bool antiAliased, float costPerColor) {
    const Image small = downsamplePixels(img);
    map.setSize(small.w, small.h);
    // EdgeImage::update runs on the downsampled pixels only when anti-aliased; the aliased k-means
    // never reads the edge image.
    const std::vector<uint8_t> edges = antiAliased ? edgeImage(small) : std::vector<uint8_t>{};
    PaletteFinderParams pp;
    pp.costPerColor = costPerColor;
    PaletteFinder finder(small, map, edges, small.w, antiAliased, pp, EdgeParams{});
    finder.findPalettesRange(2, 12);
}

}  // namespace isvg
