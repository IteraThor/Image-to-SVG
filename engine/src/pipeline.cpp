#include "isvg/pipeline.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

#include "isvg/bezier_fitter.h"
#include "isvg/boundary_tracer.h"
#include "isvg/contour_smoother.h"
#include "isvg/engine_params.h"
#include "isvg/exporter.h"
#include "isvg/pixel_segmenter.h"
#include "isvg/prng.h"
#include "isvg/subpixel_segmentation.h"
#include "isvg/subpixel_segmenter.h"

namespace isvg {

bool rescaleImage(Image& img, int maxPixels) {
    const int W = img.w, H = img.h;
    // new size
    const double fd = std::ceil(std::sqrt(static_cast<double>(W * H) / static_cast<double>(maxPixels)));
    const double s = 1.0 / fd;
    const int f = static_cast<int>(fd);
    const int outW = static_cast<int>(std::ceil(static_cast<double>(W) * s));
    const int outH = static_cast<int>(std::ceil(static_cast<double>(H) * s));
    if (W * H <= maxPixels) return false;
    Image out;
    out.w = outW;
    out.h = outH;
    out.px.assign(static_cast<size_t>(outW) * outH * 4, 0);
    std::vector<float> acc(static_cast<size_t>(outW) * 4);
    int sy = 0;
    for (int oy = 0; oy < outH; ++oy) {
        std::fill(acc.begin(), acc.end(), 0.0f);
        const int rows = std::min(H - sy, f);
        const double wy = 1.0 / static_cast<double>(rows);
        for (int y = sy; y < sy + rows; ++y) {
            int sx = 0;
            for (int ox = 0; ox < outW; ++ox) {
                const int cols = std::min(W - sx, f);
                const double wx = 1.0 / static_cast<double>(cols);
                const float wgt = static_cast<float>(wy * wx);
                for (int x = sx; x < sx + cols; ++x) {
                    const uint8_t* p = img.at(y * W + x);
                    for (int k = 0; k < 4; ++k)
                        acc[4 * ox + k] = acc[4 * ox + k] + wgt * (static_cast<float>(p[k]) * kInv255f);
                }
                sx += cols;
            }
        }
        sy += rows;
        for (int ox = 0; ox < outW; ++ox)
            for (int k = 0; k < 4; ++k)
                out.px[(static_cast<size_t>(oy) * outW + ox) * 4 + k] =
                    static_cast<uint8_t>(static_cast<int>(acc[4 * ox + k] * 255.0f + 0.5f));
    }
    img = std::move(out);
    return true;
}

VectorizeResult vectorize(const uint8_t* bgra, int w, int h, const std::string& ext, const ProgressFn& progress,
                          int maxPixels) {
    auto report = [&](const char* stage, double f) {
        if (progress) progress(stage, f);
    };
    VectorizeResult result;
    Image img;
    img.w = w;
    img.h = h;
    img.px.assign(bgra, bgra + static_cast<size_t>(w) * h * 4);
    premultiplyOnLoad(img);
    if (maxPixels > 0) rescaleImage(img, maxPixels);
    result.width = img.w;
    result.height = img.h;
    resetPrng();

    // Classification and the palette search over 2..12 colours (only its PRNG use affects the output).
    report("classify", 0.0);
    PixelMap engineMap;
    result.classification = classifyImage(img, ext, engineMap);
    const Classification& cls = result.classification;
    report("palettes", 0.0);
    findPalettesStage(img, engineMap, cls.imageType == 1, 0.15f);
    result.paramFile = engineParamFileFor(cls.imageType, cls.quality);
    const EngineParams ep = parseEngineParams(engineParamText(result.paramFile));

    // Segmentation, flexible palette.
    report("segment", 0.0);
    Palette palette;
    PixelMap map;
    map.setSize(img.w, img.h);
    PixelSegmentation seg(map, img);
    seg.setPalette(&palette);
    PixelSegmenter segmenter(seg, ep.segmenter);
    segmenter.initialize();
    if (!ep.antiAliased) {
        if (ep.useOverSegDiagExtr) {
            segmenter.overSegmentation(true);
            segmenter.clusterColorsAndGeneratePalette();
            seg.updateStatistics();
        } else {
            segmenter.overSegmentation(false);
            seg.generatePaletteFromSegments();
        }
    } else {
        segmenter.overSegmentation(false);
        report("segment", 0.3);
        SubPixelSegmentation sub(map, img, ep.subSegmentation);
        sub.setPalette(&palette);
        SubPixelSegmenter sps(sub, ep.subSegmenter);
        sps.execute();
        sps.clusterColors(0x1a);
        sub.flattenTo(seg);
        report("segment", 0.6);
        sps.flagShadowSegments();
        sps.splitAndMergeFlaggedSegmentsFromBlurryImage(2);
        sps.clusterColors(0x2e);
        sub.flattenTo(seg);
        sps.finalize();
    }
    segmenter.enforceMinNumPixels();

    report("trace", 0.0);
    auto session = std::make_shared<VectorizeSession>();
    VectorImage& vimg = session->vimg;
    buildRegions(seg, vimg.regions);
    BoundaryTracer tracer(map, vimg, ep.antiAliased);
    tracer.execute();

    report("smooth", 0.0);
    ContourSmoother smoother(vimg, map, img, palette, ep.antiAliased, ep.smoother);
    smoother.initialize(true);
    for (int p = 0; p < 3; ++p) {
        smoother.executePhase(p);
        report("smooth", (p + 1) / 3.0);
    }
    smoother.finish();
    smoother.straightenEdges();

    report("bezier", 0.0);
    session->fitter = std::make_unique<BezierFitter>(vimg, ep.bezier);
    session->fitter->execute();

    report("export", 0.0);
    session->exporter = std::make_unique<Exporter>(vimg, *session->fitter);
    session->exporter->generateCompoundCurves();
    result.svg = session->exporter->writeSvg();
    result.pathRegions = session->exporter->exportedRegions();
    result.numShapes = static_cast<int>(vimg.regions.size());
    result.numCurves = session->fitter->numCurves();
    result.session = session;
    report("export", 1.0);
    return result;
}

std::string exportSvg(const VectorizeResult& r, const std::vector<uint8_t>& background,
                      const std::vector<uint32_t>& colors) {
    return r.session->exporter->writeSvg(&background, &colors);
}

}  // namespace isvg
