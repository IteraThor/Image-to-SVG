// The whole automatic vectorization, in pipeline order:
// loadImageFromMemory -> classifyImage -> findPalettes(2, 12) -> segmentImage -> contourSmoothImage
// -> bezierFitImage -> Exporter::writeSVG (compound paths, 72 dpi).
#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "isvg/bezier_fitter.h"
#include "isvg/exporter.h"
#include "isvg/parameter_classifier.h"
#include "isvg/vector_image.h"

namespace isvg {

// The finished vectorization, kept so the SVG can be written again with shapes made background
// (the "Remove background" edit).
struct VectorizeSession {
    VectorImage vimg;
    std::unique_ptr<BezierFitter> fitter;
    std::unique_ptr<Exporter> exporter;
};

struct VectorizeResult {
    std::string svg;
    Classification classification;
    std::string paramFile;  // e.g. "EP_LOGO_HIG_AA"
    int numShapes = 0, numCurves = 0;
    int width = 0, height = 0;  // after rescaling
    std::vector<int> pathRegions;  // region index of each <path> in svg, in order
    std::shared_ptr<VectorizeSession> session;
};

// SVG with the given regions (background[region] != 0) left out, identical to the unedited
// export after "Remove background" on those shapes.
// colors: per-region fill overrides (0x01RRGGBB, 0 = unchanged; see Exporter::writeSvg).
std::string exportSvg(const VectorizeResult& r, const std::vector<uint8_t>& background,
                      const std::vector<uint32_t>& colors = {});

// stage: "classify", "palettes", "segment", "trace", "smooth", "bezier", "export"; fraction in [0, 1].
using ProgressFn = std::function<void(const char* stage, double fraction)>;

// Rescaling, applied to the premultiplied pixels: when the image has more than
// maxPixels pixels, box-average f x f blocks (f = ceil(sqrt(w*h / maxPixels))) into a
// ceil(w/f) x ceil(h/f) image. Returns false when no rescaling was needed. The web app uses 1600000
// (shrinkImageAsNeeded, the "Small" choice) and a memory-dependent cap of at most 20000000.
bool rescaleImage(Image& img, int maxPixels);

// bgra: w*h pixels, 4 bytes each in B,G,R,A order, not premultiplied (Qt ARGB32 in memory).
// ext: the file extension the classifier sees (".png", ".jpg", ...). The global PRNG is reset first,
// matching a freshly started application.
// maxPixels > 0 applies rescaleImage first.
VectorizeResult vectorize(const uint8_t* bgra, int w, int h, const std::string& ext,
                          const ProgressFn& progress = {}, int maxPixels = 0);

}  // namespace isvg
