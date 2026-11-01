// WebAssembly entry points for the web app (web/worker.js).
#include <emscripten.h>

#include <cstdlib>
#include <cstring>
#include <string>

#include "isvg/pipeline.h"

#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_PNG
#define STBI_ONLY_BMP
#define STBI_ONLY_GIF
#define STBI_ONLY_TGA
#include "../third_party/stb/stb_image.h"

EM_JS(void, isvg_report_progress, (const char* stage, double fraction), {
    if (Module.onProgress) Module.onProgress(UTF8ToString(stage), fraction);
});

namespace {
std::string g_svg, g_info, g_export;
isvg::VectorizeResult g_result;  // kept for isvg_export (remove background)
}

extern "C" {

EMSCRIPTEN_KEEPALIVE void* isvg_alloc(int bytes) { return std::malloc(bytes); }
EMSCRIPTEN_KEEPALIVE void isvg_free(void* p) { std::free(p); }

// Decodes PNG/BMP/GIF/TGA file bytes to straight (non-premultiplied) BGRA, as Qt's QImage
// (Format_ARGB32) hands them to the engine. Returns nullptr when the format is not supported here
// (the caller then decodes with the browser, e.g. for JPEG). Free the result with isvg_free.
EMSCRIPTEN_KEEPALIVE uint8_t* isvg_decode(const uint8_t* data, int len, int* w, int* h) {
    int n = 0;
    uint8_t* rgba = stbi_load_from_memory(data, len, w, h, &n, 4);
    if (!rgba) return nullptr;
    const size_t count = static_cast<size_t>(*w) * *h;
    for (size_t i = 0; i < count; ++i) std::swap(rgba[4 * i], rgba[4 * i + 2]);
    return rgba;
}

// Runs the full vectorization. Returns the SVG text (valid until the next call); isvg_info returns
// a JSON summary of the run.
// maxPixels > 0: shrink larger images first (rescaleImage).
EMSCRIPTEN_KEEPALIVE const char* isvg_vectorize(const uint8_t* bgra, int w, int h, const char* ext, int maxPixels) {
    g_result = isvg::VectorizeResult{};
    g_result =
        isvg::vectorize(bgra, w, h, ext, [](const char* stage, double f) { isvg_report_progress(stage, f); }, maxPixels);
    const isvg::VectorizeResult& r = g_result;
    g_svg = r.svg;
    const isvg::Classification& c = r.classification;
    std::string paths = "[";
    for (size_t i = 0; i < r.pathRegions.size(); ++i) paths += (i ? "," : "") + std::to_string(r.pathRegions[i]);
    paths += "]";
    g_info = "{\"type\":\"" + c.typeString + "\",\"palette\":\"" + c.paletteString +
             "\",\"classification\":[" + std::to_string(c.imageType) + "," + std::to_string(c.quality) + "," +
             std::to_string(c.paletteMode) + "," + std::to_string(c.maxColors) + "],\"params\":\"" + r.paramFile +
             "\",\"width\":" + std::to_string(r.width) + ",\"height\":" + std::to_string(r.height) +
             ",\"shapes\":" + std::to_string(r.numShapes) + ",\"curves\":" + std::to_string(r.numCurves) +
             ",\"pathRegions\":" + paths + "}";
    return g_svg.c_str();
}

EMSCRIPTEN_KEEPALIVE const char* isvg_info() { return g_info.c_str(); }

// Re-exports the last result with edits: flags[region] != 0 makes the shape background (removed);
// colors[region] = 0x01RRGGBB recolours it (0 = unchanged). n = region count; colors may be null.
EMSCRIPTEN_KEEPALIVE const char* isvg_export(const uint8_t* flags, const uint32_t* colors, int n) {
    if (!g_result.session) return "";
    g_export = isvg::exportSvg(g_result, std::vector<uint8_t>(flags, flags + n),
                             colors ? std::vector<uint32_t>(colors, colors + n) : std::vector<uint32_t>{});
    return g_export.c_str();
}

}  // extern "C"
