// Engine parameter files (embedded text, parsed by
// the parameter loader) and their selection by populateSettings.
#pragma once

#include <string>

#include "isvg/bezier_fitter.h"
#include "isvg/contour_smoother.h"
#include "isvg/palette_finder.h"
#include "isvg/parameter_classifier.h"
#include "isvg/pixel_segmenter.h"
#include "isvg/subpixel_segmentation.h"
#include "isvg/subpixel_segmenter.h"

namespace isvg {

struct EngineParamFile {
    const char* name;  // "EP_LOGO_HIG_AA"
    const char* text;
};
extern const EngineParamFile kEngineParamFiles[];  // engine_params_data.cpp (generated)

struct EngineParams {
    bool antiAliased = true;  // Shared::is_anti_aliased
    int imageTypeCode = 5;  // Shared::image_type_code
    PaletteFinderParams paletteFinder;
    EdgeParams edge;
    PixelSegmenterParams segmenter;
    bool useOverSegDiagExtr = false;  // Segmenter::use_over_seg_diag_extr
    SubPixelSegmenterParams subSegmenter;
    SubPixelSegmentationParams subSegmentation;
    ContourSmootherParams smoother;
    BezierFitterParams bezier;
};

// Text of a parameter file by name; nullptr if unknown.
const char* engineParamText(const std::string& name);
// populateSettings' ep_files table: {PHOTO,LOGO_*_AA,LOGO}_{HIG,MED,LOW}[base + quality] with base
// 3 for logo_aa (type 1), 6 for logo (type 2) and 0 otherwise.
std::string engineParamFileFor(int imageType, int quality);
// "Name ...... value [value value]" lines; floats are read with strtof, doubles with strtod.
EngineParams parseEngineParams(const char* text);

}  // namespace isvg
