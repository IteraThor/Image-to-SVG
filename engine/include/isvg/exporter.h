// Exporter: compound curve assembly and the SVG writer.
#pragma once

#include <string>
#include <utility>
#include <vector>

#include "isvg/bezier_fitter.h"
#include "isvg/vector_image.h"

namespace isvg {

class Exporter {
public:
    Exporter(VectorImage& vi, BezierFitter& bf) : vi_(vi), bf_(bf) {}

    // For each region, the closed loops (holes) formed by the curves its child regions share with
    // it, as (child region, curve reference index) pairs.
    void generateCompoundCurves();
    // SVG text as Exporter::writeSVG produces it with compound paths on (mode 1) and ptSize 72.
    // background (optional, per region): 1 = made background, i.e. left out (export settings).
    // colors (optional, per region): 0x01RRGGBB replaces the region's fill colour (an edit made after
    // vectorizing; the opacity is kept), 0 keeps it.
    std::string writeSvg(const std::vector<uint8_t>* background = nullptr,
                         const std::vector<uint32_t>* colors = nullptr);
    // Region index of each <path> writeSvg would emit with no background shapes.
    std::vector<int> exportedRegions() const;

private:
    void writeLineOrCurve(std::string& out, int region, const CurveReference& ref, int startFlag, bool hole);

    VectorImage& vi_;
    BezierFitter& bf_;
    std::vector<std::vector<std::vector<std::pair<int, int>>>> compound_;
};

}  // namespace isvg
