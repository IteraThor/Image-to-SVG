// PixelSegmentation: segment statistics over a PixelMap.
#pragma once

#include <cstdint>
#include <vector>

#include "isvg/image.h"
#include "isvg/palette.h"
#include "isvg/pixel_map.h"

namespace isvg {

// PixelSegment: pixel count, per-channel sums of byte/255, and a
// palette index carried through renumbering.
struct PixelSegment {
    int n = 0;
    float sum[4] = {0, 0, 0, 0};  // (channel order as stored in the image, i.e. B,G,R,A)
    int paletteIndex = -1;

    void clearStats() {
        n = 0;
        sum[0] = sum[1] = sum[2] = sum[3] = 0.0f;
    }
    void addPixelF(const float c[4]) {
        ++n;
        for (int k = 0; k < 4; ++k) sum[k] = sum[k] + c[k];
    }
};

class PixelSegmentation {
public:
    PixelSegmentation(PixelMap& map, const Image& image) : map_(map), image_(image) {}
    virtual ~PixelSegmentation() = default;

    void initialize();
    void updateStatistics();
    void renumber();
    void updateMaxId(int maxId);
    void compact();
    // Resolve 2×2 checkerboard saddles by relabelling one pixel.
    void performDiagonalExtrusions();

    PixelSegment& segment(int id) { return segments_[id]; }
    int numSegmentSlots() const { return static_cast<int>(segments_.size()); }
    PixelMap& map() { return map_; }
    const Image& image() const { return image_; }
    void setPalette(Palette* p) { palette_ = p; }
    Palette* palette() { return palette_; }

    // Segment colour: palette colour if it has a palette index, else the segment's own colour
    // (PixelSegmentation::getColorF.
    void getColorF(int id, float out[4]);
    // Quantise every segment colour to bytes, collect the distinct colours into the palette and
    // give each segment its palette index.
    void generatePaletteFromSegments();

protected:
    // Colour of a segment without a palette index (virtual slot 2 of the segment record): mean for
    // PixelSegment, model colour for SubPixelSegment.
    virtual void ownColorF(int id, float out[4]) {
        const PixelSegment& s = segments_[id];
        const float n = static_cast<float>(s.n);
        for (int k = 0; k < 4; ++k) out[k] = s.sum[k] / n;
    }
    // The original reaches the palette index through the virtual getPixelSegment; the sub-pixel
    // segmentation overrides it with its own segment array.
    virtual bool hasSegmentRecords() const { return !segments_.empty(); }
    virtual int& paletteIndexOf(int id) { return segments_[id].paletteIndex; }

    PixelMap& map_;
    const Image& image_;
    Palette* palette_ = nullptr;  // shared palette
    bool initialized_ = false;
    std::vector<PixelSegment> segments_;  // / count
};

}  // namespace isvg
