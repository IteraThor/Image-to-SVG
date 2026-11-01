// Palette: colours as float RGBA (stored channel order) plus a byte-colour → index map
// .
#pragma once

#include <array>
#include <cstdint>
#include <map>
#include <vector>

#include "isvg/image.h"

namespace isvg {

class Palette {
public:
    int size() const { return static_cast<int>(colors_.size()); }
    const float* color(int i) const { return colors_.at(i).data(); }
    void setNumColors(int n) { colors_.resize(n); }
    std::array<float, 4>& colorRef(int i) { return colors_.at(i); }

    void clearMap() { map_.clear(); }
    // Index of the colour (4 bytes read as a little-endian uint32); new colours get the next index.
    int tryAddColorToMap(const uint8_t c[4]) {
        const uint32_t key = c[0] | (c[1] << 8) | (c[2] << 16) | (static_cast<uint32_t>(c[3]) << 24);
        auto it = map_.find(key);
        if (it != map_.end()) return it->second;
        const int idx = static_cast<int>(map_.size());
        map_.emplace(key, idx);
        return idx;
    }
    void setNewPaletteFromMap() {
        colors_.resize(map_.size());
        for (const auto& [key, idx] : map_)
            for (int k = 0; k < 4; ++k)
                colors_[idx][k] = static_cast<float>((key >> (8 * k)) & 0xff) * kInv255f;
    }

private:
    std::vector<std::array<float, 4>> colors_;  // / count
    std::map<uint32_t, int> map_;
};

}  // namespace isvg
