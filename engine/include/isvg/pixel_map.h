// PixelMap: per-pixel segment ids with lazy merge redirects (union-find style).
// Keeps the field layout and redirect timing exact so that results are reproducible — including
// which ids survive renumbering and when stale redirects are visible — carries over.
#pragma once

#include <cstdint>
#include <vector>

namespace isvg {

class PixelMap {
public:
    void setSize(int w, int h);
    int width() const { return w_; }
    int height() const { return h_; }
    int numPixels() const { return w_ * h_; }

    // idOf follows redirects to the root id and compresses the first hop (PixelMap::idOfWithReplace).
    int idOf(int index);
    int idOf(int x, int y) { return (x < 0 || x >= w_ || y < 0 || y >= h_) ? -1 : idOf(y * w_ + x); }
    // idOfWithReplace(int,int,int,int): -1 outside the image; `from` reads as `to`.
    int idOfReplace(int x, int y, int to, int from) {
        if (x < 0 || x >= w_ || y < 0 || y >= h_) return -1;
        const int id = idOf(y * w_ + x);
        return (to >= 0 && id == from) ? to : id;
    }

    int& raw(int index) { return ids_[index]; }
    const std::vector<int32_t>& rawIds() const { return ids_; }

    // Merge `from` into `into`: later idOf lookups of `from` resolve to `into`.
    void redirect(int into, int from);
    void setMaxId(int maxId);
    void resetRedirects();
    void flushRedirectsIfDirty();
    void compact();
    // Scanline flood fill over raw ids: replace the 4-connected run of `oldId` containing `start`.
    void floodFill(int oldId, int newId, int start);

    int redirectOf(int id) const { return redirects_[id]; }
    int liveCount() const { return live_; }
    void decrementLive() { --live_; }
    int maxId() const { return maxId_; }

private:
    std::vector<int32_t> ids_;
    int w_ = 0, h_ = 0;
    int live_ = 0;
    std::vector<int32_t> redirects_;  // (size == capacity
    bool dirty_ = false;
    int maxId_ = 0;
};

}  // namespace isvg
