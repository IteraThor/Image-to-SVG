#include "isvg/pixel_map.h"

#include <algorithm>
#include <cstring>

namespace isvg {

// PixelMap::setSize
void PixelMap::setSize(int w, int h) {
    if (w != w_ || h != h_) {
        w_ = w;
        h_ = h;
        ids_.assign(std::max(0, w * h), 0);
    }
    maxId_ = -1;
    live_ = 0;
}

// PixelMap::idOfWithReplace(int,int,int) with replace disabled.
int PixelMap::idOf(int index) {
    const int first = ids_[index];
    int id = first;
    if (redirects_[first] >= 0) {
        do {
            id = redirects_[id];
        } while (redirects_[id] >= 0);
        redirects_[first] = id;
        ids_[index] = id;
    }
    return id;
}

// PixelMap::redirect
void PixelMap::redirect(int into, int from) {
    dirty_ = true;
    redirects_[from] = into;
    --live_;
}

// PixelMap::setMaxId — the redirect table is only cleared when it has to grow.
void PixelMap::setMaxId(int maxId) {
    maxId_ = maxId;
    live_ = maxId + 1;
    if (static_cast<int>(redirects_.size()) <= maxId) {
        redirects_.assign(std::max(0, maxId + 1), -1);
        dirty_ = false;
    }
}

// PixelMap::resetRedirects
void PixelMap::resetRedirects() {
    dirty_ = false;
    std::fill(redirects_.begin(), redirects_.end(), -1);
}

// PixelMap::flushRedirectsIfDirty
void PixelMap::flushRedirectsIfDirty() {
    if (!dirty_) return;
    for (int i = 0; i < w_ * h_; ++i) idOf(i);
    dirty_ = false;
    std::fill(redirects_.begin(), redirects_.end(), -1);
}

// PixelMap::compact — resize the redirect table to maxId+1, keeping contents.
void PixelMap::compact() {
    const int n = maxId_ + 1;
    if (n != static_cast<int>(redirects_.size())) redirects_.resize(std::max(0, n));
}

// PixelMap::floodFill. The original is a scanline fill; the set of pixels it
// rewrites is exactly the 4-connected component of raw value `oldId` containing `start`.
void PixelMap::floodFill(int oldId, int newId, int start) {
    if (ids_[start] != oldId || oldId == newId) return;
    std::vector<int> stack{start};
    ids_[start] = newId;
    while (!stack.empty()) {
        const int p = stack.back();
        stack.pop_back();
        const int x = p % w_, y = p / w_;
        const int nb[4] = {x + 1 < w_ ? p + 1 : -1, y + 1 < h_ ? p + w_ : -1, x > 0 ? p - 1 : -1,
                           y > 0 ? p - w_ : -1};
        for (int q : nb) {
            if (q >= 0 && ids_[q] == oldId) {
                ids_[q] = newId;
                stack.push_back(q);
            }
        }
    }
}

}  // namespace isvg
