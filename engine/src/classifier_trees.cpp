// Generated decision trees (image type and palette mode). Do not edit.
#include "isvg/parameter_classifier.h"

#include <cstdint>
#include <cstring>

#if defined(__clang__)
#pragma clang diagnostic ignored "-Wunused-label"
#endif

namespace isvg {

namespace {
constexpr uint64_t kFeat = 0x100000, kRes = 0x200000;
}  // namespace

std::string classifyImageType(const float* features, const std::string& ext) {
    uint64_t rax = 0, rbx = 0, rcx = kRes, rdx = kFeat, r8 = 0, r9 = 0, r10 = 0, r11 = 0;
    double xd[5] = {0, 0, 0, 0, 0};
    float xf[5] = {0, 0, 0, 0, 0};
    bool above = false, zero = false;
    unsigned char res[32] = {0};
    auto feat = [&](uint64_t a) { float v; std::memcpy(&v, reinterpret_cast<const char*>(features) + (a - kFeat), sizeof v); return v; };
    auto store = [&](uint64_t a, uint64_t v, int n) { if (a >= kRes && a + n <= kRes + 32) std::memcpy(res + (a - kRes), &v, n); };
    (void)rbx; (void)r8; (void)r9; (void)r10; (void)r11; (void)zero;
    xd[1] = 0.0134977;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0x48));
    above = xd[1] > xd[0];  // unordered -> false
    rbx = rcx;
    if (!above) goto L8;
    xd[2] = 0.960049;
    xd[1] = 0;
    xd[1] = static_cast<double>(feat(rdx + 0x10));
    above = xd[2] > xd[1];  // unordered -> false
    if (!above) goto L14;
    xd[3] = 0.0197072;
    xd[2] = 0;
    xd[2] = static_cast<double>(feat(rdx + 0x8));
    above = xd[3] > xd[2];  // unordered -> false
    if (!above) goto L20;
    xd[3] = 0.00219943;
    above = xd[3] > xd[0];  // unordered -> false
    if (!above) goto L24;
    xd[1] = 0.0923396;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0xc));
    above = xd[1] > xd[0];  // unordered -> false
    if (above) goto L30;
    xd[0] = 0.0111265;
    above = xd[0] > xd[2];  // unordered -> false
    if (above) goto L34;
    xd[1] = 0.000665414;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0x20));
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L40;
    xf[0] = 400368.0f;
    above = xf[0] > feat(rdx + 0x24);  // unordered -> false
    if (above) goto L34;
    xd[1] = 3.515e-07;
    xd[0] = 0;
    rax = rcx + 0x10;
    xd[0] = static_cast<double>(feat(rdx + 0x4c));
    store(rcx + 0x0, rax, 8);
    rcx = 0x2061615f6f676f6cull;
    above = xd[1] > xd[0];  // unordered -> false
    store(rbx + 0x10, rcx, 8);
    if (!above) goto L53;
L55:
    store(rbx + 0x18, 0x31ull, 1);
    store(rbx + 0x8, 0x9ull, 8);
    store(rbx + 0x19, 0x0ull, 1);
L58:
    rax = rbx;
    goto done;
L14:
    xd[1] = 5.16e-07;
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L64;
    xd[1] = 0.00130897;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0x8));
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L34;
L70:
    rax = rbx + 0x10;
    store(rbx + 0x10, 0x6f676f6cull, 4);
    store(rbx + 0x0, rax, 8);
    rax = 0x3220ull;
    store(rbx + 0x14, rax, 2);
    store(rbx + 0x8, 0x6ull, 8);
    store(rbx + 0x16, 0x0ull, 1);
L77:
    rax = rbx;
    goto done;
L8:
    xd[2] = 0.867741;
    xd[1] = 0;
    xd[1] = static_cast<double>(feat(rdx + 0x14));
    above = xd[2] > xd[1];  // unordered -> false
    if (!above) goto L85;
    xd[1] = 0.0649332;
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L89;
    xd[1] = 0.0323685;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0x3c));
    above = xd[1] > xd[0];  // unordered -> false
    if (above) goto L34;
    xf[0] = 9497.5f;
    rax = rcx + 0x10;
    above = xf[0] > feat(rdx + 0x28);  // unordered -> false
    store(rcx + 0x0, rax, 8);
    if (!above) goto L100;
L101:
    rax = rbx;
    store(rbx + 0x10, 0x746f6870ull, 4);
    store(rbx + 0x13, 0x31206f74ull, 4);
    store(rbx + 0x8, 0x7ull, 8);
    store(rbx + 0x17, 0x0ull, 1);
    goto done;
L85:
    xd[2] = 0.0124402;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0x8));
    above = xd[2] > xd[0];  // unordered -> false
    if (!above) goto L113;
    xd[1] = 0.984352;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0x10));
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L70;
    xd[1] = 0.091711;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0xc));
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L34;
    xd[1] = 0.00301763;
    xd[0] = 0;
    rax = rcx + 0x10;
    xd[0] = static_cast<double>(feat(rdx + 0x1c));
    store(rcx + 0x0, rax, 8);
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L101;
L131:
    rcx = 0x2061615f6f676f6cull;
    store(rbx + 0x10, rcx, 8);
L53:
    rax = rbx;
    store(rbx + 0x18, 0x32ull, 1);
    store(rbx + 0x8, 0x9ull, 8);
    store(rbx + 0x19, 0x0ull, 1);
    goto done;
L113:
    xd[3] = 0.241446;
    xd[2] = 0;
    xd[2] = static_cast<double>(feat(rdx + 0x0));
    xd[4] = 0;
    xd[4] = static_cast<double>(feat(rdx + 0xc));
    above = xd[3] > xd[2];  // unordered -> false
    if (!above) goto L145;
    xd[2] = 0.183605;
    above = xd[2] > xd[4];  // unordered -> false
    if (!above) goto L149;
    xd[2] = 0.0225216;
    above = xd[2] > xd[0];  // unordered -> false
    if (!above) goto L153;
    xd[0] = 0.915114;
    above = xd[0] > xd[1];  // unordered -> false
    if (above) goto L157;
    xd[1] = 0.649337;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0x10));
    above = xd[1] > xd[0];  // unordered -> false
    if (above) goto L34;
    xf[0] = 80543.5f;
    rax = rcx + 0x10;
    above = xf[0] > feat(rdx + 0x34);  // unordered -> false
    store(rcx + 0x0, rax, 8);
    if (!above) goto L131;
    goto L101;
L89:
    xd[1] = 0.20806;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0x0));
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L173;
    xd[1] = 0.158062;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0x3c));
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L157;
    xd[1] = 0.619184;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0x10));
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L157;
    xd[1] = 0.00190362;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0x1c));
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L157;
    xd[1] = 2.355e-05;
    xd[0] = 0;
    rax = rcx + 0x10;
    xd[0] = static_cast<double>(feat(rdx + 0x4));
    store(rcx + 0x0, rax, 8);
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L101;
    goto L131;
L20:
    xd[3] = 0.00291393;
    above = xd[3] > xd[0];  // unordered -> false
    if (!above) goto L199;
    xd[3] = 0.0285461;
    above = xd[3] > xd[2];  // unordered -> false
    if (!above) goto L40;
    xd[3] = 0.674069;
    above = xd[3] > xd[1];  // unordered -> false
    if (!above) goto L206;
    xd[0] = 0.0238979;
    rax = rcx + 0x10;
    store(rcx + 0x0, rax, 8);
    rcx = 0x2061615f6f676f6cull;
    above = xd[0] > xd[2];  // unordered -> false
    store(rbx + 0x10, rcx, 8);
    if (!above) goto L55;
    goto L53;
L64:
    xd[1] = 0.380912;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0x0));
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L40;
    xd[1] = 0.00114401;
    xd[0] = 0;
    rax = rcx + 0x10;
    xd[0] = static_cast<double>(feat(rdx + 0x8));
    store(rcx + 0x0, rax, 8);
    rcx = 0x2061615f6f676f6cull;
    above = xd[1] > xd[0];  // unordered -> false
    store(rbx + 0x10, rcx, 8);
    if (!above) goto L55;
    goto L53;
L153:
    xd[1] = 0.103327;
    above = xd[1] > xd[0];  // unordered -> false
    if (above) goto L157;
L40:
    rax = rbx + 0x10;
    store(rbx + 0x0, rax, 8);
L100:
    rcx = 0x2061615f6f676f6cull;
    store(rbx + 0x18, 0x31ull, 1);
    store(rbx + 0x10, rcx, 8);
    store(rbx + 0x8, 0x9ull, 8);
    store(rbx + 0x19, 0x0ull, 1);
    goto L58;
L157:
    rax = rbx + 0x10;
    store(rbx + 0x0, rax, 8);
    goto L101;
L24:
    xd[3] = 0.406237;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0x0));
    above = xd[3] > xd[0];  // unordered -> false
    if (above) goto L244;
L34:
    rax = rbx + 0x10;
    rcx = 0x2061615f6f676f6cull;
    store(rbx + 0x0, rax, 8);
    store(rbx + 0x10, rcx, 8);
    goto L53;
L199:
    xd[1] = 2.63168;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0x30));
    above = xd[1] > xd[0];  // unordered -> false
    if (above) goto L254;
    xd[0] = 0.0380126;
    above = xd[0] > xd[2];  // unordered -> false
    if (!above) goto L40;
    xd[1] = 0.400423;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0x0));
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L40;
    xd[1] = 0.136656;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0xc));
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L40;
    xd[1] = 0.000193388;
    xd[0] = 0;
    rax = rcx + 0x10;
    xd[0] = static_cast<double>(feat(rdx + 0x20));
    store(rcx + 0x0, rax, 8);
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L131;
    goto L101;
L173:
    xd[1] = 0.928628;
    rax = rcx + 0x10;
    store(rcx + 0x0, rax, 8);
    above = xd[1] > xd[0];  // unordered -> false
    if (above) goto L101;
    rcx = 0x2061615f6f676f6cull;
    store(rbx + 0x10, rcx, 8);
L282:
    store(rbx + 0x18, 0x30ull, 1);
    store(rbx + 0x8, 0x9ull, 8);
    store(rbx + 0x19, 0x0ull, 1);
    goto L77;
L145:
    xd[0] = 0.132534;
    rax = rcx + 0x10;
    store(rcx + 0x0, rax, 8);
    above = xd[0] > xd[4];  // unordered -> false
    if (!above) goto L131;
    goto L101;
L244:
    xd[0] = 0.0162755;
    rax = rcx + 0x10;
    store(rcx + 0x0, rax, 8);
    above = xd[0] > xd[2];  // unordered -> false
    if (!above) goto L101;
    goto L131;
L254:
    xd[0] = 0.0540923;
    rax = rcx + 0x10;
    store(rcx + 0x0, rax, 8);
    above = xd[0] > xd[2];  // unordered -> false
    if (!above) goto L100;
    goto L101;
L30:
    xd[1] = 0.111529;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0x40));
    above = xd[1] > xd[0];  // unordered -> false
    if (above) goto L34;
    xd[1] = 0.00792502;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0x2c));
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L34;
    xf[0] = 57.5f;
    rax = rcx + 0x10;
    above = xf[0] > feat(rdx + 0x34);  // unordered -> false
    store(rcx + 0x0, rax, 8);
    rcx = 0x2061615f6f676f6cull;
    store(rbx + 0x10, rcx, 8);
    if (!above) goto L55;
    goto L282;
L149:
    xd[1] = 0.00268996;
    xd[0] = 0;
    xd[0] = static_cast<double>(feat(rdx + 0x1c));
    above = xd[1] > xd[0];  // unordered -> false
    if (above) goto L157;
    xd[1] = 1.19829;
    xd[0] = 0;
    rax = rcx + 0x10;
    xd[0] = static_cast<double>(feat(rdx + 0x3c));
    store(rcx + 0x0, rax, 8);
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L101;
    goto L100;
L206:
    xd[1] = 0.00132773;
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L34;
    if (ext.size() == 4) { if (ext == ".bmp" || ext == ".png") goto L40; goto L34; }
    if (ext.size() != 5) goto L34;
    if (ext == ".jpeg" || ext == ".tiff") goto L40;
    goto L34;
done:
    uint64_t len;
    std::memcpy(&len, res + 8, 8);
    return std::string(reinterpret_cast<const char*>(res + 16), static_cast<size_t>(len));
}

std::string classifyPaletteMode(const double* features, const std::string& typeString) {
    uint64_t rax = 0, rbx = 0, rcx = kRes, rdx = kFeat, r8 = 0, r9 = 0, r10 = 0, r11 = 0;
    double xd[5] = {0, 0, 0, 0, 0};
    float xf[5] = {0, 0, 0, 0, 0};
    bool above = false, zero = false;
    unsigned char res[32] = {0};
    auto feat = [&](uint64_t a) { double v; std::memcpy(&v, reinterpret_cast<const char*>(features) + (a - kFeat), sizeof v); return v; };
    auto store = [&](uint64_t a, uint64_t v, int n) { if (a >= kRes && a + n <= kRes + 32) std::memcpy(res + (a - kRes), &v, n); };
    (void)rbx; (void)r8; (void)r9; (void)r10; (void)r11; (void)zero;
    xd[1] = 2.69034;
    xd[0] = feat(rdx + 0x190);
    above = xd[1] > xd[0];  // unordered -> false
    rax = rcx;
    if (!above) goto L360;
    xd[2] = 1.73555;
    xd[1] = feat(rdx + 0x198);
    above = xd[2] > xd[0];  // unordered -> false
    if (!above) goto L365;
    xd[0] = 217.677;
    above = xd[0] > xd[1];  // unordered -> false
    if (above) goto L369;
L370:
    rdx = rax + 0x10;
    store(rax + 0x8, 0x1ull, 8);
    store(rax + 0x0, rdx, 8);
L373:
    r9 = 0x32ull;
    store(rax + 0x10, r9, 2);
    goto done;
L365:
    xd[2] = 995.704;
    above = xd[2] > xd[1];  // unordered -> false
    if (!above) goto L370;
    xd[1] = 2.5;
    above = xd[1] > feat(rdx + 0x120);  // unordered -> false
    if (above) goto L370;
    xd[1] = 0.351533;
    above = xd[1] > feat(rdx + 0x138);  // unordered -> false
    if (above) goto L384;
    xd[1] = 0.938407;
    above = xd[1] > feat(rdx + 0x68);  // unordered -> false
    if (!above) goto L388;
    xd[1] = 2.52056;
    rdx = rcx + 0x10;
    store(rcx + 0x0, rdx, 8);
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L394;
    goto L396;
L360:
    xd[1] = feat(rdx + 0x140);
    xd[2] = 0.288933;
    above = xd[2] > xd[1];  // unordered -> false
    if (!above) goto L400;
    xd[0] = feat(rdx + 0x138);
    xd[2] = 0.333629;
    above = xd[2] > xd[0];  // unordered -> false
    if (!above) goto L405;
L384:
    rdx = rax + 0x10;
    store(rax + 0x8, 0x1ull, 8);
    store(rax + 0x0, rdx, 8);
L408:
    r8 = 0x33ull;
    store(rax + 0x10, r8, 2);
    goto done;
L400:
    xd[1] = 2.04511;
    above = xd[1] > feat(rdx + 0x68);  // unordered -> false
    if (!above) goto L413;
    xd[0] = 7.5;
    above = xd[0] > feat(rdx + 0x78);  // unordered -> false
    if (!above) goto L417;
    xd[0] = 2004.84;
    above = xd[0] > feat(rdx + 0x198);  // unordered -> false
    if (!above) goto L421;
    xd[0] = 6.5;
    above = xd[0] > feat(rdx + 0x120);  // unordered -> false
    if (!above) goto L425;
    goto L388;
L388:
    rdx = rax + 0x10;
    store(rax + 0x0, rdx, 8);
L396:
    rdx = 0x3331ull;
    store(rax + 0x8, 0x2ull, 8);
    store(rax + 0x10, rdx, 2);
    store(rax + 0x12, 0x0ull, 1);
    goto done;
L369:
    xd[0] = 35.25;
    above = xd[0] > feat(rdx + 0x128);  // unordered -> false
    rdx = rcx + 0x10;
    store(rcx + 0x8, 0x1ull, 8);
    store(rcx + 0x0, rdx, 8);
    if (!above) goto L373;
    goto L408;
L413:
    xd[1] = 7.5;
    above = xd[1] > feat(rdx + 0x120);  // unordered -> false
    if (!above) goto L388;
    xd[1] = 3.56775;
    above = xd[1] > xd[0];  // unordered -> false
    if (above) goto L449;
    xd[0] = 0.160751;
    above = xd[0] > feat(rdx + 0x148);  // unordered -> false
    if (!above) goto L388;
L453:
    rdx = rax + 0x10;
    rcx = 0x35ull;
    store(rax + 0x8, 0x1ull, 8);
    store(rax + 0x0, rdx, 8);
    store(rax + 0x10, rcx, 2);
    goto done;
L405:
    xd[2] = 6.5;
    above = xd[2] > feat(rdx + 0x78);  // unordered -> false
    if (!above) goto L461;
    xd[1] = 2.2514;
    above = xd[1] > feat(rdx + 0x60);  // unordered -> false
    if (!above) goto L388;
    xd[1] = 3.5;
    above = xd[1] > feat(rdx + 0x120);  // unordered -> false
    if (above) goto L384;
    xd[1] = 0.413465;
    rdx = rcx + 0x10;
    store(rcx + 0x8, 0x1ull, 8);
    store(rcx + 0x0, rdx, 8);
    above = xd[1] > xd[0];  // unordered -> false
    if (!above) goto L474;
    goto L408;
L417:
    xd[0] = 5.5;
    above = xd[0] > feat(rdx + 0x120);  // unordered -> false
    if (!above) goto L478;
    xd[0] = 0.136787;
    above = xd[0] > feat(rdx + 0x148);  // unordered -> false
    if (above) goto L453;
L461:
    rdx = rax + 0x10;
    store(rax + 0x0, rdx, 8);
L394:
    store(rax + 0x8, 0x1ull, 8);
L474:
    rcx = 0x34ull;
    store(rax + 0x10, rcx, 2);
    goto done;
L421:
    xd[0] = 986.0;
    above = xd[0] > feat(rdx + 0x0);  // unordered -> false
    rdx = rcx + 0x10;
    store(rcx + 0x8, 0x1ull, 8);
    store(rcx + 0x0, rdx, 8);
    if (!above) goto L408;
    goto L373;
L449:
    xd[1] = 3296.27;
    above = xd[1] > feat(rdx + 0x198);  // unordered -> false
    if (above) goto L388;
    goto L370;
L478:
    xd[0] = 8.5;
    above = xd[0] > feat(rdx + 0x18);  // unordered -> false
    if (!above) goto L496;
    xd[0] = 6.5;
    above = xd[0] > feat(rdx + 0x80);  // unordered -> false
    if (above) goto L388;
    xd[0] = 0.313304;
    above = xd[0] > feat(rdx + 0x150);  // unordered -> false
    rdx = rcx + 0x10;
    store(rcx + 0x8, 0x1ull, 8);
    store(rcx + 0x0, rdx, 8);
    if (!above) goto L506;
    r11 = 0x36ull;
    store(rcx + 0x10, r11, 2);
    goto done;
L496:
    xd[0] = 0.129347;
    above = xd[0] > feat(rdx + 0x160);  // unordered -> false
    rdx = rcx + 0x10;
    store(rcx + 0x8, 0x1ull, 8);
    store(rcx + 0x0, rdx, 8);
    if (!above) goto L515;
    r9 = 0x38ull;
    store(rcx + 0x10, r9, 2);
    goto done;
L425:
    xd[0] = 1.24186;
    above = xd[0] > feat(rdx + 0x50);  // unordered -> false
    rdx = rcx + 0x10;
    store(rcx + 0x0, rdx, 8);
    if (!above) goto L396;
    rdx = 0x38ull;
    store(rcx + 0x8, 0x1ull, 8);
    store(rcx + 0x10, rdx, 2);
    goto done;
L515:
    r8 = 0x39ull;
    store(rcx + 0x10, r8, 2);
    goto done;
L506:
    r10 = 0x37ull;
    store(rcx + 0x10, r10, 2);
    goto done;
    goto L388;
done:
    uint64_t len;
    std::memcpy(&len, res + 8, 8);
    return std::string(reinterpret_cast<const char*>(res + 16), static_cast<size_t>(len));
}

}  // namespace isvg
