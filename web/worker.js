// Runs the engine off the main thread.
// Messages in:  {id, name, bytes: ArrayBuffer, maxPixels, colors}  vectorize an image
//                                         (colors > 0: reduce the image to that many colors first)
//               {id, kind: "export", flags: Uint8Array, colors: Uint32Array}
//                                         re-export the last result: flags[region] = 1 removes the
//                                         shape, colors[region] = 0x01RRGGBB recolours it (0 = keep)
// Messages out: {id, type: "progress", stage, fraction} | {id, type: "done", svg, info, width, height}
//               | {id, type: "exported", svg} | {id, type: "error", message}
import createEngine from "./engine/engine.mjs";

const enginePromise = createEngine();

// PNG, BMP, GIF and TGA are decoded by stb_image inside the engine so the pixels are exactly what
// the engine expects (the browser canvas would round semi-transparent colours). Everything else
// (JPEG, WebP) is decoded by the browser without colour management.
async function decodeWithBrowser(bytes) {
  const bitmap = await createImageBitmap(new Blob([bytes]), {
    colorSpaceConversion: "none",
    premultiplyAlpha: "none",
  });
  const canvas = new OffscreenCanvas(bitmap.width, bitmap.height);
  const ctx = canvas.getContext("2d", { colorSpace: "srgb" });
  ctx.drawImage(bitmap, 0, 0);
  const rgba = ctx.getImageData(0, 0, bitmap.width, bitmap.height).data;
  for (let i = 0; i < rgba.length; i += 4) {
    const r = rgba[i];
    rgba[i] = rgba[i + 2];
    rgba[i + 2] = r;
  }
  return { width: bitmap.width, height: bitmap.height, bgra: rgba };
}

function rgbToLab(r, g, b) {
  const lin = (v) => {
    v /= 255;
    return v <= 0.04045 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4;
  };
  r = lin(r);
  g = lin(g);
  b = lin(b);
  const x = (0.4124 * r + 0.3576 * g + 0.1805 * b) / 0.95047;
  const y = 0.2126 * r + 0.7152 * g + 0.0722 * b;
  const z = (0.0193 * r + 0.1192 * g + 0.9505 * b) / 1.08883;
  const f = (t) => (t > 216 / 24389 ? Math.cbrt(t) : (24389 / 27 * t + 16) / 116);
  return [116 * f(y) - 16, 500 * (f(x) - f(y)), 200 * (f(y) - f(z))];
}

function dist2(a, b) {
  return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2;
}

// Same block averaging as the engine's rescaleImage (f x f blocks, f = ceil(sqrt(w*h / maxPixels))),
// on straight B,G,R,A bytes: colors are averaged weighted by alpha, alpha is averaged plainly.
function shrinkImage(bgra, w, h, maxPixels) {
  const f = Math.ceil(Math.sqrt((w * h) / maxPixels));
  const ow = Math.ceil(w / f), oh = Math.ceil(h / f);
  const out = new Uint8Array(ow * oh * 4);
  for (let oy = 0; oy < oh; oy++) {
    const y1 = Math.min(h, (oy + 1) * f);
    for (let ox = 0; ox < ow; ox++) {
      const x1 = Math.min(w, (ox + 1) * f);
      let sa = 0, sb = 0, sg = 0, sr = 0, n = 0;
      for (let y = oy * f; y < y1; y++)
        for (let x = ox * f; x < x1; x++) {
          const i = 4 * (y * w + x), a = bgra[i + 3];
          sb += a * bgra[i];
          sg += a * bgra[i + 1];
          sr += a * bgra[i + 2];
          sa += a;
          n++;
        }
      const o = 4 * (oy * ow + ox);
      if (sa) {
        out[o] = Math.round(sb / sa);
        out[o + 1] = Math.round(sg / sa);
        out[o + 2] = Math.round(sr / sa);
      }
      out[o + 3] = Math.round(sa / n);
    }
  }
  return { bgra: out, w: ow, h: oh };
}

// Reduces the image to `count` colors before the engine sees it. Left alone, a gradient becomes many
// bands of slightly different colors, each its own shape. The colors come from k-means over a
// 32x32x32 RGB histogram of the visible pixels (weighted by alpha, compared in CIELAB); every pixel
// then takes its cluster's color (or, on an edge, a mix of the two colors either side of it) and
// keeps its alpha, so edges against transparency stay smooth.
// bgra: straight (not premultiplied) B,G,R,A bytes of a `width` pixels wide image, changed in place.
function reduceColors(bgra, width, count) {
  const hist = new Float64Array(32768 * 4);  // per cell: weight, weighted sums of B, G, R
  for (let i = 0; i < bgra.length; i += 4) {
    const a = bgra[i + 3];
    if (a === 0) continue;
    const h = 4 * (((bgra[i + 2] >> 3) << 10) | ((bgra[i + 1] >> 3) << 5) | (bgra[i] >> 3));
    hist[h] += a;
    hist[h + 1] += a * bgra[i];
    hist[h + 2] += a * bgra[i + 1];
    hist[h + 3] += a * bgra[i + 2];
  }
  const cells = [];
  for (let c = 0; c < 32768; c++) {
    const w = hist[4 * c];
    if (!w) continue;
    const rgb = [hist[4 * c + 3] / w, hist[4 * c + 2] / w, hist[4 * c + 1] / w];
    cells.push({ c, w, rgb, lab: rgbToLab(...rgb), k: 0 });
  }
  if (!cells.length) return;

  // Deterministic k-means++ style seeds: the heaviest cell, then the cell with the largest
  // weight x squared distance to its nearest seed.
  const centers = [cells.reduce((a, b) => (b.w > a.w ? b : a)).lab];
  while (centers.length < count) {
    let best = null, bestScore = 0;
    for (const cell of cells) {
      const score = cell.w * Math.min(...centers.map((ctr) => dist2(cell.lab, ctr)));
      if (score > bestScore) {
        bestScore = score;
        best = cell;
      }
    }
    if (!best) break;  // fewer distinct colors than asked for
    centers.push(best.lab);
  }

  let palette = [];
  for (let iter = 0; iter < 40; iter++) {
    let moved = false;
    for (const cell of cells) {
      let k = 0;
      for (let j = 1; j < centers.length; j++) if (dist2(cell.lab, centers[j]) < dist2(cell.lab, centers[k])) k = j;
      if (k !== cell.k) moved = true;
      cell.k = k;
    }
    const sums = centers.map(() => [0, 0, 0, 0]);
    for (const cell of cells) {
      const s = sums[cell.k];
      s[0] += cell.w;
      for (let ch = 0; ch < 3; ch++) s[ch + 1] += cell.w * cell.rgb[ch];
    }
    palette = sums.map((s, j) => (s[0] ? [s[1] / s[0], s[2] / s[0], s[3] / s[0]] : null));
    palette.forEach((rgb, j) => {
      if (rgb) centers[j] = rgbToLab(...rgb);
    });
    if (!moved && iter > 0) break;
  }

  // Anti-aliased edge pixels keep their mix of the two colors they sit between. Snapped to the
  // nearest color (often an unrelated third one: white + purple fringe -> grey), an edge turns into a
  // noisy staircase with a 1 px outline, which the engine then traces as wobbly curves. A cell is a
  // mix of colors i and j when it is not close to any one color and the i-j blend explains it
  // better; a pixel only keeps the mix when both i and j are within 2 px of it, so a wide gradient
  // between two colors still snaps as before.
  const rgb8 = palette.map((p) => p && p.map((v) => Math.round(v)));
  const labs = rgb8.map((p) => p && rgbToLab(...p));
  const lut = new Int16Array(32768).fill(-1);
  const pairs = new Int16Array(32768).fill(-1);  // i * 256 + j
  for (const cell of cells) {
    lut[cell.c] = cell.k;
    let best = dist2(cell.lab, labs[cell.k]);
    if (best < AUTO_JOIN * AUTO_JOIN) continue;
    for (let i = 0; i < rgb8.length; i++)
      for (let j = i + 1; j < rgb8.length; j++) {
        if (!rgb8[i] || !rgb8[j]) continue;
        const t = mixFactor(cell.rgb, rgb8[i], rgb8[j]);
        if (t <= 0 || t >= 1) continue;
        const e = dist2(cell.lab, rgbToLab(...mix(rgb8[i], rgb8[j], t)));
        if (e < best) {
          best = e;
          pairs[cell.c] = i * 256 + j;
        }
      }
  }
  const cellOf = (i) => ((bgra[i + 2] >> 3) << 10) | ((bgra[i + 1] >> 3) << 5) | (bgra[i] >> 3);
  const kind = new Int16Array(bgra.length / 4);  // palette index, -1 transparent, -2 possible mix
  for (let p = 0; p < kind.length; p++) {
    const i = 4 * p;
    kind[p] = bgra[i + 3] === 0 ? -1 : pairs[cellOf(i)] >= 0 ? -2 : lut[cellOf(i)];
  }
  const height = kind.length / width;
  for (let p = 0; p < kind.length; p++) {
    if (kind[p] === -1) continue;
    const i = 4 * p, c = cellOf(i);
    let out = rgb8[lut[c]];
    if (kind[p] === -2) {
      const a = pairs[c] >> 8, b = pairs[c] & 255;
      const x = p % width, y = (p - x) / width;
      let hasA = false, hasB = false;
      for (let yy = Math.max(0, y - 2); yy <= Math.min(height - 1, y + 2); yy++)
        for (let xx = Math.max(0, x - 2); xx <= Math.min(width - 1, x + 2); xx++) {
          const k = kind[yy * width + xx];
          hasA ||= k === a;
          hasB ||= k === b;
        }
      if (hasA && hasB) {
        const t = Math.min(1, Math.max(0, mixFactor([bgra[i + 2], bgra[i + 1], bgra[i]], rgb8[a], rgb8[b])));
        out = mix(rgb8[a], rgb8[b], t).map((v) => Math.round(v));
      }
    }
    bgra[i] = out[2];
    bgra[i + 1] = out[1];
    bgra[i + 2] = out[0];
  }
  return rgb8.filter(Boolean);
}

// Where c lies along the line from a to b (0 at a, 1 at b), in RGB.
function mixFactor(c, a, b) {
  const d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
  const len2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
  return len2 ? ((c[0] - a[0]) * d[0] + (c[1] - a[1]) * d[1] + (c[2] - a[2]) * d[2]) / len2 : 0;
}

function mix(a, b, t) {
  return [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];
}

// "Auto" colors: a logo or flat illustration has a handful of dominant colors plus a thin fringe of
// anti-aliasing and compression noise, which would otherwise come out as many near-identical colors.
// Greedily groups the 32x32x32 histogram cells (heaviest first) around centers, joining a cell to the
// first center within AUTO_JOIN (CIELAB distance). If a few centers each hold a real share of the
// pixels and together hold nearly all of them, that is the color count; photos and gradients, which
// spread over many colors, return 0 (leave the image alone).
const AUTO_JOIN = 14;
const AUTO_MAX_COLORS = 12;

function autoColorCount(bgra) {
  const hist = new Float64Array(32768 * 4);
  let total = 0;
  for (let i = 0; i < bgra.length; i += 4) {
    const a = bgra[i + 3];
    if (a === 0) continue;
    const h = 4 * (((bgra[i + 2] >> 3) << 10) | ((bgra[i + 1] >> 3) << 5) | (bgra[i] >> 3));
    hist[h] += a;
    hist[h + 1] += a * bgra[i];
    hist[h + 2] += a * bgra[i + 1];
    hist[h + 3] += a * bgra[i + 2];
    total += a;
  }
  if (!total) return 0;
  const cells = [];
  for (let c = 0; c < 32768; c++) {
    const w = hist[4 * c];
    if (w) cells.push({ w, lab: rgbToLab(hist[4 * c + 3] / w, hist[4 * c + 2] / w, hist[4 * c + 1] / w) });
  }
  cells.sort((a, b) => b.w - a.w);
  const centers = [];
  for (const cell of cells) {
    const near = centers.find((c) => dist2(cell.lab, c.lab) < AUTO_JOIN * AUTO_JOIN);
    if (near) near.w += cell.w;
    else centers.push({ lab: cell.lab, w: cell.w });
  }
  const major = centers.filter((c) => c.w >= 0.015 * total);
  const covered = major.reduce((sum, c) => sum + c.w, 0);
  return major.length > 0 && major.length <= AUTO_MAX_COLORS && covered >= 0.9 * total ? major.length : 0;
}

// The engine averages the pixels inside each shape, which drifts a shape's color by a step or two
// (and blends anti-aliased edges into in-between shades). With a color target, every shape snaps
// back to the nearest chosen color. snap[region] = 0x01RRGGBB, used again when re-exporting.
let snap = null;

function snapSvg(svg, info, palette) {
  const labs = palette.map((p) => rgbToLab(...p));
  snap = new Uint32Array(info.shapes);
  let k = 0;
  return svg.replace(/<path fill="#([0-9a-f]{6})"/g, (_, hex) => {
    const n = parseInt(hex, 16);
    const lab = rgbToLab((n >> 16) & 255, (n >> 8) & 255, n & 255);
    let best = 0;
    for (let j = 1; j < labs.length; j++) if (dist2(lab, labs[j]) < dist2(lab, labs[best])) best = j;
    const [r, g, b] = palette[best];
    const out = (r << 16) | (g << 8) | b;
    snap[info.pathRegions[k++]] = 0x1000000 | out;
    return `<path fill="#${out.toString(16).padStart(6, "0")}"`;
  });
}

async function exportWithEdits(id, flags, colors) {
  const M = await enginePromise;
  const ptr = M._isvg_alloc(flags.length);
  M.HEAPU8.set(flags, ptr);
  if (snap) for (let r = 0; r < colors.length; r++) if (!colors[r]) colors[r] = snap[r] || 0;
  const cptr = M._isvg_alloc(colors.length * 4);
  new Uint32Array(M.HEAPU8.buffer, cptr, colors.length).set(colors);
  const svg = M.UTF8ToString(M._isvg_export(ptr, cptr, flags.length));
  M._isvg_free(ptr);
  M._isvg_free(cptr);
  self.postMessage({ id, type: "exported", svg });
}

self.onmessage = async (event) => {
  if (event.data.kind === "export") {
    exportWithEdits(event.data.id, event.data.flags, event.data.colors).catch((err) =>
      self.postMessage({ id: event.data.id, type: "error", message: String(err) }),
    );
    return;
  }
  const { id, name, bytes, maxPixels = 0, colors = 0 } = event.data;
  try {
    const M = await enginePromise;
    M.onProgress = (stage, fraction) => self.postMessage({ id, type: "progress", stage, fraction });

    const data = new Uint8Array(bytes);
    const src = M._isvg_alloc(data.length);
    M.HEAPU8.set(data, src);
    const wh = M._isvg_alloc(8);
    let px = M._isvg_decode(src, data.length, wh, wh + 4);
    let width, height;
    if (px) {
      width = M.HEAP32[wh >> 2];
      height = M.HEAP32[(wh >> 2) + 1];
    } else {
      const img = await decodeWithBrowser(data);
      width = img.width;
      height = img.height;
      px = M._isvg_alloc(img.bgra.length);
      M.HEAPU8.set(img.bgra, px);
    }
    M._isvg_free(src);
    M._isvg_free(wh);
    let vecW = width, vecH = height, vecMax = maxPixels, palette = null;
    snap = null;
    const target = colors > 0 ? colors : autoColorCount(M.HEAPU8.subarray(px, px + width * height * 4));
    if (target > 0) {
      // Shrink first, then reduce: the engine's own shrinking blends neighbouring colors into new
      // in-between ones, which would undo the color target.
      if (maxPixels > 0 && width * height > maxPixels) {
        const small = shrinkImage(M.HEAPU8.slice(px, px + width * height * 4), width, height, maxPixels);
        M._isvg_free(px);
        px = M._isvg_alloc(small.bgra.length);
        M.HEAPU8.set(small.bgra, px);
        vecW = small.w;
        vecH = small.h;
        vecMax = 0;
      }
      palette = reduceColors(M.HEAPU8.subarray(px, px + vecW * vecH * 4), vecW, target);
    }

    const dot = name.lastIndexOf(".");
    const ext = dot >= 0 ? name.slice(dot).toLowerCase() : "";
    const extLen = M.lengthBytesUTF8(ext) + 1;
    const extPtr = M._isvg_alloc(extLen);
    M.stringToUTF8(ext, extPtr, extLen);

    let svg = M.UTF8ToString(M._isvg_vectorize(px, vecW, vecH, extPtr, vecMax));
    const info = JSON.parse(M.UTF8ToString(M._isvg_info()));
    if (palette && palette.length) svg = snapSvg(svg, info, palette);
    M._isvg_free(px);
    M._isvg_free(extPtr);
    self.postMessage({ id, type: "done", svg, info, width, height });
  } catch (err) {
    self.postMessage({ id, type: "error", message: String(err && err.message ? err.message : err) });
  }
};
