const STAGES = {
  classify: { label: "Analyzing image", from: 0.0, to: 0.08 },
  palettes: { label: "Finding colors", from: 0.08, to: 0.15 },
  segment: { label: "Segmenting", from: 0.15, to: 0.4 },
  trace: { label: "Tracing boundaries", from: 0.4, to: 0.45 },
  smooth: { label: "Smoothing contours", from: 0.45, to: 0.95 },
  bezier: { label: "Fitting curves", from: 0.95, to: 0.98 },
  export: { label: "Writing SVG", from: 0.98, to: 1.0 },
};

const $ = (id) => document.getElementById(id);
const drop = $("drop");
const fileInput = $("file");
const work = $("work");
const panes = $("panes");
const result = $("result");
const downloadBtn = $("download");

let worker = null;
let jobId = 0;
let currentSvg = null;
let currentName = "image";
let originalUrl = null;


function startWorker() {
  if (worker) worker.terminate();
  worker = new Worker(new URL("./worker.js", import.meta.url), { type: "module" });
  worker.onmessage = onWorkerMessage;
}

function setProgress(fraction, text) {
  $("bar-fill").style.width = `${Math.round(fraction * 100)}%`;
  $("stage").textContent = text;
}

function onWorkerMessage(event) {
  const msg = event.data;
  if (msg.type === "exported") {
    if (msg.id === exportJob) saveSvg(msg.svg);
    return;
  }
  if (msg.id !== jobId) return;
  if (msg.type === "progress") {
    const s = STAGES[msg.stage] || { label: msg.stage, from: 0, to: 1 };
    setProgress(s.from + (s.to - s.from) * msg.fraction, `${s.label}…`);
  } else if (msg.type === "done") {
    converting = false;
    currentSvg = msg.svg;
    window.__lastSvg = msg.svg;  // read by tests/browser_check.mjs
    result.innerHTML = msg.svg;
    const svgEl = result.querySelector("svg");
    if (svgEl) {
      svgEl.removeAttribute("width");
      svgEl.removeAttribute("height");
    }
    const i = msg.info;
    setupEditing(i);
    if (i.width !== msg.width || i.height !== msg.height)
      $("fdims").textContent = `${msg.width} × ${msg.height}px, vectorized at ${i.width} × ${i.height}px`;
    const kind = { 0: "Photo", 1: "Logo, anti-aliased", 2: "Logo, hard edges" }[i.classification[0]] || "";
    const quality = ["high", "medium", "low"][i.classification[1]] || "";
    $("summary").textContent = `· ${kind}, ${quality} detail · ${i.shapes} shapes · ${i.curves} curves`;
    setProgress(1, `Done in ${((performance.now() - started) / 1000).toFixed(1)} s`);
    $("progress").classList.add("done");
    downloadBtn.disabled = false;
  } else if (msg.type === "error") {
    converting = false;
    setProgress(0, `Failed: ${msg.message}`);
    $("progress").classList.add("error");
  }
}

let started = 0;

// Image size limits: shrinkImageAsNeeded (the "Small" choice) and the memory cap of capImageSize,
// which is 20 megapixels on any machine with enough memory.
const SMALL_MAX_PIXELS = 1600000;
const LARGE_MAX_PIXELS = 20000000;

function imageSize(file) {
  return new Promise((resolve) => {
    const img = new Image();
    img.onload = () => resolve({ w: img.naturalWidth, h: img.naturalHeight });
    img.onerror = () => resolve(null);
    img.src = URL.createObjectURL(file);
  });
}

function askSize(pixels) {
  return new Promise((resolve) => {
    const box = $("sizechoice");
    $("megapixels").textContent = (pixels / 1e6).toFixed(1);
    box.hidden = false;
    $("progress").hidden = true;
    const done = (value) => {
      box.hidden = true;
      $("progress").hidden = false;
      resolve(value);
    };
    $("small").onclick = () => done(SMALL_MAX_PIXELS);
    $("large").onclick = () => done(LARGE_MAX_PIXELS);
  });
}

let currentFile = null;
let currentMaxPixels = 0;
let converting = false;

// maxPixels: the size choice already made for this file (re-converting with another color count).
async function handleFile(file, maxPixels = 0) {
  if (!file) return;
  if (!/^image\//.test(file.type) && !/\.(png|jpe?g|bmp|gif|webp|tga)$/i.test(file.name)) {
    alert("Please choose an image file.");
    return;
  }
  currentName = file.name.replace(/\.[^.]+$/, "") || "image";
  currentSvg = null;
  edits = null;
  $("editpanel").hidden = true;
  downloadBtn.disabled = true;
  result.innerHTML = "";
  $("summary").textContent = "";
  $("progress").classList.remove("done", "error");
  if (originalUrl) URL.revokeObjectURL(originalUrl);
  originalUrl = URL.createObjectURL(file);
  const img = $("original");
  img.src = originalUrl;
  img.onload = () => {
    $("fdims").textContent = `${img.naturalWidth} × ${img.naturalHeight}px`;
  };
  $("fname").textContent = file.name;
  drop.hidden = true;
  work.hidden = false;
  setProgress(0, "Loading engine…");

  if (!maxPixels) {
    const size = await imageSize(file);
    maxPixels = LARGE_MAX_PIXELS;
    if (size && size.w * size.h > SMALL_MAX_PIXELS) maxPixels = await askSize(size.w * size.h);
  }
  currentFile = file;
  currentMaxPixels = maxPixels;

  const bytes = await file.arrayBuffer();
  const colors = Number($("numcolors").value);
  jobId += 1;
  started = performance.now();
  // A conversion still running cannot be interrupted inside the engine; restart the worker instead.
  if (!worker || converting) startWorker();
  converting = true;
  worker.postMessage({ id: jobId, name: file.name, bytes, maxPixels, colors }, [bytes]);
}

// The Colors dropdown exists twice (upload box, edit panel); keep them the same.
for (const [from, to] of [["numcolors", "numcolors2"], ["numcolors2", "numcolors"]])
  $(from).addEventListener("change", () => ($(to).value = $(from).value));

$("reconvert").addEventListener("click", () => {
  if (!currentFile) return;
  if (hasEdits() && !confirm("Reconverting discards your edits. Continue?")) return;
  handleFile(currentFile, currentMaxPixels);
});

function reset() {
  // A running job cannot be interrupted inside the engine; restart the worker instead.
  if (converting && worker) {
    startWorker();
    converting = false;
  }
  work.hidden = true;
  drop.hidden = false;
  fileInput.value = "";
  drop.focus();
}

$("browse").addEventListener("click", () => fileInput.click());
drop.addEventListener("click", (e) => {
  if (e.target.closest(".options")) return;  // the Colors dropdown
  if (e.target === drop || e.target.closest(".drop-inner")) fileInput.click();
});
drop.addEventListener("keydown", (e) => {
  if (e.target !== drop) return;  // keys in the Colors dropdown
  if (e.key === "Enter" || e.key === " ") fileInput.click();
});
fileInput.addEventListener("change", () => handleFile(fileInput.files[0]));
["dragenter", "dragover"].forEach((t) =>
  drop.addEventListener(t, (e) => {
    e.preventDefault();
    drop.classList.add("over");
  }),
);
["dragleave", "drop"].forEach((t) =>
  drop.addEventListener(t, (e) => {
    e.preventDefault();
    drop.classList.remove("over");
  }),
);
drop.addEventListener("drop", (e) => handleFile(e.dataTransfer.files[0]));
document.addEventListener("paste", (e) => {
  const item = [...(e.clipboardData?.items || [])].find((i) => i.type.startsWith("image/"));
  if (item) handleFile(item.getAsFile());
});

$("another").addEventListener("click", reset);

// Shapes that touch share their edges exactly, but a renderer anti-aliases each <path> on its own:
// along a shared edge both shapes only partly cover the edge pixels and the background shows
// through as a hairline. Once the shapes have the same color (e.g. the bands of a gradient merged
// into one color) those hairlines look like cuts. Written as one <path>, a shared edge is drawn
// twice in opposite directions and cancels out, so the seam disappears.
const PATH_LINE = /(\r?\n)<path fill="(#[0-9a-f]{6})" opacity="([\d.]+)" d="([^"]*)" \/>/g;

function mergeSameFillPaths(svg) {
  const groups = new Map();
  for (const [, , fill, opacity, d] of svg.matchAll(PATH_LINE)) {
    const key = `${fill} ${opacity}`;
    if (!groups.has(key)) groups.set(key, []);
    groups.get(key).push(d);
  }
  const written = new Set();
  return svg.replace(PATH_LINE, (line, nl, fill, opacity) => {
    const key = `${fill} ${opacity}`;
    if (written.has(key)) return "";
    written.add(key);
    return `${nl}<path fill="${fill}" opacity="${opacity}" d="${groups.get(key).join("")}" />`;
  });
}

function saveSvg(svg) {
  svg = mergeSameFillPaths(svg);
  window.__savedSvg = svg;  // read by tests/remove_bg_check.mjs
  // The SVG text is ASCII with CRLF line endings, exactly as the engine's SVG writer produces it.
  const blob = new Blob([svg], { type: "image/svg+xml" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = `${currentName}.svg`;
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

downloadBtn.addEventListener("click", () => {
  if (!currentSvg || !edits) return;
  if (!hasEdits()) {
    saveSvg(currentSvg);
    return;
  }
  // Edited: the engine writes the SVG (same writer and number format as the unedited result).
  exportJob = ++jobId;
  const colors = new Uint32Array(edits.n);
  for (let r = 0; r < edits.n; r++)
    if (edits.color[r] && edits.color[r] !== edits.orig[r]) colors[r] = 0x1000000 | parseInt(edits.color[r].slice(1), 16);
  worker.postMessage({ id: exportJob, kind: "export", flags: edits.removed.slice(), colors });
});

// --- Editing: remove background, recolor, merge colors, paint shapes ------------------------
//
// Every edit is a list of changes {region, key: "color" | "removed", before, after} so a single
// undo/redo history covers all tools. Removing a shape works like the engine's
// the export settings (the shape is left out of the export, the hole it cut into
// its parent stays). Color edits only change fills; opacity and geometry stay as vectorized.

let edits = null;  // {n, orig[], color[], removed: Uint8Array, undo[], redo[]}
let tool = "select";
let exportJob = 0;
let previewBase = null;  // colors before a live color-picker preview

const pathsOf = () => [...result.querySelectorAll("path")];
const pathsByRegion = new Map();

function hasEdits() {
  return edits.removed.some((b) => b) || edits.color.some((c, r) => c && c !== edits.orig[r]);
}

function setupEditing(info) {
  const n = info.shapes;
  edits = { n, orig: new Array(n).fill(null), color: new Array(n).fill(null), removed: new Uint8Array(n), undo: [], redo: [] };
  pathsByRegion.clear();
  pathsOf().forEach((p, k) => {
    const r = info.pathRegions[k];
    p.dataset.region = r;
    pathsByRegion.set(r, p);
    edits.orig[r] = edits.color[r] = p.getAttribute("fill").toLowerCase();
  });
  $("editpanel").hidden = false;
  setTool("select");
  refresh();
}

function applyChanges(changes, useAfter) {
  for (const c of changes) edits[c.key][c.region] = useAfter ? c.after : c.before;
}

function commit(changes) {
  changes = changes.filter((c) => c.before !== c.after);
  if (!changes.length) return;
  applyChanges(changes, true);
  edits.undo.push(changes);
  edits.redo = [];
  refresh();
}

function undo() {
  const ch = edits?.undo.pop();
  if (!ch) return;
  applyChanges(ch, false);
  edits.redo.push(ch);
  refresh();
}

function redo() {
  const ch = edits?.redo.pop();
  if (!ch) return;
  applyChanges(ch, true);
  edits.undo.push(ch);
  refresh();
}

// Colors in use, most shapes first (removed shapes are not counted).
function palette() {
  const groups = new Map();
  for (const [r] of pathsByRegion) {
    const c = edits.color[r];
    if (!groups.has(c)) groups.set(c, { color: c, regions: [], visible: 0 });
    const g = groups.get(c);
    g.regions.push(r);
    if (!edits.removed[r]) g.visible++;
  }
  return [...groups.values()].filter((g) => g.visible > 0).sort((a, b) => b.visible - a.visible);
}

function refresh() {
  for (const [r, p] of pathsByRegion) {
    p.setAttribute("fill", edits.color[r]);
    p.classList.toggle("removed", edits.removed[r] === 1);
  }
  renderUnderlay();
  $("undo").disabled = edits.undo.length === 0;
  $("redo").disabled = edits.redo.length === 0;
  $("resetedits").disabled = !hasEdits();
  renderPalette();
  updateHint();
}

// The preview keeps one <path> per shape so shapes stay clickable. To hide the seams between
// same-color shapes (see mergeSameFillPaths), an underlay below them draws every color used by
// more than one opaque shape as a single path; the seams then show that same color instead of the
// background.
const SVG_NS = "http://www.w3.org/2000/svg";
let underlayFrame = 0;

function renderUnderlay() {
  const svgEl = result.querySelector("svg");
  if (!svgEl) return;
  let layer = svgEl.querySelector("g.underlay");
  if (!layer) {
    layer = document.createElementNS(SVG_NS, "g");
    layer.setAttribute("class", "underlay");
    svgEl.prepend(layer);
  }
  const byFill = new Map();
  for (const [r, p] of pathsByRegion) {
    if (edits.removed[r] || p.getAttribute("opacity") !== "1.00") continue;
    const fill = p.getAttribute("fill");
    if (!byFill.has(fill)) byFill.set(fill, []);
    byFill.get(fill).push(p.getAttribute("d"));
  }
  const merged = [];
  for (const [fill, ds] of byFill) {
    if (ds.length < 2) continue;
    const path = document.createElementNS(SVG_NS, "path");
    path.setAttribute("fill", fill);
    path.setAttribute("d", ds.join(""));
    merged.push(path);
  }
  layer.replaceChildren(...merged);
}

function scheduleUnderlay() {
  if (underlayFrame) return;
  underlayFrame = requestAnimationFrame(() => {
    underlayFrame = 0;
    renderUnderlay();
  });
}

function updateHint() {
  const removed = edits.removed.reduce((a, b) => a + b, 0);
  const text = {
    select: "",
    remove: "Click a shape to remove it, click again to bring it back.",
    fill: "Click a shape to fill it with the color around it.",
    paint: "Click a shape to paint it.",
  }[tool];
  $("edithint").textContent = [text, removed ? `${removed} removed.` : ""].filter(Boolean).join(" ");
  $("summary").hidden = tool !== "select";
}

const checked = new Set();

function renderPalette() {
  const list = $("colorlist");
  const groups = palette();
  for (const c of [...checked]) if (!groups.some((g) => g.color === c)) checked.delete(c);
  $("colorcount").textContent = `(${groups.length})`;
  list.replaceChildren(
    ...groups.map((g) => {
      const li = document.createElement("li");
      li.dataset.color = g.color;
      li.innerHTML = `
        <input type="checkbox" class="pick" aria-label="Select ${g.color} for merging">
        <label class="swatch" title="Change this color everywhere">
          <input type="color" value="${g.color}">
        </label>
        <button type="button" class="name" title="Paint with this color">${g.color}</button>
        <span class="count">${g.visible} shape${g.visible === 1 ? "" : "s"}</span>`;
      li.querySelector(".swatch").style.background = g.color;
      const box = li.querySelector(".pick");
      box.checked = checked.has(g.color);
      box.addEventListener("change", () => {
        box.checked ? checked.add(g.color) : checked.delete(g.color);
        $("merge").disabled = checked.size < 2;
      });
      const picker = li.querySelector('input[type="color"]');
      picker.addEventListener("input", () => previewRecolor(g, picker.value));
      picker.addEventListener("change", () => finishRecolor(g, picker.value));
      li.querySelector(".name").addEventListener("click", () => {
        setPaintColor(g.color);
        setTool("paint");
      });
      li.addEventListener("mouseenter", () => highlightColor(g.color));
      li.addEventListener("mouseleave", () => highlightColor(null));
      return li;
    }),
  );
  $("merge").disabled = checked.size < 2;
}

function highlightColor(color) {
  result.classList.toggle("spot", color !== null);
  for (const [r, p] of pathsByRegion) p.classList.toggle("spotlit", color !== null && edits.color[r] === color);
}

function previewRecolor(group, value) {
  if (!previewBase) previewBase = { color: group.color };
  for (const r of group.regions) pathsByRegion.get(r).setAttribute("fill", value);
  scheduleUnderlay();
}

function finishRecolor(group, value) {
  previewBase = null;
  value = value.toLowerCase();
  commit(group.regions.map((r) => ({ region: r, key: "color", before: edits.color[r], after: value })));
}

$("merge").addEventListener("click", () => {
  const groups = palette().filter((g) => checked.has(g.color));
  if (groups.length < 2) return;
  const target = groups[0].color;  // palette is sorted: the color used by the most shapes wins
  const changes = [];
  for (const g of groups.slice(1))
    for (const r of g.regions) changes.push({ region: r, key: "color", before: edits.color[r], after: target });
  checked.clear();
  commit(changes);
});

$("resetedits").addEventListener("click", () => {
  const changes = [];
  for (let r = 0; r < edits.n; r++) {
    if (edits.color[r] !== edits.orig[r]) changes.push({ region: r, key: "color", before: edits.color[r], after: edits.orig[r] });
    if (edits.removed[r]) changes.push({ region: r, key: "removed", before: 1, after: 0 });
  }
  commit(changes);
});

// The color that most of a shape's outline borders: points just outside the outline are looked up
// in the other visible shapes (topmost first). Returns null when the background borders it most,
// undefined when it cannot be told.
function surroundingColor(region) {
  const own = pathsByRegion.get(region);
  const others = [...pathsByRegion].filter(([r]) => r !== region && !edits.removed[r]).reverse();
  const len = own.getTotalLength();
  if (!(len > 0)) return undefined;
  const offset = 1.5;
  const samples = Math.min(200, Math.max(16, Math.round(len / 2)));
  const votes = new Map();
  for (let i = 0; i < samples; i++) {
    const at = (i / samples) * len, eps = Math.min(0.5, len / samples / 2);
    const p = own.getPointAtLength(at), a = own.getPointAtLength(Math.max(0, at - eps)), b = own.getPointAtLength(Math.min(len, at + eps));
    const tx = b.x - a.x, ty = b.y - a.y, tl = Math.hypot(tx, ty);
    if (!tl) continue;
    const nx = -ty / tl, ny = tx / tl;
    let q = null;
    for (const side of [1, -1]) {
      const pt = new DOMPoint(p.x + side * nx * offset, p.y + side * ny * offset);
      if (!own.isPointInFill(pt)) {
        q = pt;
        break;
      }
    }
    if (!q) continue;
    let key = null;
    for (const [r, path] of others) {
      if (path.isPointInFill(q)) {
        key = edits.color[r];
        break;
      }
    }
    votes.set(key, (votes.get(key) || 0) + 1);
  }
  let best, bestN = 0;
  for (const [key, n] of votes) if (n > bestN) [best, bestN] = [key, n];
  return bestN ? best : undefined;
}

function setPaintColor(hex) {
  $("paintcolor").value = hex;
  $("painthex").textContent = hex;
}
$("paintcolor").addEventListener("input", () => setPaintColor($("paintcolor").value.toLowerCase()));

function setTool(name) {
  tool = name;
  document.querySelectorAll(".tools .tool").forEach((b) => {
    b.classList.toggle("on", b.dataset.tool === name);
    b.setAttribute("aria-pressed", String(b.dataset.tool === name));
  });
  $("paintrow").hidden = name !== "paint";
  result.classList.toggle("editing", name !== "select");
  result.dataset.tool = name;
  if (edits) updateHint();
}
document.querySelectorAll(".tools .tool").forEach((b) => b.addEventListener("click", () => setTool(b.dataset.tool)));

$("undo").addEventListener("click", undo);
$("redo").addEventListener("click", redo);
result.addEventListener("click", (e) => {
  if (!edits || tool === "select") return;
  const path = e.target.closest("path");
  if (!path || path.dataset.region === undefined) return;
  const r = Number(path.dataset.region);
  if (tool === "remove") {
    commit([{ region: r, key: "removed", before: edits.removed[r], after: edits.removed[r] ^ 1 }]);
  } else if (tool === "fill") {
    const color = surroundingColor(r);
    if (color === undefined) return;
    const changes = [];
    if (color === null) {
      changes.push({ region: r, key: "removed", before: edits.removed[r], after: 1 });  // only background around it
    } else {
      changes.push({ region: r, key: "color", before: edits.color[r], after: color });
      if (edits.removed[r]) changes.push({ region: r, key: "removed", before: 1, after: 0 });
    }
    commit(changes);
  } else if (tool === "paint") {
    const changes = [{ region: r, key: "color", before: edits.color[r], after: $("paintcolor").value.toLowerCase() }];
    if (edits.removed[r]) changes.push({ region: r, key: "removed", before: 1, after: 0 });  // painting brings it back
    commit(changes);
  }
});
document.addEventListener("keydown", (e) => {
  if (!edits || !(e.ctrlKey || e.metaKey) || e.target.closest("input")) return;
  const k = e.key.toLowerCase();
  if (k === "z" && !e.shiftKey) {
    e.preventDefault();
    undo();
  } else if (k === "y" || (k === "z" && e.shiftKey)) {
    e.preventDefault();
    redo();
  }
});

document.querySelectorAll(".seg button").forEach((b) =>
  b.addEventListener("click", () => {
    document.querySelectorAll(".seg button").forEach((x) => x.classList.toggle("on", x === b));
    panes.className = `panes ${b.dataset.view}`;
    requestAnimationFrame(() => {
      const c = canvases.find((k) => contentOf(k));
      if (c) clampView(contentOf(c));
      applyView();
    });
  }),
);

// --- Zoom and pan, shared by the Original and Vector views ---------------------------------------
//
// One view {s, x, y} (scale, and translation in px) is written to both canvases as CSS variables,
// so whatever is zoomed or dragged in one pane is mirrored in the other.

const view = { s: 1, x: 0, y: 0 };
const MAX_ZOOM = 40;
const canvases = [...document.querySelectorAll(".canvas")];

function applyView() {
  for (const c of canvases) {
    c.style.setProperty("--zs", view.s);
    c.style.setProperty("--zx", `${view.x}px`);
    c.style.setProperty("--zy", `${view.y}px`);
  }
  $("zoomlevel").textContent = `${Math.round(view.s * 100)}%`;
}

// The content element of a canvas and its untransformed layout box (SVG elements have no offsetWidth).
function contentOf(canvas) {
  const el = canvas.querySelector(":scope > img, :scope > svg");
  if (!el) return null;
  const cs = getComputedStyle(el);
  const w = parseFloat(cs.width);
  const h = parseFloat(cs.height);
  if (!(w > 0 && h > 0)) return null;
  return { w, h, left: Math.max(0, (canvas.clientWidth - w) / 2), top: Math.max(0, (canvas.clientHeight - h) / 2) };
}

// Zoom to scale s keeping the point (px, py), in canvas coordinates, still under the cursor.
function zoomAt(canvas, s, px, py) {
  const el = contentOf(canvas);
  if (!el) return;
  s = Math.min(MAX_ZOOM, Math.max(1, s));
  const qx = px - el.left;
  const qy = py - el.top;
  view.x = qx - ((qx - view.x) / view.s) * s;
  view.y = qy - ((qy - view.y) / view.s) * s;
  view.s = s;
  clampView(el);
  applyView();
}

// Keep the zoomed content covering its pane (no panning off into empty space).
function clampView(el) {
  view.x = Math.min(0, Math.max(el.w * (1 - view.s), view.x));
  view.y = Math.min(0, Math.max(el.h * (1 - view.s), view.y));
}

function resetView() {
  view.s = 1;
  view.x = view.y = 0;
  applyView();
}

function zoomButton(factor) {
  const c = canvases.find((k) => contentOf(k));
  if (!c) return;
  zoomAt(c, view.s * factor, c.clientWidth / 2, c.clientHeight / 2);
}
$("zoomin").addEventListener("click", () => zoomButton(1.5));
$("zoomout").addEventListener("click", () => zoomButton(1 / 1.5));
$("zoomfit").addEventListener("click", resetView);

let drag = null;
let dragged = false;

for (const canvas of canvases) {
  canvas.addEventListener("wheel", (e) => {
    e.preventDefault();
    const r = canvas.getBoundingClientRect();
    zoomAt(canvas, view.s * Math.exp(-e.deltaY * 0.0015), e.clientX - r.left, e.clientY - r.top);
  }, { passive: false });

  canvas.addEventListener("pointerdown", (e) => {
    if (e.button !== 0) return;
    e.preventDefault();  // no native image/text drag
    drag = { id: e.pointerId, sx: e.clientX, sy: e.clientY, vx: view.x, vy: view.y };
    dragged = false;
  });
  canvas.addEventListener("pointermove", (e) => {
    if (!drag || drag.id !== e.pointerId) return;
    const dx = e.clientX - drag.sx;
    const dy = e.clientY - drag.sy;
    if (!dragged && Math.hypot(dx, dy) < 4) return;
    if (!dragged) {
      dragged = true;
      canvas.setPointerCapture(e.pointerId);
      canvas.classList.add("panning");
    }
    const el = contentOf(canvas);
    if (!el) return;
    view.x = drag.vx + dx;
    view.y = drag.vy + dy;
    clampView(el);
    applyView();
  });
  const endDrag = (e) => {
    if (!drag || drag.id !== e.pointerId) return;
    canvas.classList.remove("panning");
    drag = null;
  };
  canvas.addEventListener("pointerup", endDrag);
  canvas.addEventListener("pointercancel", endDrag);
  canvas.addEventListener("dblclick", resetView);
  canvas.addEventListener("dragstart", (e) => e.preventDefault());
}

// A drag that panned must not also count as a click with the Remove/Paint tools.
result.addEventListener("click", (e) => {
  if (dragged) {
    e.stopImmediatePropagation();
    dragged = false;
  }
}, true);

window.addEventListener("resize", () => {
  const c = canvases.find((k) => contentOf(k));
  if (c) clampView(contentOf(c));
  applyView();
});

// A new image starts unzoomed.
$("original").addEventListener("load", resetView);

startWorker();
