#include "isvg/exporter.h"

#include <cstdio>

namespace isvg {

void Exporter::generateCompoundCurves() {
    const int nr = static_cast<int>(vi_.regions.size());
    const auto& refs = bf_.references();
    std::vector<std::vector<int>> children(nr);
    for (int u = 0; u < nr; ++u)
        if (vi_.regions[u].parent != -1) children[vi_.regions[u].parent].push_back(u);

    compound_.assign(nr, {});
    for (int p = 0; p < nr; ++p) {
        // Runs of consecutive child curves bordering p.
        std::vector<std::vector<std::pair<int, int>>> runs;
        for (int c : children[p]) {
            std::vector<std::pair<int, int>> run;
            const int cnt = static_cast<int>(refs[c].size());
            for (int i = 0; i < cnt; ++i) {
                if (refs[c][i].neighbor == p) {
                    run.emplace_back(c, i);
                } else if (!run.empty()) {
                    runs.push_back(run);
                    run.clear();
                }
            }
            if (!run.empty()) runs.push_back(run);
        }
        const int k = static_cast<int>(runs.size());
        if (k == 0) continue;
        // Link each run to the (last) run that starts where it ends, then follow the links.
        std::vector<int> next(k, -1);
        for (int i = 0; i < k; ++i) {
            const CurveReference& last = refs[runs[i].front().first][runs[i].back().second];
            for (int j = 0; j < k; ++j) {
                const CurveReference& first = refs[runs[j].front().first][runs[j].front().second];
                if (first.startNode == last.endNode) next[i] = j;
            }
        }
        for (int s = 0; s < k; ++s) {
            if (next[s] < 0) continue;
            std::vector<std::pair<int, int>> loop;
            int cur = s, nx = next[s];
            while (nx >= 0) {
                loop.insert(loop.end(), runs[cur].begin(), runs[cur].end());
                next[cur] = -1;
                cur = nx;
                nx = next[cur];
            }
            compound_[p].push_back(std::move(loop));
        }
    }
}

namespace {

void appendFixed2(std::string& out, double v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.2f", v);
    out += buf;
}

void appendPoint(std::string& out, const double p[2]) {
    appendFixed2(out, p[0]);
    out += ' ';
    appendFixed2(out, p[1]);
}

}  // namespace

void Exporter::writeLineOrCurve(std::string& out, int region, const CurveReference& ref, int startFlag, bool hole) {
    if (ref.type == 1 || ref.type == 2) {
        const BezierCurve& c = bf_.curves()[ref.index];
        const bool reversed = (ref.type == 1) == hole;
        double pts[4][2];
        for (int i = 0; i < 4; ++i) {
            const int s = reversed ? 3 - i : i;
            pts[i][0] = c.p[s][0];
            pts[i][1] = c.p[s][1];
        }
        if (startFlag == 1 || startFlag == 2) {
            out += " M ";
            appendPoint(out, pts[0]);
        }
        out += " C ";
        appendPoint(out, pts[1]);
        out += ' ';
        appendPoint(out, pts[2]);
        out += ' ';
        appendPoint(out, pts[3]);
        return;
    }
    if (ref.type != 0) return;
    const Region& reg = vi_.regions[region];
    const Node& a = vi_.nodes[reg.nodes[ref.index]];
    const Node& b = vi_.nodes[reg.nodes[(ref.index + 1) % static_cast<int>(reg.nodes.size())]];
    const double pa[2] = {a.x, a.y}, pb[2] = {b.x, b.y};
    const double* from = hole ? pb : pa;
    const double* to = hole ? pa : pb;
    if (startFlag == 1 || startFlag == 2) {
        out += " M ";
        appendPoint(out, from);
    }
    out += " L ";
    appendPoint(out, to);
}

std::vector<int> Exporter::exportedRegions() const {
    std::vector<int> ids;
    for (int u = 0; u < static_cast<int>(vi_.regions.size()); ++u)
        if (vi_.regions[u].color[3] > 10) ids.push_back(u);
    return ids;
}

std::string Exporter::writeSvg(const std::vector<uint8_t>* background, const std::vector<uint32_t>* colors) {
    // computeImAndPtSizes(72): image size from the node extent; points = pixels at 72 dpi.
    int imW = 0, imH = 0;
    for (const Node& nd : vi_.nodes) {
        if (static_cast<double>(imW) < nd.x) imW = static_cast<int>(nd.x + 0.5);
        if (static_cast<double>(imH) < nd.y) imH = static_cast<int>(nd.y + 0.5);
    }
    const int ptW = static_cast<int>((72.0f * static_cast<float>(imW)) / 72.0f);
    const int ptH = static_cast<int>((72.0f * static_cast<float>(imH)) / 72.0f);

    std::string out;
    out.reserve(1 << 20);
    out += "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
    out += "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";
    out += "<svg width=\"" + std::to_string(ptW) + "pt\" height=\"" + std::to_string(ptH) + "pt\" viewBox=\"0 0 " +
           std::to_string(imW) + " " + std::to_string(imH) + "\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">";

    const auto& refs = bf_.references();
    const int nr = static_cast<int>(vi_.regions.size());
    for (int u = 0; u < nr; ++u) {
        const Region& reg = vi_.regions[u];
        if (!(reg.color[3] > 10)) continue;  // SvgWriter::includeRegion
        // the export settings: shapes made background are skipped; the holes they cut
        // into their parent stay, so the area becomes transparent.
        if (background && u < static_cast<int>(background->size()) && (*background)[u]) continue;
        char col[16];
        if (colors && u < static_cast<int>(colors->size()) && ((*colors)[u] >> 24) != 0)
            std::snprintf(col, sizeof(col), "%06x", static_cast<unsigned>((*colors)[u] & 0xffffff));
        else
            std::snprintf(col, sizeof(col), "%02x%02x%02x", reg.color[2], reg.color[1], reg.color[0]);
        out += "\n<path fill=\"#";
        out += col;
        out += "\" opacity=\"";
        if (reg.color[3] < 0xfb) {
            char a[32];
            std::snprintf(a, sizeof(a), "%0.2f", static_cast<double>(reg.color[3]) / 255.0);
            out += a;
        } else {
            out += "1.00";
        }
        out += "\" d=\"";
        for (int i = 0; i < static_cast<int>(refs[u].size()); ++i) writeLineOrCurve(out, u, refs[u][i], i == 0 ? 1 : 0, false);
        // Holes: each loop is walked backwards (list prev pointers) with every piece reversed.
        for (const auto& loop : compound_[u])
            for (size_t i = loop.size(); i-- > 0;)
                writeLineOrCurve(out, loop[i].first, refs[loop[i].first][loop[i].second], i + 1 == loop.size() ? 2 : 0,
                                 true);
        out += " Z\" />";
    }
    out += "\n</svg>\n";
    // The original writes through a text-mode ofstream on Windows: every '\n' becomes "\r\n".
    std::string crlf;
    crlf.reserve(out.size() + out.size() / 64);
    for (char ch : out) {
        if (ch == '\n') crlf += '\r';
        crlf += ch;
    }
    return crlf;
}

}  // namespace isvg
