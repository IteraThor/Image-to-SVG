#include "isvg/engine_params.h"

#include <cstdlib>
#include <cstring>
#include <map>
#include <sstream>
#include <vector>

namespace isvg {

const char* engineParamText(const std::string& name) {
    for (const EngineParamFile* f = kEngineParamFiles; f->name; ++f)
        if (name == f->name) return f->text;
    return nullptr;
}

std::string engineParamFileFor(int imageType, int quality) {
    static const char* const kFiles[9] = {"EP_PHOTO_HIG",   "EP_PHOTO_MED",   "EP_PHOTO_LOW",
                                          "EP_LOGO_HIG_AA", "EP_LOGO_MED_AA", "EP_LOGO_LOW_AA",
                                          "EP_LOGO_HIG",    "EP_LOGO_MED",    "EP_LOGO_LOW"};
    const int base = imageType == 1 ? 3 : imageType == 2 ? 6 : 0;
    return kFiles[base + quality];
}

namespace {

using Values = std::map<std::string, std::vector<std::string>>;

Values tokenize(const char* text) {
    Values v;
    std::istringstream in(text);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const size_t sp = line.find(' ');
        if (sp == std::string::npos || line.find("::") == std::string::npos) continue;
        const std::string name = line.substr(0, sp);
        std::istringstream rest(line.substr(sp));
        std::string tok;
        std::vector<std::string> vals;
        while (rest >> tok)
            if (tok.find_first_not_of('.') != std::string::npos) vals.push_back(tok);
        v[name] = vals;
    }
    return v;
}

struct Reader {
    const Values& v;
    const std::vector<std::string>* get(const std::string& name) const {
        auto it = v.find(name);
        return it == v.end() || it->second.empty() ? nullptr : &it->second;
    }
    void i(const std::string& name, int& out, int idx = 0) const {
        if (auto* s = get(name)) out = std::atoi((*s)[idx].c_str());
    }
    void b(const std::string& name, bool& out) const {
        int t = out;
        i(name, t);
        out = t != 0;
    }
    void f(const std::string& name, float& out, int idx = 0) const {
        if (auto* s = get(name)) out = std::strtof((*s)[idx].c_str(), nullptr);
    }
    void d(const std::string& name, double& out, int idx = 0) const {
        if (auto* s = get(name)) out = std::strtod((*s)[idx].c_str(), nullptr);
    }
};

}  // namespace

EngineParams parseEngineParams(const char* text) {
    const Values v = tokenize(text);
    const Reader r{v};
    EngineParams p;
    r.b("Shared::is_anti_aliased", p.antiAliased);
    r.i("Shared::image_type_code", p.imageTypeCode);
    r.f("PaletteFinder::cost_per_color", p.paletteFinder.costPerColor);
    r.b("PaletteFinder::use_16_hist_bins", p.paletteFinder.use16HistBins);
    r.d("EdgeImage::edge_weight_scale", p.edge.scale);
    r.d("EdgeImage::edge_weight_floor", p.edge.floor);

    PixelSegmenterParams& s = p.segmenter;
    r.i("Segmenter::max_iterations", s.max_iterations);
    r.f("Segmenter::lambda_pre", s.lambda_pre);
    r.f("Segmenter::lambda_initial", s.lambda_initial);
    r.f("Segmenter::lambda_final", s.lambda_final);
    r.i("Segmenter::min_num_pixels", s.min_num_pixels);
    r.i("Segmenter::num_iter_btwn_renum", s.num_iter_btwn_renum);
    r.f("Segmenter::rise_fraction", s.rise_fraction);
    r.f("Segmenter::size_penalty_weight", s.size_penalty_weight);
    r.f("Segmenter::size_penalty_power", s.size_penalty_power);
    r.b("Segmenter::use_over_seg_diag_extr", p.useOverSegDiagExtr);

    r.f("SubPixelSegmenter::lambda_f", p.subSegmenter.lambda_f);
    r.f("SubPixelSegmenter::beach_lambda_f", p.subSegmenter.beach_lambda_f);
    r.f("SubPixelSegmenter::color_cluster_margin", p.subSegmenter.color_cluster_margin);
    r.i("SubPixelSegmenter::use_connected_24_iteration", p.subSegmenter.use_connected_24_iteration);
    r.f("SubPixelSegmentation::self_weight_eps", p.subSegmentation.self_weight_eps);
    r.f("SubPixelSegmentation::self_weight_penalty_low", p.subSegmentation.self_weight_penalty_low);
    r.f("SubPixelSegmentation::self_weight_penalty_high", p.subSegmentation.self_weight_penalty_high);

    ContourSmootherParams& c = p.smoother;
    for (int k = 0; k < 3; ++k) {
        r.i("ContourSmoother::measurement_types", c.measurementTypes[k], k);
        r.i("ContourSmoother::optimizer_type", c.optimizerType[k], k);
        r.i("ContourSmoother::prior_types", c.priorTypes[k], k);
        r.d("ContourSmoother::prior_strengths", c.priorStrengths[k], k);
        r.d("ContourSmoother::length_penalty_weights", c.lengthPenaltyWeights[k], k);
        r.i("ContourSmoother::do_puncture_corners", c.doPunctureCorners[k], k);
        r.i("ContourSmoother::is_enabled", c.isEnabled[k], k);
    }
    r.d("ContourSmoother::length_barrier_weight", c.lengthBarrierWeight);
    r.d("ContourSmoother::air_pressure_weight", c.airPressureWeight);
    r.d("ContourSmoother::perturbation_range", c.perturbationRange);
    r.d("ContourSmoother::puncture_flatness_thresh", c.punctureFlatnessThresh);
    r.d("ContourSmoother::puncture_corner_thresh", c.punctureCornerThresh);
    r.d("ContourSmoother::anti_inv_pot_meas_scale", c.antiInvPotMeasScale);
    r.d("ContourSmoother::anti_inv_pot_prior_scale", c.antiInvPotPriorScale);
    for (int k = 0; k < 3; ++k) {
        const std::string ph = "ContourSmoother::phase_" + std::to_string(k) + ".cg_";
        OptimizerSettings& o = c.phase[k];
        r.i(ph + "line_search_type", o.lineSearchType);
        r.i(ph + "iter_for_restart", o.iterForRestart);
        r.i(ph + "min_iter", o.minIter);
        r.i(ph + "max_iter", o.maxIter);
        r.i(ph + "knock_out_count_down", o.knockOutCountDown);
        r.i(ph + "use_eps_from_step", o.useEpsFromStep);
        r.i(ph + "quad_line_srch_max_iter", o.quadMaxIter);
        r.d(ph + "min_eps", o.minEps);
        r.d(ph + "max_eps", o.maxEps);
        r.d(ph + "step_fraction", o.stepFraction);
        r.d(ph + "rel_tol", o.relTol);
        r.d(ph + "abs_tol", o.absTol);
        r.d(ph + "quad_line_srch_rel_tol", o.quadRelTol);
        r.d(ph + "max_step_size", o.maxStepSize);
        r.i(ph + "num_blind_steps", o.numBlindSteps);
        r.i(ph + "num_seeing_steps", o.numSeeingSteps);
        r.d(ph + "blind_blend_latest", o.blindBlendLatest);
    }
    r.i("BezierFitter::max_iterations", p.bezier.maxIterations);
    r.d("BezierFitter::stat_thresh_initial", p.bezier.statThreshInitial);
    r.d("BezierFitter::stat_thresh_final", p.bezier.statThreshFinal);
    return p;
}

}  // namespace isvg
