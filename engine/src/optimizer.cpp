#include "isvg/optimizer.h"

#include <cmath>
#include <cstdlib>

namespace isvg {

void Optimizer::setup(StateContainer& s) {
    sc_ = &s;
    n_ = s.size();
    if (n_ < 2) std::abort();
    state_ = s.statePointers();
    grad_ = s.gradientPointers();
    saved_.assign(n_, 0.0);
    dir_.assign(n_, 0.0);
}

void Optimizer::initializeRun(bool ownDir) {
    iter_ = 0;
    step_ = st_.maxEps;
    resets_ = 0;
    lineIters_ = 0;
    knockOut_ = 0;
    stepSum_ = 0;
    ownDir_ = ownDir;
    blindCounter_ = 0;
    blindPeriod_ = st_.numSeeingSteps + st_.numBlindSteps;
}

void Optimizer::storeState() {
    for (int i = 0; i < n_; ++i) saved_[i] = *state_[i];
}

void Optimizer::recallState() {
    for (int i = 0; i < n_; ++i) *state_[i] = saved_[i];
}

void Optimizer::move(double step) {
    for (int i = 0; i < n_; ++i) *state_[i] = dir_[i] * step + *state_[i];
    sc_->afterMove();
}

void Optimizer::resetBestYet() {
    f_ = sc_->findPotential();
    fPrev_ = f_;
}

void Optimizer::forceCGReset(bool keepBest) {
    beta_ = 0;
    sc_->afterCGReset(keepBest);
    for (int i = 0; i < n_; ++i) dir_[i] = -*grad_[i];
    ++resets_;
}

double Optimizer::doQuadraticLineSearch() {
    // Two equal steps along the direction, a parabola through the three potentials, a third step
    // to its minimum; keep the best of the four positions.
    const double maxStep = st_.maxStepSize;
    double step0 = step_;
    if (0.0 < maxStep && maxStep < maxAbsDir_ * step0) {
        step0 = maxStep / maxAbsDir_;
        step_ = step0;
    }
    int it = 0;
    stepSum_ = 0;
    for (;;) {
        if (st_.storeState != 0) storeState();
        double t = 0.0;
        const double f0 = f_;
        move(step_);
        const double f1 = sc_->findPotential();
        f_ = f1;
        move(step_);
        const double f2 = sc_->findPotential();
        double h = step_;
        f_ = f2;
        const double a = ((f2 - (f1 + f1)) + f0) / (h * (h + h));
        double f3;
        bool moved3 = false;
        if (a > 0.0) {
            t = -((f2 - f0) / (h + h) + h * (a + a)) / (a + a);
            if (0.0 < maxStep) {
                if (maxStep < maxAbsDir_ * t) t = maxStep / maxAbsDir_;
                if (maxStep < (-t) * maxAbsDir_) t = (-maxStep) / maxAbsDir_;
            }
            if (t != 0.0) {
                move(t);
                f3 = sc_->findPotential();
                h = step_;
                f_ = f3;
                moved3 = true;
            }
        }
        if (!moved3) f3 = 1.0 + f0;
        double fm = f1;
        if (f0 <= f1) fm = f0;
        int best = f1 < f0 ? 1 : 0;
        bool keepCurrent;
        if (fm <= f2) {
            keepCurrent = f3 < fm;
        } else {
            best = 2;
            keepCurrent = !(f2 <= f3);
        }
        double f, moved;
        if (keepCurrent) {
            f = f_;
            moved = h + h + t;
        } else {
            moved = static_cast<double>(best) * h;
            if (st_.storeState == 0) {
                move(-((h + h + t) - moved));
            } else {
                recallState();
                move(moved);
            }
            f = sc_->findPotential();
            h = step_;
            f_ = f;
        }
        ++it;
        const double grown = h * st_.stepGrowth;
        const double s = (st_.minEps <= grown) ? grown : st_.minEps;
        const double prevSum = stepSum_;
        stepSum_ = moved + prevSum;
        step_ = (s <= st_.maxEps) ? s : st_.maxEps;
        if (f0 == 0.0 || (f0 - f) / f0 <= st_.quadRelTol || st_.quadMaxIter <= it) {
            if (st_.useEpsFromStep == 0) {
                step_ = step0;
            } else {
                const double e = std::fabs(moved + prevSum) * st_.stepFraction;
                const double lo = (st_.minEps <= e) ? e : st_.minEps;
                step_ = (lo <= st_.maxEps) ? lo : st_.maxEps;
            }
            lineIters_ += it;
            return f - f0;
        }
    }
}

double Optimizer::doBlindLineSearch() {
    // After a few "seeing" (quadratic) searches per period, take blind steps of the running
    // average length; a blind step that increases the potential is undone and ends the run.
    const double f0 = f_;
    int c = blindCounter_ % blindPeriod_;
    blindCounter_ = c;
    if (st_.numSeeingSteps <= c) {
        storeState();
        double step = blindStep_;
        const double maxStep = st_.maxStepSize;
        if (0.0 < maxStep) {
            if (maxStep < maxAbsDir_ * step) step = maxStep / maxAbsDir_;
            if (maxStep < (-step) * maxAbsDir_) step = (-maxStep) / maxAbsDir_;
        }
        move(step * 0.8);
        const double f = sc_->findPotential();
        f_ = f;
        if (f0 < f) {
            recallState();
            f_ = f0;
            c = 0;
            blindCounter_ = 0;
        } else {
            c = blindCounter_;
        }
        if (st_.numSeeingSteps <= c) {
            blindCounter_ = c + 1;
            return f_ - f0;
        }
    }
    doQuadraticLineSearch();
    if (iter_ == 0) blindStep_ = stepSum_;
    else blindStep_ = (1.0 - st_.blindBlendLatest) * blindStep_ + st_.blindBlendLatest * stepSum_;
    blindCounter_ = blindCounter_ + 1;
    return f_ - f0;
}

double Optimizer::doLineSearch() {
    if (0.0 < st_.maxStepSize && ownDir_) {
        maxAbsDir_ = std::fabs(dir_[0]);
        for (int i = 1; i < n_; ++i) {
            const double v = std::fabs(dir_[i]);
            if (maxAbsDir_ < v) maxAbsDir_ = v;
        }
    }
    // Types 0 and >5: quadratic; 5: blind (types 1-4 are not used by any parameter file).
    if (st_.lineSearchType == 5) return doBlindLineSearch();
    if (st_.lineSearchType == 0 || st_.lineSearchType > 5) return doQuadraticLineSearch();
    std::abort();
}

void Optimizer::doBasicCG() {
    f_ = sc_->findPotential();
    sc_->findGradient();
    double gg = 0.0;
    for (int i = 0; i < n_; ++i) {
        const double g = *grad_[i];
        dir_[i] = -g;
        gg = gg + g * g;
    }
    iter_ = 0;
    fPrev_ = f_;
    double cur = f_;
    if (st_.maxIter <= 0) return;
    for (;;) {
        fPrev_ = cur;
        doLineSearch();
        sc_->endOfIteration();
        const double diff = fPrev_ - f_;
        const double rel = (fPrev_ > 0.0) ? diff / fPrev_ : 0.0;
        if (st_.relTol > rel || (st_.absTol > diff && !(diff < 0.0))) {
            ++knockOut_;
            if (knockOut_ > st_.knockOutCountDown && iter_ > st_.minIter) return;
        } else {
            knockOut_ = 0;
        }
        sc_->findGradient();
        double gg2 = 0.0;
        for (int i = 0; i < n_; ++i) gg2 = gg2 + *grad_[i] * *grad_[i];
        beta_ = (1e-30 > std::fabs(gg)) ? 0.0 : gg2 / gg;
        double acc = 0.0;
        dir_[0] = beta_ * dir_[0] - *grad_[0];
        acc = acc - dir_[0] * *grad_[0];
        for (int i = 1; i < n_; ++i) {
            dir_[i] = dir_[i] * beta_ - *grad_[i];
            acc = acc - dir_[i] * *grad_[i];
        }
        if (0.0 >= acc || iter_ % st_.iterForRestart == 0) forceCGReset(true);
        ++iter_;
        if (iter_ >= st_.maxIter) return;
        cur = (f_ < fPrev_) ? f_ : fPrev_;
        gg = gg2;
    }
}

void Optimizer::doCG() {
    if (n_ == 0) std::abort();
    initializeRun(true);
    doBasicCG();
    f_ = sc_->findPotential();
}

}  // namespace isvg
