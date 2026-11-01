// Optimizer: nonlinear conjugate gradient
// (Fletcher–Reeves) with the "blind" line search (type 5) used by the ContourSmoother.
#pragma once

#include <vector>

namespace isvg {

// OptimizerSettings (0xd8 bytes; copied to Optimizer+8 per phase). Offsets are settings-relative.
struct OptimizerSettings {
    int iterForRestart = 5;
    int minIter = 10;
    int maxIter = 30;
    int knockOutCountDown = 10;
    int storeState = 1;  // (quadratic search restores the saved state on reject)
    int useEpsFromStep = 1;
    double minEps = 1e-6;
    double maxEps = 0.01;
    double stepGrowth = 0.5;
    double stepFraction = 0.1;
    double fdEps = 9.999999747378752e-05;  // (finite differences; unused here)
    int numBlindSteps = 50;
    int numSeeingSteps = 2;
    double blindBlendLatest = 0.5;
    double relTol = 1e-5;
    double absTol = 1e-5;
    int lineSearchType = 5;
    double quadRelTol = 0.001;
    int quadMaxIter = 1;
    double maxStepSize = 0.3;
};

// StateContainer: the problem being optimised (ContourSmoother implements it).
class StateContainer {
public:
    virtual ~StateContainer() = default;
    virtual int size() = 0;
    virtual double findPotential() = 0;
    virtual std::vector<double*> statePointers() = 0;
    virtual std::vector<double*> gradientPointers() = 0;
    virtual void findGradient() = 0;
    virtual void endOfIteration() {}
    virtual void afterMove() {}
    virtual void afterCGReset(bool /*keepBest*/) {}
};

class Optimizer {
public:
    void setup(StateContainer& s);
    void setSettings(const OptimizerSettings& st) { st_ = st; }
    void doCG();
    void forceCGReset(bool keepBest);
    void resetBestYet();
    int iteration() const { return iter_; }

private:
    void initializeRun(bool ownDir);
    void doBasicCG();
    double doLineSearch();
    double doBlindLineSearch();
    double doQuadraticLineSearch();
    void move(double step);
    void storeState();
    void recallState();

    OptimizerSettings st_;
    StateContainer* sc_ = nullptr;
    int n_ = 0;
    int iter_ = 0;
    std::vector<double*> state_;
    std::vector<double> saved_;
    std::vector<double*> grad_;
    std::vector<double> dir_;
    double f_ = 0;
    double fPrev_ = 0;
    double stepSum_ = 0;
    int resets_ = 0;
    double beta_ = 0;
    int lineIters_ = 0;
    int knockOut_ = 0;
    double maxAbsDir_ = 0;
    double step_ = 0;
    int blindPeriod_ = 0;
    double blindStep_ = 0;
    int blindCounter_ = 0;
    bool ownDir_ = false;
};

}  // namespace isvg
