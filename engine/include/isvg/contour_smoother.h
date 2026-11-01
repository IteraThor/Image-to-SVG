// ContourSmoother: moves the traced pixel-corner nodes to
// sub-pixel positions by minimising (contour priors + rendering error of the GenerativeModel) with
// conjugate gradients, in three phases.
#pragma once

#include <vector>

#include "isvg/generative_model.h"
#include "isvg/optimizer.h"
#include "isvg/palette.h"
#include "isvg/vector_image.h"
#include "isvg/weights_finder.h"

namespace isvg {

struct ContourSmootherParams {  // EP_LOGO_HIG_AA defaults
    int measurementTypes[3] = {1, 1, 1};
    int priorTypes[3] = {3, 3, 3};  // (0 quadratic, 1 sqrt, 3 angular)
    double priorStrengths[3] = {0.1, 0.02, 0.9};
    double lengthPenaltyWeights[3] = {0.5, 0.1, 0.01};
    int doPunctureCorners[3] = {0, 1, 0};
    int isEnabled[3] = {1, 1, 1};
    int optimizerType[3] = {0, 0, 0};
    double airPressureWeight = 5.0;
    double perturbationRange = 0.1;
    double lengthBarrierWeight = 0.1;
    double punctureFlatnessThresh = 0.4;
    double punctureCornerThresh = 0.5;
    double antiInvPotMeasScale = 1.5;
    double antiInvPotPriorScale = 1.5;
    int inversionCheckPeriod = 10;
    OptimizerSettings phase[3];  // + 0xd8·p
    ContourSmootherParams() {
        phase[1].iterForRestart = 25;
        phase[1].minIter = 50;
        phase[1].maxIter = 225;
        phase[2].iterForRestart = 25;
        phase[2].minIter = 20;
        phase[2].maxIter = 300;
        phase[2].relTol = 2e-5;
        phase[2].absTol = 0.01;
    }
};

class ContourSmoother : public StateContainer {
public:
    ContourSmoother(VectorImage& vi, PixelMap& map, const Image& img, const Palette& palette, bool antiAliased,
                    const ContourSmootherParams& p);

    void execute(bool initPositions);
    void initialize(bool initPositions);
    void executePhase(int phase);  // executeCurrentPhase
    void finish();  // closing flagInversionNodes of execute
    // Puts runs of nodes that are straight up to noise exactly on their line (and the node where two
    // such runs meet on the lines' intersection), so straight edges do not come out wavy.
    void straightenEdges();
    GenerativeModel& model() { return gm_; }

    // StateContainer
    int size() override { return 2 * static_cast<int>(vi_.nodes.size()); }
    double findPotential() override;
    std::vector<double*> statePointers() override;
    std::vector<double*> gradientPointers() override;
    void findGradient() override;
    void endOfIteration() override;
    void afterMove() override;
    void afterCGReset(bool keepBest) override;

private:
    void initializeNodePositions();
    void initializeNodePosition(int i);
    void punctureCorners();
    // Smooths out short wiggles (pixel staircases) on the curved stretches, keeping nodes on straight
    // lines, corners, junctions and border nodes where they are; called by straightenEdges.
    void fairCurves(const std::vector<char>& onLine);
    double prior(int a, int b, int c);  // pot. member pointer
    void priorGradient(int a, int b, int c, double o1[2], double o2[2]);
    double quadraticPrior(int a, int b, int c) const;
    void quadraticPriorGradient(int a, int b, int c, double o1[2], double o2[2]) const;
    double angularPrior(int a, int b, int c) const;
    void angularPriorGradient(int a, int b, int c, double o1[2], double o2[2]) const;

    VectorImage& vi_;
    PixelMap& map_;
    const Image& img_;
    const Palette& palette_;
    bool antiAliased_;
    ContourSmootherParams p_;
    Optimizer opt_;
    GenerativeModel gm_;
    WeightsFinder wf_;
    std::vector<double> grad_;
    int phase_ = 0;
    int stepsThis_ = 1, stepsNext_ = 1;
    const double M_[4] = {-0.7071067811865475, -0.7071067811865475, -0.7071067811865475, 0.7071067811865475};
};

}  // namespace isvg
