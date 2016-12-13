#include "XCPotential.h"
#include "XCFunctional.h"
#include "Potential.h"
#include "TelePrompter.h"

#include "constants.h"

extern MultiResolutionAnalysis<3> *MRA;

using namespace std;
using namespace Eigen;

void XCPotential::setup(double prec) {
    XCOperator::setup(prec);
    calcDensity();
    setupXCInput();
    setupXCOutput();
    evaluateXCFunctional();
    calcEnergy();
    calcPotential();
    clearXCInput();
    clearXCOutput();
    this->potential[0].setup(prec);
    this->potential[1].setup(prec);
    this->potential[2].setup(prec);
}

void XCPotential::clear() {
    this->energy = 0.0;
    this->density.clear();
    this->potential[0].clear();
    this->potential[1].clear();
    this->potential[2].clear();
    if (this->gradient != 0) {
        if (this->gradient[0] != 0) this->gradient[0]->clear();
        if (this->gradient[1] != 0) this->gradient[1]->clear();
        if (this->gradient[2] != 0) this->gradient[2]->clear();
    }
    this->gradient = deletePtrArray<Density>(3, &this->gradient);
    XCOperator::clear();
}

void XCPotential::calcPotential() {
    if (this->xcOutput == 0) MSG_ERROR("XC output not initialized");

    bool lda = this->functional->isLDA();
    bool gga = this->functional->isGGA();
    bool spin = this->functional->isSpinSeparated();

    Timer timer;
    if (lda) {
        if (not spin) {
            calcPotentialLDA(Paired);
        } else {
            calcPotentialLDA(Alpha);
            calcPotentialLDA(Beta);
        }
    } else if (gga) {
        if (not spin) {
            calcPotentialGGA(Paired);
        } else {
            calcPotentialGGA(Alpha);
            calcPotentialGGA(Beta);
        }
    } else {
        MSG_FATAL("Invalid functional type");
    }
    timer.stop();
    double t = timer.getWallTime();
    int n = 0;
    n += this->potential[0].getNNodes();
    n += this->potential[1].getNNodes();
    n += this->potential[2].getNNodes();
    TelePrompter::printTree(0, "XC potential", n, t);
}

void XCPotential::calcPotentialLDA(int spin) {
    if (spin == Paired) {
        if (this->xcOutput[1] == 0) MSG_ERROR("Invalid XC output");
        this->potential[0].setReal(this->xcOutput[1]);
        this->xcOutput[1] = 0;
    } else if (spin == Alpha) {
        if (this->xcOutput[1] == 0) MSG_ERROR("Invalid XC output");
        this->potential[1].setReal(this->xcOutput[1]);
        this->xcOutput[1] = 0;
    } else if (spin == Beta) {
        if (this->xcOutput[2] == 0) MSG_ERROR("Invalid XC output");
        this->potential[2].setReal(this->xcOutput[2]);
        this->xcOutput[2] = 0;
    } else {
        MSG_FATAL("Invalid spin");
    }
}

void XCPotential::calcPotentialGGA(int spin) {
    FunctionTreeVector<3> xc_funcs;
    FunctionTreeVector<3> dRho_a;
    FunctionTreeVector<3> dRho_b;

    Density &rho_x = *this->gradient[0];
    Density &rho_y = *this->gradient[0];
    Density &rho_z = *this->gradient[0];

    if (spin == Paired) {
        if (this->xcOutput[1] == 0) MSG_ERROR("Invalid XC output");
        if (this->xcOutput[2] == 0) MSG_ERROR("Invalid XC output");
        xc_funcs.push_back(this->xcOutput[1]);
        xc_funcs.push_back(this->xcOutput[2]);
        xc_funcs.push_back(0);
        dRho_a.push_back(&rho_x.total());
        dRho_a.push_back(&rho_y.total());
        dRho_a.push_back(&rho_z.total());
        dRho_b.push_back(0);
        dRho_b.push_back(0);
        dRho_b.push_back(0);

        FunctionTree<3> *V = calcPotentialGGA(xc_funcs, dRho_a, dRho_b);
        this->potential[0].setReal(V);

        xc_funcs.clear();
        dRho_a.clear();
        dRho_b.clear();
    }
    if (spin == Alpha) {
        if (this->xcOutput[1] == 0) MSG_ERROR("Invalid XC output");
        if (this->xcOutput[3] == 0) MSG_ERROR("Invalid XC output");
        if (this->xcOutput[4] == 0) MSG_ERROR("Invalid XC output");
        xc_funcs.push_back(this->xcOutput[1]);
        xc_funcs.push_back(this->xcOutput[3]);
        xc_funcs.push_back(this->xcOutput[4]);
        dRho_a.push_back(&rho_x.alpha());
        dRho_a.push_back(&rho_y.alpha());
        dRho_a.push_back(&rho_z.alpha());
        dRho_b.push_back(&rho_x.beta());
        dRho_b.push_back(&rho_y.beta());
        dRho_b.push_back(&rho_z.beta());

        FunctionTree<3> *V = calcPotentialGGA(xc_funcs, dRho_a, dRho_b);
        this->potential[1].setReal(V);

        xc_funcs.clear();
        dRho_a.clear();
        dRho_b.clear();
    }
    if (spin == Beta) {
        if (this->xcOutput[2] == 0) MSG_ERROR("Invalid XC output");
        if (this->xcOutput[4] == 0) MSG_ERROR("Invalid XC output");
        if (this->xcOutput[5] == 0) MSG_ERROR("Invalid XC output");
        xc_funcs.push_back(this->xcOutput[2]);
        xc_funcs.push_back(this->xcOutput[5]);
        xc_funcs.push_back(this->xcOutput[4]);
        dRho_a.push_back(&rho_x.beta());
        dRho_a.push_back(&rho_y.beta());
        dRho_a.push_back(&rho_z.beta());
        dRho_b.push_back(&rho_x.alpha());
        dRho_b.push_back(&rho_y.alpha());
        dRho_b.push_back(&rho_z.alpha());

        FunctionTree<3> *V = calcPotentialGGA(xc_funcs, dRho_a, dRho_b);
        this->potential[2].setReal(V);

        xc_funcs.clear();
        dRho_a.clear();
        dRho_b.clear();
    }
}

FunctionTree<3>* XCPotential::calcPotentialGGA(FunctionTreeVector<3> &xc_funcs,
                                               FunctionTreeVector<3> &dRho_a,
                                               FunctionTreeVector<3> &dRho_b) {
    if (xc_funcs[0] == 0) MSG_ERROR("Invalid XC output");

    FunctionTreeVector<3> funcs;
    funcs.push_back(1.0, xc_funcs[0]);

    FunctionTree<3> *tmp_1 = 0;
    if (xc_funcs[1] != 0) {
        tmp_1 = calcGradDotPotDensVec(*xc_funcs[1], dRho_a);
        funcs.push_back(-2.0, tmp_1);
    }

    FunctionTree<3> *tmp_2 = 0;
    if (xc_funcs[2] != 0) {
        tmp_2 = calcGradDotPotDensVec(*xc_funcs[2], dRho_b);
        funcs.push_back(-1.0, tmp_2);
    }

    GridGenerator<3> G(this->max_scale);
    MWAdder<3> add(-1.0, this->max_scale);

    FunctionTree<3> *V = new FunctionTree<3>(*MRA);
    G(*V, funcs);
    add(*V, funcs, 0);
    funcs.clear(false);

    if (tmp_1 != 0) delete tmp_1;
    if (tmp_2 != 0) delete tmp_2;

    return V;
}

