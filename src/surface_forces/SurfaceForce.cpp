#include "surface_forces/SurfaceForce.h"

#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"

#include "qmoperators/one_electron/NablaOperator.h"
#include <vector>
#include "qmoperators/one_electron/NuclearGradientOperator.h"

#include "chemistry/Molecule.h"
#include "chemistry/Nucleus.h"
#include "chemistry/PhysicalConstants.h"
#include <fstream>
#include <unsupported/Eigen/CXX11/Tensor>

#include "mrdft/Factory.h"
#include "mrdft/MRDFT.h"
#include "qmoperators/two_electron/XCOperator.h"
#include "mrdft/Functional.h"

#include "qmoperators/one_electron/HessianOperator.h"
#include "tensor/RankOneOperator.h"

#include "surface_forces/lebvedev.h"
#include <string>
#include <iostream>
#include <filesystem>

#include <MRCPP/Timer>

extern mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

using namespace Eigen;
using namespace mrchem;
using nlohmann::json;

namespace surface_force {

MatrixXd nuclearEfield(const MatrixXd &nucPos, const VectorXd &nucCharge, const VectorXd &nucSmoothing, const MatrixXd gridPos) {
    int nGrid = gridPos.rows();
    int nNuc = nucPos.rows();
    MatrixXd Efield = MatrixXd::Zero(nGrid, 3);
    MatrixXd r = MatrixXd::Zero(nGrid, 3);
    double temp;
    for (int i = 0; i < nNuc; ++i) {
        for (int j = 0; j < nGrid; j++)
        {
            r.row(j) = nucPos.row(i) - gridPos.row(j);
            temp = r.row(j).norm();
            temp = temp * temp * temp;
            Efield.row(j) += nucCharge(i) * r.row(j) / temp;
        }
    }
    return Efield;
}

mrcpp::ComplexFunction calcPotential(Density &rho, mrcpp::PoissonOperator &poisson, double prec) {
    mrcpp::ComplexFunction V(false);
    V.alloc(mrchem::NUMBER::Real);
    mrcpp::apply(prec, V.real(), poisson, rho.real());
    return V;
}

MatrixXd electronicEfield(mrchem::OrbitalVector &negEfield, const MatrixXd &gridPos) {
    int nGrid = gridPos.rows();
    MatrixXd Efield = MatrixXd::Zero(nGrid, 3);
    for (int i = 0; i < nGrid; i++)
    {
        std::array<double, 3> pos = {gridPos(i, 0), gridPos(i, 1), gridPos(i, 2)};
        Efield(i, 0) = -negEfield[0].real().evalf(pos);
        Efield(i, 1) = -negEfield[1].real().evalf(pos);
        Efield(i, 2) = -negEfield[2].real().evalf(pos);
    }
    return Efield;
}

std::vector<Eigen::Matrix3d> maxwellStress(const Molecule &mol, mrchem::OrbitalVector &negEfield, const MatrixXd &gridPos, double prec){
    int nGrid = gridPos.rows();
    int nNuc = mol.getNNuclei();

    Eigen::MatrixXd nucPos(nNuc, 3);
    Eigen::VectorXd nucCharge(nNuc);
    Eigen::VectorXd nucSmoothing(nNuc);
    for (int i = 0; i < nNuc; i++)
    {
        std::array<double, 3> coord = mol.getNuclei()[i].getCoord();
        nucPos(i, 0) = coord[0];
        nucPos(i, 1) = coord[1];
        nucPos(i, 2) = coord[2];
        nucCharge(i) = mol.getNuclei()[i].getCharge();
        double tmp = 0.00435 * prec / std::pow(nucCharge(i), 5.0);
        nucSmoothing(i) =  std::cbrt(tmp);
    }

    MatrixXd Efield = electronicEfield(negEfield, gridPos) + nuclearEfield(nucPos, nucCharge, gridPos, nucSmoothing);

    std::vector<Eigen::Matrix3d> stress(nGrid);
    for (int i = 0; i < nGrid; i++) {
        for (int i1 = 0; i1 < 3; i1++) {
            for (int i2 = 0; i2 < 3; i2++) {
                stress[i](i1, i2) = Efield(i, i1) * Efield(i, i2);
            }
        }
        for (int i1 = 0; i1 < 3; i1++){
            stress[i](i1, i1) = stress[i](i1, i1) - 0.5 * (Efield(i, 0) * Efield(i, 0) + Efield(i, 1) * Efield(i, 1) + Efield(i, 2) * Efield(i, 2));
        }
        for (int i1 = 0; i1 < 3; i1++){
            for (int i2 = 0; i2 < 3; i2++){
                stress[i](i1, i2) *= 1.0 / (4 * M_PI);
            }
        }
    }
    return stress;
}

std::vector<Matrix3d> xcStress(const Molecule &mol, const Density &rho, std::shared_ptr<XCOperator> XC_p, const MatrixXd &gridPos, double prec){
    int nGrid = gridPos.rows();

    std::vector<Matrix3d> stress(nGrid);

    Eigen::MatrixXd rhoGrid(nGrid, 1);
    std::array<double, 3> pos;
    for (int i = 0; i < nGrid; i++)
    {
        pos[0] = gridPos(i, 0);
        pos[1] = gridPos(i, 1);
        pos[2] = gridPos(i, 2);
        rhoGrid(i, 0) = rho.real().evalf(pos);
    }

    VectorXd xcGrid(nGrid);
    VectorXd vxcGrid(nGrid);

    for (int i = 0; i < nGrid; i++){
        pos[0] = gridPos(i, 0);
        pos[1] = gridPos(i, 1);
        pos[2] = gridPos(i, 2);
        xcGrid(i) = std::get<1>(XC_p->potential->xc_vec[0])->evalf(pos);
        vxcGrid(i) = std::get<1>(XC_p->potential->xc_vec[1])->evalf(pos);
        for (int i1 = 0; i1 < 3; i1++) {
            for (int i2 = 0; i2 < 3; i2++) {
                stress[i](i1, i2) = 0.0;
            }
        }
        for (int i1 = 0; i1 < 3; i1++)
        {
            stress[i](i1, i1) = xcGrid(i) - vxcGrid(i) * rhoGrid(i);
        }
        
    }

    return stress;
}

std::vector<Matrix3d> kineticStress(const Molecule &mol, OrbitalVector &Phi, std::vector<OrbitalVector> &nablaPhi, OrbitalVector &hessRho, double prec, const MatrixXd &gridPos){

    // original formula for kinetic stress:
    // sigma_ij = 0.5 \sum_k phi_k del_i del_j phi_k - (del_i phi_k) (del_j phi_k)
    // using the product rule for the first term:
    // sigma_ij = 0.5 * del_i del_j rho - 2 * \sum_k (del_i phi_k) (del_j phi_k)
    // That way, only second derivatives of density is needed which speeds things up. quite a lot.

    int nGrid = gridPos.rows();
    int nOrbs = Phi.size();

    double orbVal;
    std::vector<Matrix3d> stress(nGrid);
    for (int i = 0; i < nGrid; i++) {
        stress[i] = Matrix3d::Zero();
    }
    std::array<double, 3> pos;
    double n1, n2, n3;
    int occ;
    for (int iOrb = 0; iOrb < Phi.size(); iOrb++) {
        occ = Phi[iOrb].occ();
    
        for (int i = 0; i < nGrid; i++) {
            pos[0] = gridPos(i, 0);
            pos[1] = gridPos(i, 1);
            pos[2] = gridPos(i, 2);
            n1 = nablaPhi[iOrb][0].real().evalf(pos);
            n2 = nablaPhi[iOrb][1].real().evalf(pos);
            n3 = nablaPhi[iOrb][2].real().evalf(pos);
            stress[i](0, 0) -= occ * n1 * n1;
            stress[i](1, 1) -= occ * n2 * n2;
            stress[i](2, 2) -= occ * n3 * n3;
            stress[i](0, 1) -= occ * n1 * n2;
            stress[i](0, 2) -= occ * n1 * n3;
            stress[i](1, 2) -= occ * n2 * n3;
        }
    }
    // loop over grid
    for (int i = 0; i < nGrid; i++) {
        pos[0] = gridPos(i, 0);
        pos[1] = gridPos(i, 1);
        pos[2] = gridPos(i, 2);
        stress[i](0, 0) += 0.25 * hessRho[0].real().evalf(pos);
        stress[i](1, 1) += 0.25 * hessRho[1].real().evalf(pos);
        stress[i](2, 2) += 0.25 * hessRho[2].real().evalf(pos);
        stress[i](0, 1) += 0.25 * hessRho[3].real().evalf(pos);
        stress[i](0, 2) += 0.25 * hessRho[4].real().evalf(pos);
        stress[i](1, 2) += 0.25 * hessRho[5].real().evalf(pos);
        // symmetrize stress tensor
        stress[i](1, 0) = stress[i](0, 1);
        stress[i](2, 0) = stress[i](0, 2);
        stress[i](2, 1) = stress[i](1, 2);
    }
    return stress;
}

/**
 * Calculates the distance to the nearest neighbor for each point in the given position matrix.
 *
 * @param pos The matrix containing the positions of the points. Shape (nPoints, 3).
 * @return A vector containing the distances to the nearest neighbor for each point.
 */
VectorXd distanceToNearestNeighbour(MatrixXd pos){
    int n = pos.rows();
    VectorXd dist(n);
    double temp;
    for (int i = 0; i < n; i++){
        dist(i) = (pos.row(i) - pos.row((i + 1) % n)).norm();
        for (int j = 0; j < n; j++){
            if (i != j){
                temp = (pos.row(i) - pos.row(j)).norm();
                if (temp < dist(i)){
                    dist(i) = temp;
                }
            }
        }
    }
    return dist;
}

/**
 * Calculates the forces using surface integrals for a given molecule and orbital vector.
 *
 * @param mol The molecule for which to calculate the forces.
 * @param Phi The orbital vector obtained from the SCF calculation.
 * @param prec The precision value used in the calculation.
 * @param json_fock The JSON object containing the Fock matrix settings.
 * @return The matrix of forces, shape (nAtoms, 3).
 */
Eigen::MatrixXd surface_forces(mrchem::Molecule &mol, mrchem::OrbitalVector &Phi, double prec, const json &json_fock) {

    // setup density
    mrchem::Density rho(false);
    mrchem::density::compute(prec, rho, Phi, DensityType::Total);

    // setup operators and potentials
    auto poisson_op = mrcpp::PoissonOperator(*mrchem::MRA, prec);
    int derivOrder = 1;
    auto mrcd = std::make_shared<mrcpp::BSOperator<3>>(*mrchem::MRA, derivOrder);
    mrchem::NablaOperator nabla(mrcd);
    nabla.setup(prec);
    double abs_prec = prec / mrchem::orbital::get_electron_number(Phi);
    mrcpp::ComplexFunction pot = calcPotential(rho, poisson_op, abs_prec);
    mrchem::OrbitalVector negEfield = nabla(pot);

    // set up operators for kinetic stress:
    int derivOrder1 = 1;
    auto D1 = std::make_shared<mrcpp::BSOperator<3>>(*mrchem::MRA, derivOrder1);
    int derivOrder2 = 2;
    auto D2 = std::make_shared<mrcpp::BSOperator<3>>(*mrchem::MRA, derivOrder2);
    mrchem::HessianOperator hess(D1, D2, prec);
    hess.setup(prec);

    std::vector<mrchem::OrbitalVector> nablaPhi;
    mrchem::OrbitalVector hessRho = hess(rho);
    for (int i = 0; i < Phi.size(); i++) {
        nablaPhi.push_back(nabla(Phi[i]));
    }

    // setup xc stuff:
    int order = 0;
    bool shared_memory = json_fock["xc_operator"]["shared_memory"];
    auto json_xcfunc = json_fock["xc_operator"]["xc_functional"];
    auto xc_spin = json_xcfunc["spin"];
    auto xc_cutoff = json_xcfunc["cutoff"];
    auto xc_funcs = json_xcfunc["functionals"];
    auto xc_order = order + 1;
    auto funcVectorShared = std::make_shared<mrcpp::MPI_FuncVector>(Phi);

    mrdft::Factory xc_factory(*MRA);
    xc_factory.setSpin(xc_spin);
    xc_factory.setOrder(xc_order);
    xc_factory.setDensityCutoff(xc_cutoff);
    for (const auto &f : xc_funcs) {
        auto name = f["name"];
        auto coef = f["coef"];
        xc_factory.setFunctional(name, coef);
    }
    std::unique_ptr<mrdft::MRDFT> mrdft_p = xc_factory.build();
    std::shared_ptr<XCOperator> XC_p = std::make_shared<XCOperator>(mrdft_p, funcVectorShared, shared_memory);
    XC_p->potential->setup(prec);

    
    int numAtoms = mol.getNNuclei();
    int numOrbitals = Phi.size();

    Eigen::MatrixXd forces = Eigen::MatrixXd::Zero(numAtoms, 3);
    Vector3d center;
    std::array<double, 3> coord;

    Eigen::MatrixXd posmatrix = Eigen::MatrixXd::Zero(numAtoms, 3);
    for (int i = 0; i < numAtoms; i++) {
        coord = mol.getNuclei()[i].getCoord();
        posmatrix(i, 0) = coord[0];
        posmatrix(i, 1) = coord[1];
        posmatrix(i, 2) = coord[2];
    }
    VectorXd dist = distanceToNearestNeighbour(posmatrix);

    std::filesystem::path p = __FILE__;
    std::filesystem::path parent_dir = p.parent_path();
    std::string filename = parent_dir.string() + "/lebvedev.txt";
    std::string tinyPoints = parent_dir.string() + "/lebvedev_tiny.txt";

    double tinyRadius = 0.15;
    double radius = 0.6;
    int nRad = 11;
    VectorXd radii(nRad);
    double step;
    for (int i = 0; i < nRad; i++) {
        step = (1.0 * i) / (1.0 * nRad);
        radii(i) = -0.2 + 0.4 * step;
    }

    for (int iAtom = 0; iAtom < numAtoms; iAtom++) {
        radius = dist(iAtom) *.5;
        coord = mol.getNuclei()[iAtom].getCoord();
        center << coord[0], coord[1], coord[2];
        LebedevIntegrator tintegrator(tinyPoints, tinyRadius, center);
        MatrixXd tinyPos = tintegrator.getPoints();
        VectorXd tinyWeights = tintegrator.getWeights();
        for (int iTiny = 0; iTiny < tintegrator.n; iTiny++){
            // Vector3d x = center + tinyPos.row(iTiny).transpose();
            // std::cerr << "x " << x(0) << " " << x(1) << " " << x(2) << "\n";
            // std::cerr << "center " << center << "\n";
            // std::cerr << "tinyPos " << tinyPos.row(iTiny) << "\n";
            LebedevIntegrator integrator(filename, radius, tinyPos.row(iTiny).transpose());
            MatrixXd gridPos = integrator.getPoints();
            VectorXd weights = integrator.getWeights();
            MatrixXd normals = integrator.getNormals();
            std::vector<Matrix3d> xStress = xcStress(mol, rho, XC_p, gridPos, prec);
            std::vector<Matrix3d> kstress = kineticStress(mol, Phi, nablaPhi, hessRho, prec, gridPos);
            std::vector<Matrix3d> mstress = maxwellStress(mol, negEfield, gridPos, prec);
            std::vector<Matrix3d> stress(integrator.n);
            for (int i = 0; i < integrator.n; i++){
                stress[i] = xStress[i] + kstress[i] + mstress[i];
                forces.row(iAtom) -= stress[i] * normals.row(i).transpose() * weights(i) * tinyWeights(iTiny) / (4.0 * M_PI * tinyRadius * tinyRadius);
            }
        }
    }

    for (int iAtom = 0; iAtom < numAtoms; iAtom++) {
        std::cerr << "forces " << forces(iAtom, 0) << " " << forces(iAtom, 1) << " " << forces(iAtom, 2) << std::endl;
    }

    hess.clear();
    nabla.clear();
    XC_p->potential->clear();
    return forces;
}

} // namespace surface_force
