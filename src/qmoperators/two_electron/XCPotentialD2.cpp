#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "XCPotential.h"
#include "XCPotentialD2.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"
#include "utils/print_utils.h"

using mrcpp::FunctionTree;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

XCPotentialD2::XCPotentialD2(std::unique_ptr<mrdft::MRDFT> &F,
                             std::shared_ptr<OrbitalVector> Phi,
                             std::shared_ptr<OrbitalVector> X,
                             std::shared_ptr<OrbitalVector> Y,
                             bool mpi_shared)
        : XCPotential(F, Phi, mpi_shared)
        , orbitals_x(X)
        , orbitals_y(Y) {
    densities.push_back(Density(false)); // rho_0 total
    densities.push_back(Density(false)); // rho_0 alpha
    densities.push_back(Density(false)); // rho_0 beta
    densities.push_back(Density(false)); // rho_1 total
    densities.push_back(Density(false)); // rho_1 alpha
    densities.push_back(Density(false)); // rho_1 beta
}

/** @brief Prepare the operator for application
 *
 * @param[in] prec Apply precision
 *
 * Sequence of steps required to compute the XC potentials:
 *
 * 1) Compute density
 * 2) Setup xcfun input functions (gradients etc.)
 * 3) Evaluate xcfun
 * 4) Compute XC energy by integrating energy density
 * 5) Compute XC potential(s) from xcfun output functions
 * 6) Remove excess grid nodes based on precision
 * 7) Add extra grid nodes based on precision
 * 8) Clear internal functions in XCFunctional (density grid is kept)
 *
 */
mrcpp::FunctionTreeVector<3> XCPotentialD2::setupDensities(double prec, mrcpp::FunctionTree<3> &grid) {
    mrcpp::FunctionTreeVector<3> dens_vec;
    if (not this->mrdft->functional().isSpin()) {
        { // Unperturbed total density
            Timer timer;
            Density &rho = getDensity(DensityType::Total, 0);
            if (not rho.hasReal()) {
                rho.alloc(NUMBER::Real);
                mrcpp::copy_grid(rho.real(), grid);
                density::compute(prec, rho, *orbitals, DensityType::Total);
            }
            print_utils::qmfunction(2, "XC rho_0", rho, timer);
            dens_vec.push_back(std::make_tuple(1.0, &rho.real()));
        }
        { // Perturbed total density
            Timer timer;
            Density &rho = getDensity(DensityType::Total, 1);
            if (not rho.hasReal()) {
                rho.alloc(NUMBER::Real);
                mrcpp::copy_grid(rho.real(), grid);
                density::compute(prec, rho, *orbitals, *orbitals_x, *orbitals_y, DensityType::Total);
            }
            print_utils::qmfunction(2, "XC rho_1", rho, timer);
            dens_vec.push_back(std::make_tuple(1.0, &rho.real()));
        }
    } else {
        { // Unperturbed alpha density
            Timer timer;
            Density &rho = getDensity(DensityType::Alpha, 0);
            if (not rho.hasReal()) {
                rho.alloc(NUMBER::Real);
                mrcpp::copy_grid(rho.real(), grid);
                density::compute(prec, rho, *orbitals, DensityType::Alpha);
            }
            print_utils::qmfunction(2, "XC rho_0 (alpha)", rho, timer);
            dens_vec.push_back(std::make_tuple(1.0, &rho.real()));
        }
        { // Unperturbed beta density
            Timer timer;
            Density &rho = getDensity(DensityType::Beta, 0);
            if (not rho.hasReal()) {
                rho.alloc(NUMBER::Real);
                mrcpp::copy_grid(rho.real(), grid);
                density::compute(prec, rho, *orbitals, DensityType::Beta);
            }
            print_utils::qmfunction(2, "XC rho_0 (beta)", rho, timer);
            dens_vec.push_back(std::make_tuple(1.0, &rho.real()));
        }
        { // Perturbed alpha density
            Timer timer;
            Density &rho = getDensity(DensityType::Alpha, 1);
            if (not rho.hasReal()) {
                rho.alloc(NUMBER::Real);
                mrcpp::copy_grid(rho.real(), grid);
                density::compute(prec, rho, *orbitals, *orbitals_x, *orbitals_y, DensityType::Alpha);
            }
            print_utils::qmfunction(2, "XC rho_1 (alpha)", rho, timer);
            dens_vec.push_back(std::make_tuple(1.0, &rho.real()));
        }
        { // Perturbed beta density
            Timer timer;
            Density &rho = getDensity(DensityType::Beta, 1);
            if (not rho.hasReal()) {
                rho.alloc(NUMBER::Real);
                mrcpp::copy_grid(rho.real(), grid);
                density::compute(prec, rho, *orbitals, *orbitals_x, *orbitals_y, DensityType::Beta);
            }
            print_utils::qmfunction(2, "XC rho_1 (beta)", rho, timer);
            dens_vec.push_back(std::make_tuple(1.0, &rho.real()));
        }
    }
    return dens_vec;
}

} // namespace mrchem
