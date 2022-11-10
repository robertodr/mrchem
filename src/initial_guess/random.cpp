/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2022 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#include "random.h"

#include <random>

#include <Eigen/Dense>

#include <MRCPP/MWFunctions>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "parallel.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/one_electron/PositionOperator.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

namespace initial_guess {
namespace random {
bool setup(OrbitalVector &XorY, OrbitalVector &Phi, double prec, int seed) {
    if (XorY.size() == 0) return false;

    XorY = orbital::deep_copy(Phi);

    mrcpp::print::separator(0, '~');
    print_utils::text(0, "Calculation     ", "Compute initial orbitals");
    print_utils::text(0, "Method          ", "Random linear combination of multipoles");
    print_utils::text(0, "Precision       ", print_utils::dbl_to_str(prec, 5, true));

    auto rng_seed = (seed == 0) ? std::random_device()() : seed;
    print_utils::text(0, "PRNG seed       ", std::to_string(rng_seed));

    if (orbital::size_singly(Phi)) {
        print_utils::text(0, "Restricted    ", "False");
    } else {
        print_utils::text(0, "Restricted    ", "True");
    }
    mrcpp::print::separator(0, '~', 2);

    // Separate alpha/beta from paired orbitals
    auto XorY_a = orbital::disjoin(XorY, SPIN::Alpha);
    auto XorY_b = orbital::disjoin(XorY, SPIN::Beta);

    // random number generator
    auto prng = std::default_random_engine(rng_seed);
    // uniform distribution in [-1.0, 1.0) from which we will sample the coefficients
    auto dist = std::uniform_real_distribution<double>(0.0, 1.0);
    auto sampler = [&]() { return dist(prng); };

    auto success = true;
    {
        Eigen::Vector3d coeffs = Eigen::Vector3d::NullaryExpr(sampler);
        success &= apply_random(XorY, coeffs, prec);
    }

    {
        Eigen::Vector3d coeffs = Eigen::Vector3d::NullaryExpr(sampler);
        success &= apply_random(XorY_a, coeffs, prec);
    }

    {
        Eigen::Vector3d coeffs = Eigen::Vector3d::NullaryExpr(sampler);
        success &= apply_random(XorY_b, coeffs, prec);
    }

    // Collect orbitals into one vector
    for (auto &orb : XorY_a) XorY.push_back(orb);
    for (auto &orb : XorY_b) XorY.push_back(orb);

    // orthogonalize WRT unpertubed set
    orbital::orthogonalize(prec, XorY, Phi);

    return success;
}

bool apply_random(OrbitalVector &XorY, const Eigen::Vector3d &coeffs, double prec) {
    if (XorY.size() == 0) return true;

    bool success = true;

    // TODO would using a random gauge origin make sense here?
    PositionOperator Op;

    for (auto &orb : XorY) {
        Op.setup(prec);

        auto op_orb = Op(orb);

        orbital::linear_combination(orb, coeffs, op_orb, prec);

        Op.clear();
    }

    return success;
}
} // namespace random
} // namespace initial_guess
} // namespace mrchem
