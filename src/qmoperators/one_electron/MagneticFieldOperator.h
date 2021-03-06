#pragma once

#include "H_B_dip.h"
#include "RankZeroTensorOperator.h"

/** @class MagneticFieldOperator
 *
 * @brief External magnetic field operator
 *
 * An external magnetic field interacts with the molecular dipole
 * moment. The operator is simply implemented as a scalar product of
 * the dipole moment operator with the magnetic field vector.
 *
 */

namespace mrchem {

class MagneticFieldOperator final : public ExternalFieldOperator {
public:
    MagneticFieldOperator(const Eigen::Vector3d &f,
                          std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
                          const mrcpp::Coord<3> &o)
            : field(f)
            , dipole(D, o) {
        RankZeroTensorOperator &d_x = this->dipole[0];
        RankZeroTensorOperator &d_y = this->dipole[1];
        RankZeroTensorOperator &d_z = this->dipole[2];

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &HMF = (*this);
        HMF = f[0] * d_x + f[1] * d_y + f[2] * d_z;
        HMF.name() = "B . mu_B";
    }

private:
    Eigen::Vector3d field;
    H_B_dip dipole;
};

} // namespace mrchem
