#ifndef MAGNETIZABILITY_H
#define MAGNETIZABILITY_H

#include <Eigen/Core>

class AngularMomentumOperator;
class DMOperator;

class Magnetizability {
public:
    Magnetizability() { }
    virtual ~Magnetizability() { }

    Eigen::MatrixXd get() const { return this->diamagnetic + this->paramagnetic; }
    Eigen::MatrixXd getDiamagnetic() const { return this->diamagnetic; }
    Eigen::MatrixXd getParamagnetic() const { return this->paramagnetic; }

    void compute(DMOperator &h_ij, OrbitalVector &orbs) { NOT_IMPLEMENTED_ABORT; }
    void compute(int j, AngularMomentumOperator &h_i, OrbitalVector &orbs_0, OrbitalVector &orbs_1) { NOT_IMPLEMENTED_ABORT; }

    friend std::ostream& operator<<(std::ostream &o, const Magnetizability &mag) {
        double w_au = 0.0;  // Only static magnetizability
        double isoPMau = mag.paramagnetic.trace()/3.0;
        double isoDMau = mag.diamagnetic.trace()/3.0;
        double isoTMau = isoDMau + isoPMau;

        double isoDMsi = isoDMau*78.9451185;
        double isoPMsi = isoPMau*78.9451185;
        double isoTMsi = isoTMau*78.9451185;

        int oldPrec = TelePrompter::setPrecision(10);
        o<<"                                                            "<<std::endl;
        o<<"================== Magnetizability tensor =================="<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Frequency:        (au)       " << std::setw(30) << w_au     <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"-------------------- Isotropic averages --------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Total            (au)        " << std::setw(30) << isoTMau  <<std::endl;
        o<<" Diamagnetic      (au)        " << std::setw(30) << isoDMau  <<std::endl;
        o<<" Paramagnetic     (au)        " << std::setw(30) << isoPMau  <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Total            (SI)        " << std::setw(30) << isoTMsi  <<std::endl;
        o<<" Diamagnetic      (SI)        " << std::setw(30) << isoDMsi  <<std::endl;
        o<<" Paramagnetic     (SI)        " << std::setw(30) << isoPMsi  <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"-------------------------- Total ---------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<< mag.diamagnetic + mag.paramagnetic                            <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"----------------------- Diamagnetic ------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<< mag.diamagnetic                                               <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"----------------------- Paramagnetic -----------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<< mag.paramagnetic                                              <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                                                            "<<std::endl;
        TelePrompter::setPrecision(oldPrec);
        return o;
    }
protected:
    Eigen::MatrixXd diamagnetic;
    Eigen::MatrixXd paramagnetic;
};

#endif // MAGNETIZABILITY_H