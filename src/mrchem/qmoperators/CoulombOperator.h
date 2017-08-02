#ifndef COULOMBOPERATOR_H
#define COULOMBOPERATOR_H

#include "QMPotential.h"
#include "Density.h"

class CoulombOperator : public QMPotential {
public:
    CoulombOperator(PoissonOperator &P, OrbitalVector &phi)
        : poisson(&P),
          orbitals(&phi),
#ifdef HAVE_MPI
          density(false, true){//Use shared memory. 
#else
          density(false) {
#endif
    }
    virtual ~CoulombOperator() { }

protected:
    PoissonOperator *poisson;   // Pointer to external object
    OrbitalVector *orbitals;    // Pointer to external object
    Density density;            // Density that defines the potential
};

#endif // COULOMBOPERATOR_H
