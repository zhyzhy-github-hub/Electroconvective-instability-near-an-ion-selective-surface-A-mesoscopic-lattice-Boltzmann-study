#pragma once
#include "lattice1.h"
#include "poisson1.h"

void WriteBinary_All(const size_t step, const size_t NX_, const size_t NY_,
                     const double *ff0, const double *ff1,
                     const double *gp0, const double *gp1,
                     const double *gn0, const double *gn1,
                     const double *fV0, const double *fV1,
                     const double *rho, const double *ux, const double *uy,
                     const double *phi, const double *Ex, const double *Ey,
                     const double *Cp, const double *Cn);
void ReStartBinary_All(const size_t &step, const size_t NX_, const size_t NY_,
                       const double *ff0, const double *ff1,
                       const double *gp0, const double *gp1,
                       const double *gn0, const double *gn1,
                       const double *fV0, const double *fV1,
                       const double *rho, const double *ux, const double *uy,
                       const double *phi, const double *Ex, const double *Ey,
                       const double *Cp, const double *Cn);