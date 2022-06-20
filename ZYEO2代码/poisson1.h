#ifndef _POISSON1_H
#define _POISSON1_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>

#include "lattice1.h"

const double PI {2 * asin(1.0)};

const double dt_phi = 1.0;
const double dx_phi = 1.0;
const double c_phi_ = dx_phi / dt_phi;


const int C_xQ5[Q5] {0, 1, 0, -1,  0};
const int C_yQ5[Q5] {0, 0, 1,  0, -1};
const int OppQ5[Q5] = {0, 3, 4, 1, 2};

const double w0 = 1.0 / 3.0;
const double wi = 1.0 / 6.0;
// The Cs2 of poisson equation is defined by omega_i, Cs2_phi = alpha
const double alpha = 1.0 / 3.0;

const double w0_h = 0.0 / 3.0;
const double wi_h = 1.0 / 4.0;
const double w[Q5]   {w0  , wi  , wi  , wi  , wi  };
const double w_h[Q5] {w0_h, wi_h, wi_h, wi_h, wi_h};
const double tau_phi = 1.36;
const double omega_phi = 1.0 / tau_phi;

// const double D = (0.5 - tau_phi) * dt / 3;
const double D_phi = (0.5 - tau_phi) * dt_phi * alpha * c_phi_ * c_phi_;
extern double k_zeta;

const double C2q_to_eps = C2q / eps;

const size_t MAXSTEP_Poi {100000};
const double rela_max_Poi_err = 1.0e-6;

void Init_popu_poi(double *f0, double *f1, double *rho);

void poisson(double *f0, double *f1, double *phi,
             const double *Cp, const double *Cn,
             const int errflag, double &errsum);
void Poisson_solver(size_t &tot, double *f0, double *f1, double *phi, 
             const double *Cp, const double *Cn, 
             double *Ex, double *Ey, 
             double &phiall, double &rela_err);



inline void feqD2Q5_Poi(double* feq, double u) {
    feq[0] = (w0 - 1.0) * u;
    
    feq[1] = wi * u;
    feq[2] = feq[1]; 
    feq[3] = feq[1]; 
    feq[4] = feq[1]; 
    // feq[2] = wi * u;
    // feq[3] = wi * u;
    // feq[4] = wi * u;
}
inline void Anti_BB_Poi(double* f1, double* f0){
    // I don't know what it is for!!!
    f1[1] = -f0[3];
    f1[2] = -f0[4];
    f1[3] = -f0[1];
    f1[4] = -f0[2];
    // for(size_t k = 0; k < Q5; ++k) {
    //     f1[k] = -f1[OppQ5[k]];
    // }
}

inline void NEE_Q5(double* f0, double* f1, double* feq0, double* feq1){
    for(size_t k = 0; k < Q5; ++k){
        f0[k] = feq0[k] + f1[k] - feq1[k];
    }
}

inline void ABB_Poi(double* f1, double TD){
    // It is wrong formula!!!
    for(size_t k = 0; k < Q5; ++k) {
        f1[k] = -f1[OppQ5[k]] + 1.0 * w[k] * TD;
    }
}

inline void TopBC_Inamuro_D_Poi(double* f0, const double TD) {
    // double rho_temp = TD;
    // f0[4] = rho_temp * (1 - w0) - f0[3] - f0[1] - f0[2];
    f0[4] = TD * (1 - w0) - f0[3] - f0[1] - f0[2];
}
inline void BottomBC_Inamuro_D_Poi(double* f0, const double TD) {
    // double rho_temp = TD;
    // f0[2] = rho_temp * (1 - w0) - f0[3] - f0[1] - f0[4];
    f0[2] = TD * (1 - w0) - f0[3] - f0[1] - f0[4];
}
inline void RightBC_Inamuro_D_Poi(double* f0, const double TD) {
    // double rho_temp = TD;
    // f0[3] = rho_temp * (1 - w0) - f0[4] - f0[1] - f0[2];
    f0[3] = TD * (1 - w0) - f0[4] - f0[1] - f0[2];
}
inline void LeftBC_Inamuro_D_Poi(double* f0, const double TD) {
    f0[1] = TD * (1 - w0) - f0[3] - f0[4] - f0[2];
}

inline void BottomLeftBC_Inamuro_D_Poi(double* f0, double TD) {
    double rho_temp = TD;
    double temp = rho_temp * (1 - w0) - f0[4] - f0[3];
    f0[1] = temp * 0.5;
    f0[2] = temp * 0.5;
}
inline void BottomRightBC_Inamuro_D_Poi(double* f0, double TD) {
    double rho_temp = TD;
    double temp = rho_temp * (1 - w0) - f0[4] - f0[1];
    f0[3] = temp * 0.5;
    f0[2] = temp * 0.5;
}
inline void TopLeftBC_Inamuro_D_Poi(double* f0, double TD) {
    double rho_temp = TD;
    double temp = rho_temp * (1 - w0) - f0[2] - f0[3];
    f0[1] = temp * 0.5;
    f0[4] = temp * 0.5;
}
inline void TopRightBC_Inamuro_D_Poi(double* f0, double TD) {
    double rho_temp = TD;
    double temp = rho_temp * (1 - w0) - f0[2] - f0[1];
    f0[3] = temp * 0.5;
    f0[4] = temp * 0.5;
}

#endif