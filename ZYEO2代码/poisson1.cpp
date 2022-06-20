// #pragma once
#include "poisson1.h"

void Init_scalar_poi(double* phi){
    for (size_t j = 0; j < NY; ++j) {
        for (size_t i = 0; i < NX; ++i) {
            const size_t index = scalar_index(i, j);
            phi[index] = 0.0;
            //if(j == 0){
            //    phi[index] = zeta1;
            //}
            //if(j == NY - 1){
            //    phi[index] = zeta2;
            //}
        }
    }
}

void Init_popu_poi(double* f0, double* f1, double* rho){
    for (size_t j = 0; j < NY; ++j) {
        for (size_t i = 0; i < NX; ++i) {
            const size_t index = scalar_index(i, j);
            const size_t f_index = fQ5_index(i, j, 0);
            double rho_ = rho[index];
            double feq[Q5];
            feqD2Q5_Poi(feq, rho_);
            for(size_t k = 0; k < Q5; ++k){
                f0[f_index + k] = feq[k];
                f1[f_index + k] = feq[k];
            }
        }
    }
}
void StreamingQ5(const size_t i, const size_t j, double* f0, const double* f1) {
    for(size_t k = 0; k < Q5; ++k) {
        const int ip = (NX + i - C_xQ5[k]) % NX;
        const int jp = (NY + j - C_yQ5[k]) % NY;
        const size_t f_pop_index = fQ5_index(ip, jp, k);
        f0[k] = f1[f_pop_index];
    }
}

inline void CollisionQ5_Poi(double* f1, double* f0, 
                const double R_, const double rho,
                const double omega_, const double dt_) {
    double feq[Q5];
    feqD2Q5_Poi(feq, rho);
    for(size_t k = 0; k < Q5; ++k) {
        double OMEGA = dt_ * w_h[k] * R_ * D_phi;
        //double OMEGA = dt * w_h[k] * R;
        f1[k] = f0[k] - omega_ * (f0[k] - feq[k]) + OMEGA;
    }
}

inline void CollisionQ5_Poi_R(double* f1, const double* f0, 
                        const double R_, const double rho_, 
                        const double omega_, const double dt_) {
    double feq[Q5];
    feqD2Q5_Poi(feq, rho_);
    double g1 = (f0[1] - feq[1]) - (f0[3] - feq[3]);
    double g2 = (f0[2] - feq[2]) - (f0[4] - feq[4]);
    for(size_t k = 0; k < Q5; ++k) {
        double OMEGA = dt_ * w_h[k] * R_ * D_phi;
        double fneq = w[k] * (C_xQ5[k] * g1 + C_yQ5[k] * g2) / alpha;
        f1[k] = feq[k] + (1 - omega_) * fneq + OMEGA;
    }
}
inline double Comp_Macro(double* f0, double G){
    double temp = 1.0 / (1.0 - w[0]) * (f0[1] + f0[2] + f0[3] + f0[4]);
    // double temp = 1.0 / (1.0 - w[0]) * (f0[1] + f0[2] + f0[3] + f0[4])
    //     + w_h[0] / (1.0 - w[0]) * tau_phi * dt * G;
    return temp;
}
void poisson(double* f0, double* f1, double* phi, 
             const double* Cp, const double* Cn, 
             const int errflag, double& errsum){
    for (size_t j = 0; j < NY; ++j) {
        for (size_t i = 0; i < NX; ++i) {
            const size_t index = scalar_index(i, j);
            size_t f_index = fQ5_index(i, j, 0);
            //double R = k * sinh(kk * phi[index] * C_phi) / C_phi * C_l * C_l;
            double R = -1.0 * (Cp[index] - Cn[index]) * C2q_to_eps;
            //double R = 0;
            double temp = Comp_Macro(f0 + f_index, R);
            if(errflag == 1){
                errsum += fabs(temp - phi[index]);
            }
            phi[index] = temp;

            CollisionQ5_Poi(f1 + f_index, f0 + f_index, R, phi[index], omega_phi, dt_phi);
            // CollisionQ5_Poi_R(f1 + f_index, f0 + f_index, R, phi[index], omega_phi, dt_phi);
        }
    }
    #pragma omp for collapse(2) schedule(static,8)
    for (size_t j = 0; j < NY; ++j) {
        for (size_t i = 0; i < NX; ++i) {
            size_t f_index = fQ5_index(i, j, 0);
            StreamingQ5(i, j, f0 + f_index, f1);
        }
    }

    // -------- Boundary condiction of Poisson equation -----------
    for (size_t i = 0; i < NX; ++i){
        // --------- Inamuro Boundary condiction ---------------
        // const size_t i_B = scalar_index(i, 0);
        // const size_t f_B = fQ5_index(i, 0, 0);
        // BottomBC_Inamuro_D_Poi(f0 + f_B, k_delta_phi * Phi);
        // const size_t i_T = scalar_index(i, NY - 1);
        // const size_t f_T = fQ5_index(i, NY - 1, 0);
        // TopBC_Inamuro_D_Poi(f0 + f_T, 0);

        // ---------- NBE Boundary scheme -----------
        double feq0[Q5], feq1[Q5];

        const size_t j_B1 = scalar_index(i, 1);
        const size_t jf_B1 = fQ5_index(i, 1, 0);
        const size_t jf_B0 = fQ5_index(i, 0, 0);
        feqD2Q5_Poi(feq1, phi[j_B1]);
        feqD2Q5_Poi(feq0, k_delta_phi * Phi);
        NEE_Q5(f0 + jf_B0, f0 + jf_B1, feq0, feq1);

        const size_t j_T1 = scalar_index(i, NY - 2);
        const size_t jf_T1 = fQ5_index(i, NY - 2, 0);
        const size_t jf_T0 = fQ5_index(i, NY - 1, 0);
        feqD2Q5_Poi(feq1, phi[j_T1]);
        feqD2Q5_Poi(feq0, 0);
        NEE_Q5(f0 + jf_T0, f0 + jf_T1, feq0, feq1);
    }
    // for (size_t j = 0; j < NY; ++j){
    //     const size_t i_L = scalar_index(0, j);
    //     const size_t f_L    = fQ5_index(0, j, 0);
        // LeftBC_Inamuro_D_Poi(f0 + f_L, k_zeta * zeta2);

    //     const size_t i_R = scalar_index(NX - 1, j);
    //     const size_t f_R =    fQ5_index(NX - 1, j, 0);
    //     RightBC_Inamuro_D_Poi(f0 + f_R, k_zeta * zeta2);
    // }
    // double half_phi_BC = (zeta2 * k_zeta + 0) * 0.5;
    // BottomLeftBC_Inamuro_D_Poi( f0 + fQ5_index(0     , 0     , 0), half_phi_BC);
    // TopLeftBC_Inamuro_D_Poi(    f0 + fQ5_index(0     , NY - 1, 0), half_phi_BC);
    // BottomRightBC_Inamuro_D_Poi(f0 + fQ5_index(NX - 1, 0     , 0), half_phi_BC);
    // TopRightBC_Inamuro_D_Poi(   f0 + fQ5_index(NX - 1, NY - 1, 0), half_phi_BC);
}

void Comp_E(double* f0, double* Ex, double* Ey){
    for (size_t j = 0; j < NY; ++j) {
        for (size_t i = 0; i < NX; ++i) {
            double sumx = 0.0;
            double sumy = 0.0;
            const size_t index = scalar_index(i, j);
            size_t f_index = fQ5_index(i, j, 0);
            for(size_t k = 0; k < Q5; ++k){
                f_index = fQ5_index(i, j, k);
                sumx += C_xQ5[k] * f0[f_index];
                sumy += C_yQ5[k] * f0[f_index];
            }
            Ex[index] = 1.0 / tau_phi / dt_phi / alpha * sumx;
            Ey[index] = 1.0 / tau_phi / dt_phi / alpha * sumy;
            //Ex[index] = 500 / C_E;
        }
    }
    // for (size_t j = 0; j < NY; ++j){
    //     const size_t i0 = scalar_index(0, j);
    //     const size_t i1 = scalar_index(1, j);
    //     Ex[i0] = Ex[i1];
    //     const size_t iN0 = scalar_index(NX - 1, j);
    //     const size_t iN1 = scalar_index(NX - 2, j);
    //     Ex[iN0] = Ex[iN1];
    // }
}
void Efield4(double* phi, double* Ex, double* Ey, const double dx) {
    const double dy = dx;
    const double twelve = 1.0 / 12.0;
    const double twothr = 2.0 / 3.0;
    const double elesix = 11.0 / 6.0;
    const double thrtwo = 3.0 / 2.0;
    const double onethr = 1.0 / 3.0;

    for(size_t i = 0; i < NX; ++i){
        for(size_t j = 0; j < NY; ++j){
            size_t ip1 = (i + 1) % NX;
            size_t jp1 = (j + 1) % NY;
            size_t ip2 = (i + 2) % NX;
            size_t jp2 = (j + 2) % NY;
            size_t im1 = (NX + i - 1) % NX;
            size_t jm1 = (NY + j - 1) % NY;
            size_t im2 = (NX + i - 2) % NX;
            size_t jm2 = (NY + j - 2) % NY;
    //for(size_t i = 2; i < NX - 2; ++i){
    //    for(size_t j = 2; j < NY - 2; ++j){
    //        size_t ip1 = (i + 1);
    //        size_t jp1 = (j + 1);
    //        size_t ip2 = (i + 2);
    //        size_t jp2 = (j + 2);
    //        size_t im1 = (i - 1);
    //        size_t jm1 = (j - 1);
    //        size_t im2 = (i - 2);
    //        size_t jm2 = (j - 2);
            size_t in = scalar_index(i, j);
            size_t inL1 = scalar_index(im1, j  );
            size_t inR1 = scalar_index(ip1, j  );
            size_t inL2 = scalar_index(im2, j  );
            size_t inR2 = scalar_index(ip2, j  );
            size_t inU1 = scalar_index(i  , jp1);
            size_t inD1 = scalar_index(i  , jm1);
            size_t inU2 = scalar_index(i  , jp2);
            size_t inD2 = scalar_index(i  , jm2);
            Ex[in] = -(-twelve * phi[inR2] + twothr * phi[inR1] - twothr * phi[inL1] + twelve * phi[inL2]) / dx;
            Ey[in] = -(-twelve * phi[inU2] + twothr * phi[inU1] - twothr * phi[inD1] + twelve * phi[inD2]) / dy;
        }
    }

    for(size_t j = 0; j < NY; ++j){
        size_t j0 = scalar_index(0, j);
        size_t j1 = scalar_index(1, j);
        size_t j2 = scalar_index(2, j);
        size_t j3 = scalar_index(3, j);
        size_t j4 = scalar_index(4, j);
        Ex[j0] = -(-elesix * phi[j0] + 3 * phi[j1] - thrtwo * phi[j2] + onethr * phi[j3]) / dx;
        Ex[j1] = -(-elesix * phi[j1] + 3 * phi[j2] - thrtwo * phi[j3] + onethr * phi[j4]) / dx;


        size_t jb0 = scalar_index(NX - 1 - 0, j);
        size_t jb1 = scalar_index(NX - 1 - 1, j);
        size_t jb2 = scalar_index(NX - 1 - 2, j);
        size_t jb3 = scalar_index(NX - 1 - 3, j);
        size_t jb4 = scalar_index(NX - 1 - 4, j);

        Ex[jb0] = -(elesix * phi[jb0] - 3 * phi[jb1] + thrtwo * phi[jb2] -onethr * phi[jb3]) / dx;
        Ex[jb1] = -(elesix * phi[jb1] - 3 * phi[jb2] + thrtwo * phi[jb3] -onethr * phi[jb4]) / dx;
    }

    for(size_t i = 0; i < NX; ++i){
        size_t i0 = scalar_index(i, 0);
        size_t i1 = scalar_index(i, 1);
        size_t i2 = scalar_index(i, 2);
        size_t i3 = scalar_index(i, 3);
        size_t i4 = scalar_index(i, 4);
        Ey[i0] = -(-elesix * phi[i0] + 3 * phi[i1] - thrtwo * phi[i2] + onethr * phi[i3]) / dy;
        Ey[i1] = -(-elesix * phi[i1] + 3 * phi[i2] - thrtwo * phi[i3] + onethr * phi[i4]) / dy;

        size_t ib0 = scalar_index(i, NY - 1 - 0);
        size_t ib1 = scalar_index(i, NY - 1 - 1);
        size_t ib2 = scalar_index(i, NY - 1 - 2);
        size_t ib3 = scalar_index(i, NY - 1 - 3);
        size_t ib4 = scalar_index(i, NY - 1 - 4);
        Ey[ib0] = -(elesix * phi[ib0] - 3 * phi[ib1] + thrtwo * phi[ib2] -onethr * phi[ib3]) / dy;
        Ey[ib1] = -(elesix * phi[ib1] - 3 * phi[ib2] + thrtwo * phi[ib3] -onethr * phi[ib4]) / dy;
    }
}
void Poisson_solver(size_t &tot, double *f0, double *f1, double *phi, 
             const double *Cp, const double *Cn, 
             double *Ex, double *Ey, 
             double &phiall, double &rela_err){
    //Init_popu_Poi(f0, f1, phi);
    size_t tot_step = 0;
    bool errFlag;
    double errSum;
    //double rela_err;

    tot = 0;
    for(size_t n = 0; n < MAXSTEP_Poi; ++n){
        errSum = 0.0;
        errFlag = true;
        //if(n% 100 == 0){
        //    errFlag = true;
        //    errSum = 0.0;
        //}else {
        //        errFlag = false;
        //}
        poisson(f0, f1, phi, Cp, Cn, errFlag, errSum);
        //if (n % 100 == 0){
        //    if (n % 1000 ==0){
        //        phiall = Compute_all(phi, NX * NY);
        //        rela_err = errSum / abs_all(phi, NX * NY);
        //        //rela_err = errSum;
        //        cout << "      -------Poi_Solver----------------------------\n"
        //             << "       step = " << std::setw(6) << n << ", "
        //             << std::fixed  << std::setprecision(10)
        //             << "phiall = " << phiall << " "
        //             << "errSum = " << errSum << " "
        //             << "Iteration Error = " <<  rela_err << endl;
        //    }
        //}
        if(errFlag){
            phiall = Compute_all(phi, NX * NY);
            rela_err = errSum / abs_all(phi, NX * NY);
            if (rela_err < rela_max_Poi_err){
                // cout << "    | -------Poi_Solver----------------------------\n"
                //     << "    |   step: " << std::setw(6) << n << ", "
                //     << std::fixed  << std::setprecision(10)
                //     << "errSum: " << errSum << " "
                //     << "phiall: " << phiall << " "
                //     //<< "errSum: " << errSum << " "
                //     << "Iteration Error: " <<  rela_err << endl;
                break;
            }
        }
        tot_step = n;
    }
    tot = tot_step;
    //ofstream Exline("line_Ex.txt");
    //for(size_t i = 0; i < NX; ++i){
    //    size_t index_in = scalar_index(i, 4);
    //    Exline << std:: setw(8) << static_cast<double>(i) / NX << "  "
    //         << std::fixed << std::setprecision(7)
    //         << Ex[index_in] << " " << phi[index_in] << " "
    //         << endl;
    //}
    //Exline.close();
    // Comp_E(f0, Ex, Ey);
    Efield4(phi, Ex, Ey, dx); 
}