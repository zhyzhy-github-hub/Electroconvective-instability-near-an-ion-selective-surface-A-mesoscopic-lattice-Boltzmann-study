#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include "lattice1.h"
#include "poisson1.h"
#include "appendix1.h"

using namespace std;
bool Restart_Simulation = 0;
double k_zeta = 1.0;
double k_delta_phi = -40;
int FLOW = 0;
double C_Chemical_potential = exp(k_delta_phi * Phi0 / Phi0);

double Comp_I(size_t n, const double *Cp, const double *Cn,
              const double *ux, const double *uy,
              const double *phi, const double *Ey,
              // const double dx_,
              double &U_max, double *cu)
{
    size_t J_Start = 5;
    double J = 0.0;
    double J1 = 0.0;
    double J2 = 0.0;
    double u_abs_max = 0.0;
    double ux_abs_max = 0.0;
    double uy_abs_max = 0.0;
    for (size_t j = J_Start; j < NY - 1; ++j)
    {
        for (size_t i = 0; i < NX; ++i)
        {
            const size_t j1 = scalar_index(i, j);
            const size_t j2 = scalar_index(i, j - 1);
            const size_t j3 = scalar_index(i, j + 1);
            const size_t ju = scalar_index(i, j - J_Start);
            double u_abs = (ux[ju] * ux[ju] + uy[ju] * uy[ju]);
            double ux_abs = ux[ju] * ux[ju];
            double uy_abs = uy[ju] * uy[ju];

            u_abs_max = (u_abs > u_abs_max) ? u_abs : u_abs_max;
            ux_abs_max = (ux_abs > ux_abs_max) ? ux_abs : ux_abs_max;
            uy_abs_max = (uy_abs > uy_abs_max) ? uy_abs : uy_abs_max;

            // if (u_abs > u_abs_max)
            // {
            //     u_abs_max = u_abs;
            // }
            // if (ux_abs > ux_abs_max)
            // {
            //     ux_abs_max = ux_abs;
            // }
            // if (uy_abs > uy_abs_max)
            // {
            //     uy_abs_max = uy_abs;
            // }

            double J1_temp = -(Diff_P_p * (Cp[j3] - Cp[j2]) * 0.5 + KP_p * Cp[j1] * (-Ey[j1]) * C_Phi) * C_C / C_l + uy[j1] * Cp[j1] * C_C * C_u;
            double J2_temp = -(Diff_N_p * (Cn[j3] - Cn[j2]) * 0.5 + KN_p * Cn[j1] * (-Ey[j1]) * C_Phi) * C_C / C_l + uy[j1] * Cn[j1] * C_C * C_u;
            // double J1_temp = -(Diff_P_p * (Cp[j3] - Cp[j2]) * 0.5 * dx_ + KP_p * Cp[j1] * (phi[j1] - phi[j2]) * C_Phi) * C_C / C_l + uy[j1] * Cp[j1] * C_C * C_u;
            // double J2_temp = -(Diff_N_p * (Cn[j3] - Cn[j2]) * 0.5 * dx_ + KN_p * Cn[j1] * (phi[j1] - phi[j2]) * C_Phi) * C_C / C_l + uy[j1] * Cn[j1] * C_C * C_u;
            //double J1_temp = -(Diff_P_p * (Cp[j1] - Cp[j2]) + KP_p * Cp[j1] * (-Ey[j1]) * C_Phi) * C_C / C_l + uy[j1] * Cp[j1] * C_C * C_u;
            //double J2_temp = -(Diff_N_p * (Cn[j1] - Cn[j2]) + KN_p * Cn[j1] * (-Ey[j1]) * C_Phi) * C_C / C_l + uy[j1] * Cn[j1] * C_C * C_u;
            J += (J1_temp - J2_temp) * C_l;
            J1 += (J1_temp)*C_l;
            J2 += (J2_temp)*C_l;
        }
    }
    U_max = sqrt(u_abs_max) * C_u;
    double Ux_max = sqrt(ux_abs_max) * C_u;
    double Uy_max = sqrt(uy_abs_max) * C_u;
    // double cu[3];
    cu[0] = J / (NY - J_Start);
    cu[1] = J1 / (NY - J_Start);
    cu[2] = J2 / (NY - J_Start);
    //return (J );

    std::ofstream write("line_jp_jn_j_uMax.txt", std::ios::app);
    write << std::fixed << std::setprecision(12) << n * C_t;
    write << std::fixed << std::setprecision(16);
    write << " " << cu[1];
    write << " " << cu[2];
    write << " " << cu[0];
    write << " " << U_max;
    write << " " << Ux_max;
    write << " " << Uy_max;
    write << std::endl;

    // }
    // write << endl;
    write.close();

    return U_max;
    // return (J / (NY - J_Start));
}


int main(int argc, char **argv)
{

    cout << "C_Chemical_potential = " << C_Chemical_potential << endl;
    // // omp_set_num_threads(2);
    // int id;

    // #pragma omp parallel private(id)
    //     {
    //     id = omp_get_thread_num();
    //     printf("%d: Hello World!\n", id);
    //     }

    // *************** print some parameters about this parameters****************
    ofstream fout("AParameter.txt");
    Print_parameter(cout);
    Print_parameter(fout);
    Print_lattice_para(cout);
    Print_lattice_para(fout);
    Print_lattice_Relaxation_Matrix(cout);
    Print_lattice_Relaxation_Matrix(fout);
    fout.close();

    if (argc == 4)
    {
        k_delta_phi = std::atof(argv[1]);
        Restart_Simulation = std::atoi(argv[2]);
        FLOW = std::atoi(argv[3]);
    }
    else
    {
        cout << "Please input zeta(-) and Restart" << std::endl;
    }

    cout << "input zeta " << k_delta_phi << ", " << Restart_Simulation << ", " << FLOW << endl;

    C_Chemical_potential = exp(k_delta_phi * Phi0 / Phi0);
    cout << "C_Chemical_potential = " << C_Chemical_potential << endl;

    double *f0 = (double *)malloc(Mesh_Size_Population_Q9);
    double *f1 = (double *)malloc(Mesh_Size_Population_Q9);
    double *gp0 = (double *)malloc(Mesh_Size_Population_Q9);
    double *gp1 = (double *)malloc(Mesh_Size_Population_Q9);
    double *gn0 = (double *)malloc(Mesh_Size_Population_Q9);
    double *gn1 = (double *)malloc(Mesh_Size_Population_Q9);
    double *ux = (double *)malloc(Mesh_Size_Scalar);
    double *uy = (double *)malloc(Mesh_Size_Scalar);
    double *Ex = (double *)malloc(Mesh_Size_Scalar);
    double *Ey = (double *)malloc(Mesh_Size_Scalar);
    double *rho = (double *)malloc(Mesh_Size_Scalar);
    double *Cp = (double *)malloc(Mesh_Size_Scalar);
    double *Cn = (double *)malloc(Mesh_Size_Scalar);
    double *phi = (double *)malloc(Mesh_Size_Scalar);
    double *fV0 = (double *)malloc(Mesh_Size_Population_Q5);
    double *fV1 = (double *)malloc(Mesh_Size_Population_Q5);

    Initialization(rho, ux, uy, Cp, Cn, Ex, Ey, phi);
    Init_f(f0, f1, ux, uy, rho);
    Init_D(gp0, gp1, gn0, gn1, Cp, Cn, ux, uy);
    Init_popu_poi(fV0, fV1, phi);
    outputTec(0, rho, ux, uy, Cp, Cn, phi, Ex, Ey);

    size_t last_step = 1;
    if (Restart_Simulation)
    {
        ReStartBinary_All(last_step, NX, NY, f0, f1, gp0, gp1, gn0, gn1,
                          fV0, fV1, rho, ux, uy, phi, Ex, Ey, Cp, Cn);
    }

    size_t n;
    bool err_flag_NSE, err_flag_Poi;
    double rela_err_NSE;
    double errsum_Poi;
    double phi_all, phi_ite_err;
    size_t poi_n;

    for (n = last_step; n < 10000000; ++n)
    {
        if (n * C_t > 0.040)
            break;

        if (FLOW)
        {
            for (size_t i_NSE = 0; i_NSE < D_to_nu; ++i_NSE)
            {
                // //     NSE(f0, f1, rho, ux, uy, err_flag_NSE, rela_err_NSE);

                NSE(f0, f1, rho, ux, uy, dt_nu, S_init_nu,
                    Cp, Cn, Ex, Ey, err_flag_NSE, rela_err_NSE);
            }
        }
        PN_ADE(gp0, gp1, gn0, gn1, Cp, Cn, ux, uy, phi, Ex, Ey,
               S_init_p, S_init_n, KP, KN, Dp, Dn, NX, NY, dt_D);

        Poisson_solver(poi_n, fV0, fV1, phi,
                       Cp, Cn, Ex, Ey, phi_all, phi_ite_err);

        // ---------------------------------------------------
        // poisson(fV0, fV1, phi, Cp, Cn, err_flag_Poi, errsum_Poi);
        // ---------------------------------------------------

        if (n % 10 == 0)
        {
            double rho_all = Compute_all(rho, NX * NY);
            double uy_all = Compute_all(uy, NX * NY);
            double ux_all = Compute_all(ux, NX * NY);
            double Cp_all = Compute_all(Cp, NX * NY);
            double Cn_all = Compute_all(Cn, NX * NY);
            cout << "n = " << n << " , "
                 << "t = " << n * C_t << " ,"
                 << "ux_all = " << ux_all << " , "
                 << "uy_all = " << uy_all << " , "
                 << "rho_all = " << rho_all << " , "
                 << "Cp_all = " << Cp_all << " , "
                 << "Cn_all = " << Cn_all << " , "
                 << "phi_all = " << phi_all << " , "
                 << "poi_n   = " << poi_n << " , "
                 << endl;

            double umax;
            double cu[3];
            double current = Comp_I(n, Cp, Cn, ux, uy, phi, Ey, umax, cu);

            //** write_data(n, umax, current);

            if (n % 10000 == 0)
            {
                WriteBinary_All(n, NX, NY, f0, f1, gp0, gp1, gn0, gn1,
                                fV0, fV1, rho, ux, uy, phi, Ex, Ey, Cp, Cn);
            }
            if (n % (20000) == 0)
            {
                outputTec(n, rho, ux, uy, Cp, Cn, phi, Ex, Ey);
            }
        }
    }

    free(f0);
    free(f1);
    free(gp0);
    free(gp1);
    free(gn0);
    free(gn1);
    free(ux);
    free(uy);
    free(Ex);
    free(Ey);
    free(rho);
    free(Cp);
    free(Cn);
    free(phi);
    free(fV0);
    free(fV1);
    return 0;
}
