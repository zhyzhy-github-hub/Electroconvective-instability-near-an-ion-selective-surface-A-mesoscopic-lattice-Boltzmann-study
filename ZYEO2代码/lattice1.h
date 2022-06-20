#ifndef __LATTICE1_H
#define __LATTICE1_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>

using std::cout;
using std::endl;
using std::vector;

const double half = 1.0 / 2.0, sixth = 1.0 / 6.0;
const double twothirds = 2.0 / 3.0, quarter = 1 / 4.0;

const double one_ninth = 1.0 / 9.0;
const double one_sixth = 1.0 / 6.0;
const double one_fourth = 1.0 / 4.0;
const double one_thirtysixth = 1.0 / 36.0;
const double one_twelfth = 1.0 / 12.0;
const double one_eighteenth = 1.0 / 18.0;

/* ---------------------------------------------------------------------*/
/*************** Some invariant physical constant ***********************/
const double k_B = 1.380649e-23;          // Boltzmann constant ( J/K )
const double NA = 6.02214076e23;          // Avogadro constant ( 1/mol )
const double e = 1.6021766341e-19;        // Elementary charge ( C, Coulomb )
const double epsilon_0 = 8.854187817e-12; // permittivity of free space ( F/m, C/(V*m) )
const double F = 96485.3415;              //Faraday constant ( C/molar )
const double R = 8.3144621;               // Universal gas constant ( J/(K*mol) )
// const double F         = NA * e; //Faraday constant ( C/molar )
// const double R         = k_B * NA; // Universal gas constant ( J/(K*mol) )
const double e_math = 2.7182818284; //
/*----------------------------------------------------------------------*/

// ****************** define height and width (um - 1e-6 -> m)
const double H_p = 5 * 1e-6;
const double W_p = 0.008 * H_p;
// ****************** mesh ******************************
const unsigned int NY = 161;
// const unsigned int NX = NY / H_p * W_p + 1;
const unsigned int NX = 1 + (NY - 1) * 4;

// ****************** define velocity (m/s)
const double U_p = 0.0;
// ****************** define relative permittivity
const double epsilon_p = 6.95e-10; // C^2/Jm
const double epsilon_r_p = epsilon_p / epsilon_0;
//const double epsilon = 6.93e-10;

// ****************** define consentration (mM = molar/m^3, 1M = 1mol/L)
const double C_unit_p = 0.0001e3;
const double C_P = C_unit_p;
// ****************** define Surface charge C (C/m^2) Now we don't neet it
//const double delta_p = -0.001;
// ****************** define Temperature T (K)
const double T_p = 273 + 20;
// const double T_p = 273 + 0;
// ****************** define Diffusivity D+, D- (m^2/s)
const double Diff_P_p = 1.00e-9;
const double Diff_N_p = 1.00e-9;
// ****************** define liquid density rho (kg/m^3) (dilute solution)
const double rho0_p = 1000;
// ****************** define Kinetic and Dynamic viscoity (m^2/s, Pa*s)
const double nu_p = 1.0e-6;
const double mu_p = rho0_p * nu_p;

// ****************** define Voltage( V )
const double V_p = 0.0;
const double zeta_p = -50e-3;
const double zeta1_p = -50e-3;
const double zeta2_p = -50e-3;
// ****************** surface charge
const double temp_delta = sqrt(8 * NA * C_P * epsilon_p * k_B * T_p);
const double delta_p = temp_delta * sinh(e * zeta_p / 2 / k_B / T_p);

// ****************** Debye length  EDL, z=1, Potassium  Chloride Solution(KCl)
const double lambda_D = sqrt(epsilon_p * R * T_p / 2 / F / F / C_unit_p);
const double kappa_p = 1.0 / lambda_D;

// ******************define electric field intensity
const double Ex1_p = 0.0 * 500; // V/m
const double U_eof = epsilon_p * zeta_p * Ex1_p / mu_p * (1.0 / cosh(kappa_p * H_p / 2) - 1);

// ****************** ion mobility(+) ( m^2/(s*V) )
const double KP_p = Diff_P_p * e / k_B / T_p;
const double KN_p = -Diff_N_p * e / k_B / T_p;

// ******************  Conversion factor from ion concentration (mol/m^3) to charge density (C/m^3)
const double C2q_p = F;

// ******************  Pressure
//const double grad_P_p = 500000; // Pa/m
const double grad_P_p = 0 * 500000; // Pa/m
const double dP_p = grad_P_p * W_p; // Pascal
//const double dP_p =  0.01 * 100;// Pascal
//const double grad_P_p = dP_p / W_p; // Pa/m
const double u_poi_max = grad_P_p * H_p * H_p / 8 / mu_p;

// ################## characteristic length
const double l0 = H_p;

// ################## characteristic consentration
const double C0 = C_P;
// ################## characteristic Diffusivity
const double Diff0 = Diff_P_p;
// ################## characteristic times
const double t0 = lambda_D * l0 / Diff0;
// ################## characteristic Deybe length
const double lambda_D0 = lambda_D / l0;
// ################## characteristic voltage
const double Phi0 = R * T_p / F;
// ################## characteristic velocity
//const double U0 = Diff0 / l0; // diffusion dominated
const double U0 = epsilon_0 * Phi0 * Phi0 / mu_p / l0; // convective

// ****************** Dimensionless parameters
const double Pe_p = U_eof * l0 / Diff0;
const double Sc_p = nu_p / Diff0;
const double Re_p = Pe_p / Sc_p;
//const double Re_p = U0 * l0 / nu_p;

const unsigned int Q9 = 9;
const unsigned int Q5 = 5;
const size_t Mesh_Size = sizeof(double) * NX * NY;
const size_t Mesh_Size_Scalar = Mesh_Size;
const size_t Mesh_Size_Population_Q5 = sizeof(double) * NX * NY * Q5;
const size_t Mesh_Size_Population_Q9 = sizeof(double) * NX * NY * Q9;

const double Cs = 1 / sqrt(3), Cs2 = pow(Cs, 2);

const int C_xQ9[Q9]{0, 1, 0, -1, 0, 1, -1, -1, 1};
const int C_yQ9[Q9]{0, 0, 1, 0, -1, 1, 1, -1, -1};

const double w9c = 4.0 / 9.0, w9s = 1.0 / 9.0, w9d = 1.0 / 36.0;
const double w9[Q9]{w9c, w9s, w9s, w9s, w9s, w9d, w9d, w9d, w9d};
const int OppQ9[Q9]{0, 3, 4, 1, 2, 7, 8, 5, 6};

/* ********lattice units****************************** */
const double Pe = Pe_p, Sc = Sc_p, Re = Re_p;
//const double C_l = l0 / NX;
const double C_l = l0 / (NY - 1);
// Unit lattice units, everyone is 1 ***************************
const double dx = 1.0;
const double dy = dx;
const double dt = dx;
const double c = dx / dt;
const double rho0 = 1.0, C_rho = rho0_p / rho0;

// Different lattice time step and lattice acoustic velocity **********
// ----------------- Advection diffusion equation ----------------
const double dx_D = dx;
const double dt_D = dt;
const double c_D = dx_D / dt_D;
const double Cs_D = 1 * c_D / sqrt(3), Cs2_D = pow(Cs_D, 2);
const double Cs2_D_inv = 1.0 / Cs2_D;
const double p_D = 1 * c_D;

const double C_xQ9_D[Q9]{0, p_D, 0, -p_D, 0, p_D, -p_D, -p_D, p_D};
const double C_yQ9_D[Q9]{0, 0, p_D, 0, -p_D, p_D, p_D, -p_D, -p_D};
// ------------------ Navier-Stokes equations ----------------------
const int D_to_nu = 100;
const double dx_nu = dx;
const double dt_nu = dt / D_to_nu;
const double c_nu = dx_nu / dt_nu;
const double c_nu_2 = c_nu * c_nu;
const double c_nu_3 = c_nu_2 * c_nu;
const double c_nu_4 = c_nu_2 * c_nu_2;
const double c_nu_1_inv = 1.0 / c_nu;
const double c_nu_2_inv = 1.0 / c_nu_2;
const double c_nu_3_inv = 1.0 / c_nu_3;
const double c_nu_4_inv = 1.0 / c_nu_4;

const double Cs_nu = 1 * c_nu / sqrt(3), Cs2_nu = pow(Cs_nu, 2);
const double p_nu = 1 * c_nu;
const double C_xQ9_nu[Q9]{0, p_nu, 0, -p_nu, 0, p_nu, -p_nu, -p_nu, p_nu};
const double C_yQ9_nu[Q9]{0, 0, p_nu, 0, -p_nu, p_nu, p_nu, -p_nu, -p_nu};

// positive charges
const double tau_p = 0.9;
const double Dp = Cs2_D * (tau_p - 0.5) * dt_D;
const double C_Dp = Diff_P_p / Dp;
const double C_t = pow(C_l, 2) / C_Dp;
const double delta_t = 1.0 / C_t;
const double delta_tD = t0 / C_t;
const double t_tD = 1.0 / t0;
const double C_u = C_l / C_t;
// negative charges
const double Dn = Diff_N_p / C_Dp;
const double tau_n = Dn / Cs2_D / dt_D + 0.5;
const double C_Dn = Diff_N_p / Dn;
const double C_t_Cn = pow(C_l, 2) / C_Dn;
// fluid
// const double nu = 0.1 * (NY - 1) / 1000;
const double nu = Sc * Dp;
const double C_nu = nu_p / nu;
// const double tau_f = nu * 3 + 0.5;
const double tau_f = nu / (Cs2_nu * dt_nu) + 0.5;

const double C_t_nu = pow(C_l, 2) / C_nu;

//const double u_poi_max = grad_P_p * H_p * H_p / 8 / mu_p;
const double u_poi = u_poi_max / C_u;
const double grad_P = 8 * nu * u_poi / (NY - 1) / (NY - 1);
const double C_Press = C_rho * C_u * C_u;
const double dP = grad_P * (NX - 1);
//const double dP = dP_p / C_Press;
//const double grad_P = grad_P_p / (C_P / C_l);
//const double grad_P = dP / (NX - 1);
const double rho_outlet = rho0;
const double rho_inlet = rho_outlet + 3 * (NX - 1) * grad_P;
//const double rho_inlet = 3 * (NX - 1) * gradP + rho_outlet;

// consentration
const double CP_bulk = 1.0;
const double C_C = C0 / CP_bulk;
const double CN_bulk = C0 / C_C;
const double Phi = 1.0;
const double C_Phi = Phi0 / Phi;
const double C_K = C_l * C_l / C_t / C_Phi;
const double KP = KP_p / C_K;
const double KN = KN_p / C_K;
const double C_C2q = C_rho * pow(C_l, 2) / pow(C_t, 2) / C_Phi / C_C;

const double C2q = C2q_p / C_C2q;     // Faraday constant
const double C_n = C_C * pow(C_l, 3); // molar
const double C_q = C_C2q * C_n;       // charge C
const double C_Q = C_C * C_C2q;       // charge density C/m^3
const double C_e = C_q * C_Phi;       // energy Joule
const double C_e_elec = C_q * C_Phi;
const double C_e_mech = C_rho * pow(C_l, 5) / pow(C_t, 2);
const double C_eps = C_q / C_Phi / C_l; //
const double C_E = C_Phi / C_l;         // electric field strength

const double eps = epsilon_p / C_eps;
const double delta = delta_p / (C_Q * C_l);
const double zeta1 = zeta1_p / C_Phi;
const double zeta2 = zeta2_p / C_Phi;
const double zeta_delta = Phi0 / C_Phi;
const double Ex1 = Ex1_p / C_E;

//const double temp_delta = delta_p / sqrt(8 * NA * C0 * epsilon_p * k_B * T_p);
//const double phi_sur = asinh(temp_delta) * 2 * k_B * T_p / e;
const double C_plus_sur = C0 * exp(-zeta_p * e / k_B / T_p);
const double C_minu_sur = C0 * exp(zeta_p * e / k_B / T_p);

extern double k_delta_phi;
extern double C_Chemical_potential;

const double C_RT = C_e / C_n;
const double RT = R * T_p / C_RT;
const double F_star = F / C_C2q;

const double omega_f = 1.0 / tau_f;
const double omega_p = 1.0 / tau_p;
const double omega_n = 1.0 / tau_n;

const double Magic_f = 1.0 / 4.0;
const double Magic_Cp = 1.0 / 4.0;
const double Magic_Cn = 1.0 / 4.0;

// const double omegaP_f  = 1.0 / tau_f;
// const double omegaP_Cp = 1.0 / tau_p;
// const double omegaP_Cn = 1.0 / tau_n;
// const double omegaM_f  = 1.0 / (0.5 + Magic_f * omegaP_f * 2 / (2 - omegaP_f));
// const double omegaM_Cp = 1.0 / (0.5 + Magic_Cp * omegaP_Cp * 2 / (2 - omegaP_Cp));
// const double omegaM_Cn = 1.0 / (0.5 + Magic_Cn * omegaP_Cn * 2 / (2 - omegaP_Cn));

// ------------ ADE relaxtion parameters ------ matrix --------
const double sD_0 = 0.0;
const double SD_add = 0.5;
const double sD_1 = 1.00 + SD_add;
const double sD_2 = 1.00 + SD_add;
const double sD_4 = 1.00 + SD_add;
const double sD_6 = 1.00 + SD_add;
const double sD_7 = 1.00 + SD_add;
const double sD_8 = 1.00 + SD_add;
const double sD_3 = omega_p;
const double sD_5 = omega_p;
const double S_init_p[Q9] = {sD_0, sD_1, sD_2, omega_p, sD_4, omega_p, sD_6, sD_7, sD_8};
const double S_init_n[Q9] = {sD_0, sD_1, sD_2, omega_n, sD_4, omega_n, sD_6, sD_7, sD_8};

// ------------ NSE relaxtion parameters ------ matrix --------
const double s_nu = omega_f;
const double s_e = s_nu;
// const double s_q = 1.8;
const double s_q = 8 * (2 - s_nu) / (8 - s_nu);
// const double s_eps = 2.0 / (6 * nu + 1);
const double s_eps = s_nu;

const double snu_0 = 1.0;
const double snu_1 = s_e;
const double snu_2 = s_eps;
const double snu_3 = 1.0;
const double snu_4 = s_q;
const double snu_5 = 1.0;
const double snu_6 = s_q;
const double snu_7 = s_nu;
const double snu_8 = s_nu;
const double S_init_nu[Q9] = {snu_0, snu_1, snu_2, snu_3, snu_4,
                              snu_5, snu_6, snu_7, snu_8};

const size_t Max_Iteration = floor(delta_tD) * 10;
const size_t NMESSAGE = floor(delta_tD) * 0.001;
const size_t NSAVE = NMESSAGE * 10;

void Print_parameter(std::ostream &out);
void Print_lattice_para(std::ostream &out);
void Print_lattice_Relaxation_Matrix(std::ostream &out);

void Initialization(double *rho, double *ux, double *uy, double *Cp, double *Cn, double *Ex, double *Ey, double *V);
void Init_equilibrium(double *f0, double *f1, double *gp0, double *gp1, double *gn0, double *gn1, double *rho, double *ux, double *uy, double *Cp, double *Cn, double *Ex, double *Ey);
void outputTec(int m, vector<double> &X, vector<double> &Y, double *rho, double *ux, double *uy, double *Cp, double *Cn, double *V, double *Ex, double *Ey);
void outputTec(int m, double *rho, double *ux, double *uy, double *Cp, double *Cn, double *V, double *Ex, double *Ey);

void NSE(double *f0, double *f1, double *rho, double *ux, double *uy, double *Cp, double *Cn, double *Ex, double *Ey);
void NSE(double *f0, double *f1, double *rho, double *ux, double *uy, const bool err_flag, double &rela_err);

void NSE(double *f0, double *f1,
         double *rho, double *ux, double *uy,
         const double dt_nu_, const double *S_nu_,
         const double *Cp, const double *Cn,
         const double *Ex, const double *Ey,
         const bool err_flag, double &rela_err);

void PN_ADE(double *fp0, double *fp1, double *fn0, double *fn1, double *Cp, double *Cn, double *ux, double *uy, double *phi, double *Ex, double *Ey, const double omega_pp, const double omega_pm, const double omega_np, const double omega_nm, const double Kp, const double Kn, const double Dp, const double Dn, const double C_bulk);
void PN_ADE(double *fp0, double *fp1, double *fn0, double *fn1,
            double *Cp, double *Cn,
            const double *ux, const double *uy,
            const double *phi, const double *Ex, const double *Ey,
            const double *S_p, const double *S_n,
            const double Kp, const double Kn, const double Dp, const double Dn,
            const size_t NX0, const size_t NY0, const double dt_);

void WriteBinary(size_t step, double *rho, double *ux, double *uy, double *phi, double *Ex, double *Ey, double *Cp, double *Cn);
void Restart(size_t &step, double *rho, double *ux, double *uy, double *phi, double *Ex, double *Ey, double *Cp, double *Cn);

void Initialization_Perturb(double *rho, double *ux, double *uy);
void Init_f(double *f0, double *f1, double *ux, double *uy, double *rho);
void Init_D(double *fp0, double *fp1, double *fn0, double *fn1,
            const double *Cp, const double *Cn,
            const double *ux, const double *uy);

inline size_t scalar_index(unsigned int x, unsigned int y)
{
    return (NX * y + x);
}
inline size_t fQ5_index(unsigned int x, unsigned int y, unsigned int d)
{
    return (Q5 * (NX * y + x) + d);
}
inline size_t fQ9_index(unsigned int x, unsigned int y, unsigned int d)
{
    return (Q9 * (NX * y + x) + d);
}
inline double Compute_all(double *f, size_t max)
{
    double f_all{0};
    for (size_t i = 0; i < max; ++i)
    {
        f_all += f[i];
    }
    return f_all;
}
inline double abs_all(double *f, size_t max)
{
    double f_all{0};
    for (size_t i = 0; i < max; ++i)
    {
        f_all += fabs(f[i]);
    }
    return f_all;
}

#endif
