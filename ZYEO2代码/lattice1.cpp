#include "lattice1.h"
void Print_lattice_para(std::ostream& out){
    out << "/* ******** lattice units , non-dimensionalisation ********** */\n";
    out << "******dx           = " << dx      << endl;
    out << "******dt           = " << dt      << endl;
    out << "******c            = " << c      << endl;
    out << "******dx_nu        = " << dx_nu      << endl;
    out << "******dt_nu        = " << dt_nu      << endl;
    out << "******c_nu         = " << c_nu      << endl;
    out << "******dx_D         = " << dx_D      << endl;
    out << "******dt_D         = " << dt_D      << endl;
    out << "******c_D          = " << c_D      << endl;
    out << "---------------------------------------------------------" << endl;
    out << "******C_Chemical_potential  = " << C_Chemical_potential  << endl;
    out << "---------------------------------------------------------" << endl;
    out << "******rho0         = " << rho0      << endl;
    out << "******rho_in       = " << rho_inlet      << endl;
    out << "******rho_out      = " << rho_outlet    << endl;
    out << "******Dp           = " << Dp         << endl;
    out << "******Dn           = " << Dn         << endl;
    out << "******delta_t      = " << delta_t    << endl;
    out << "******delta_tDiff  = " << delta_tD << endl;
    out << "******t0<-->tDiff  = " << t_tD << endl;
    out << "******nu           = " << nu         << endl;
    out << "******Phi          = " << Phi         << endl;
    out << "******kp           = " << KP         << endl;
    out << "******Kn           = " << KN         << endl;
    out << "******C2q(C -> q)  = " << C2q        << endl;
    out << "******eps          = " << eps        << endl;
    out << "******delta        = " << delta      << endl;
    out << "******Ex0          = " << delta / eps      << endl;
    out << "******RT           = " << RT         << endl;
    out << "******F_star       = " << F_star      << endl;
    //out << "******phi_sur      = " << phi_sur     << endl;
    out << "******zeta         = " << zeta_p     << endl;
    out << "******temp_delta   = " << temp_delta     << endl;
    out << "******C_plus_sur   = " << C_plus_sur     << endl;
    out << "******C_minu_sur   = " << C_minu_sur     << endl;
    out << "******delta        = " << delta << endl;
    out << "******epsilon      =  " << eps << endl;
    out << "******velocity_eof =  " << U_eof << endl;
    out << "******velocity_p   =  " << u_poi_max << endl;
    out << "******u_poi   =  " << u_poi << endl;

    out << "/****************** conversion factor **********************/\n";
    out << "******C_l (m)        = " << C_l << endl;
    out << "******C_rho (kg/m^3) = " << C_rho << endl;
    out << "******C_t (s)        = " << C_t << endl;
    out << "******C_t_Cn (s)     = " << C_t_Cn << endl;
    out << "******C_t_nu (s)     = " << C_t_nu << endl;

    out << "******C_q (C)        = " << C_q << endl;
    out << "******C_n (mol)      = " << C_n << endl;
    out << "******C_Phi (V)      = " << C_Phi << endl;
    out << "******C_u (m/s)      = " << C_u << endl;
    out << "******C_Dp (m^2/s)   = " << C_Dp << endl;
    out << "******C_nu (m^2/s)   = " << C_nu << endl;
    out << "******C_C (C/m^3)    = " << C_C  << endl;
    out << "******C_K (m^2/(s*V) = " << C_K << endl;
    out << "******C_C2q (C/mol)  = " << C_C2q << endl;
    out << "******C_e (J = C*V)  = " << C_e << endl;
    out << "******C_e_elec (J = C*V)  = " << C_e_elec << endl;
    out << "******C_e_mach (J = C*V)  = " << C_e_mech << endl;
    out << "******C_E (V/m)           = " << C_E << endl;
    out << "******C_Press (Pa)        = " << C_Press << endl;
    out << "******C_RT (J/mol)        = " << C_RT << endl;

    out << "/********************* relaxtion time **********************/\n";
    out << "******tau_p = " << tau_p << endl;
    out << "******tau_n = " << tau_n << endl;
    out << "******tau_f = " << tau_f << endl;
    // out << "******omegaP_f  = " << omegaP_f  << endl;
    // out << "******omegaP_Cp = " << omegaP_Cp << endl;
    // out << "******omegaP_Cn = " << omegaP_Cn << endl;
    // out << "******omegaM_f  = " << omegaM_f  << endl;
    // out << "******omegaM_Cp = " << omegaM_Cp << endl;
    // out << "******omegaM_Cn = " << omegaM_Cn << endl;




    out << "/******************* Step of PRINT Message *****************/\n";
    out << "Mesh number -- NX x NY ---" << NX << " x " << NY << endl;
    out << "---Max_Iteration = (delta_tD) *  10   " << Max_Iteration << endl;
    out << "---NMESSAGE      = (delta_tD) * 0.001 " << NMESSAGE      << endl;
    out << "---NSAVE         = NMESSAGE * 10      " << NSAVE         << endl;
}
void Print_lattice_Relaxation_Matrix(std::ostream& out){
    out << "----------------------- Relaxation Matrix ----------------" << endl;
    out << std::fixed << std::setprecision(6);
    out << "+------------- NSE -------------+\n";
    out << "|" << S_init_nu[0] << " | " << S_init_nu[1] << " | " << S_init_nu[2] << " |" << endl;
    out << "|" << S_init_nu[3] << " | " << S_init_nu[4] << " | " << S_init_nu[5] << " |" << endl;
    out << "|" << S_init_nu[6] << " | " << S_init_nu[7] << " | " << S_init_nu[8] << " |" << endl;
    out << "+-------------------------------+\n";
    out << "+                               +\n";
    out << "+-------------------------------+\n";
    out << "+------------ ADE P ------------+\n";
    out << "|" << S_init_p[0] << " | " << S_init_p[1] << " | " << S_init_p[2] << " |" << endl;
    out << "|" << S_init_p[3] << " | " << S_init_p[4] << " | " << S_init_p[5] << " |" << endl;
    out << "|" << S_init_p[6] << " | " << S_init_p[7] << " | " << S_init_p[8] << " |" << endl;
    out << "+-------------------------------+\n";
    out << "+                               +\n";
    out << "+-------------------------------+\n";
    out << "+------------ ADE N ------------+\n";
    out << "|" << S_init_n[0] << " | " << S_init_n[1] << " | " << S_init_n[2] << " |" << endl;
    out << "|" << S_init_n[3] << " | " << S_init_n[4] << " | " << S_init_n[5] << " |" << endl;
    out << "|" << S_init_n[6] << " | " << S_init_n[7] << " | " << S_init_n[8] << " |" << endl;
    out << "+-------------------------------+\n";
    out << "-------------------- Relaxation Matrix END *** -------------" << endl;
}

void Print_parameter(std::ostream& out){
    out << "/------------Physical parameters ----------------/\n";
    out << "******height and width (um -> m) <--> " << H_p << " x " << W_p << endl;
    out << "******velocity (m/s)             <--> " << U_p << endl;
    out << "******permittivity (C/(V*m), eps)<--> " << epsilon_p << endl;
    out << "******consentration (molar/m^3)  <--> " << C_unit_p << endl;
    out << "******Surface charge C (C/m^2)   <--> " << delta_p << endl;
    out << "******Temperature T (K)          <--> " << T_p << endl;
    out << "******Diffusivity D+, D- (m^2/s) <--> " << Diff_P_p << ", " << Diff_N_p << endl;
    out << "******liquid density rho (kg/m^3)<--> " << rho0_p << endl;
    out << "******viscoity (m^2/s, Pa*s)     <--> " << nu_p << ", " << mu_p << endl;
    out << "******Voltage( V )               <--> " << V_p << endl;
    out << "******Debye length (m)           <--> " << lambda_D << endl;
    out << "******(+)ion mobility  m^2/(s*V) <--> " << KP_p << endl;
    out << "******(-)ion mobility  m^2/(s*V) <--> " << KN_p << endl;
    out << "******mol/m^3 C->q C/m^3 (C/mol) <--> " << C2q_p << endl;

    out << "\n/------------Characteristic parameters-------------/\n";
    out << "## characteristic length               l0 <---> " << l0 << endl;
    out << "## characteristic consentration        C0 <---> " << C0 << endl;
    out << "## characteristic Diffusivity       Diff0 <---> " << Diff0 << endl;
    out << "## characteristic times                t0 <---> " << t0 << endl;
    out << "## characteristic Deybe length  lambda_D0 <---> " << lambda_D0 << endl;
    out << "## characteristic velocity             U0 <---> " << U0 << endl;
    out << "## characteristic voltage            Phi0 <---> " << Phi0 << endl;

    out << "\n/----------- Dimensionless parameters -------------/\n";
    out << "** Peclet number (ratio of convective to diffusion)       Pe = " << Pe_p << endl;
    out << "** Schmidt number (Momentum diffusion to ionic diffusion) Sc = " << Sc_p << endl;
    out << "** Reynolds number (inertial forces to viscous forces)    Re = " << Re_p << endl;
}
// NS equation Boundary BC  -----Zou/He------
inline void BottomBC_NEBB(double* f0, const double dt_, const double c_,
                double& rho, double ux = 0.0, double uy = 0.0, 
                double fx = 0.0, double fy = 0.0) {
    
    const double rho_t = (f0[0] + f0[1] + f0[3] 
                + 2 * (f0[4] + f0[7] + f0[8]) - fy * dt_ / 2 / c_) 
                * c_ / (c_ - uy);
    
    double f1_3 = half * (f0[1] - f0[3]);

    f0[2] = f0[4] + (twothirds * rho_t * uy - sixth * fy * dt_) / c_;
    f0[5] = f0[7] + rho_t * (half * ux + sixth * uy) / c_ 
            - f1_3 - dt_ * (quarter * fx + sixth * fy) / c_;
    f0[6] = f0[8] - rho_t * (half * ux - sixth * uy) / c_ 
            + f1_3 + dt_ * (quarter * fx - sixth * fy) / c_;
}

inline void TopBC_NEBB(double* f0, const double dt_, const double c_, 
                double& rho, double ux = 0.0, double uy = 0.0, 
                double fx = 0.0, double fy = 0.0) {

    const double rho_t = (f0[0] + f0[1] + f0[3] 
                + 2 * (f0[2] + f0[5] + f0[6]) + fy * dt_ / 2 / c_) 
                * c_ / (c_ + uy);
    
    double f1_3 = half * (f0[1] - f0[3]);

    f0[4] = f0[2] - (twothirds * rho_t * uy - sixth * fy * dt_) / c_;
    f0[7] = f0[5] - rho_t * (half * ux + sixth * uy) / c_ 
            + f1_3 + dt_ * (quarter * fx + sixth * fy) / c_;
    f0[8] = f0[6] + rho_t * (half * ux - sixth * uy) / c_
            - f1_3 - dt_ * (quarter * fx - sixth * fy) / c_;
}

inline void LeftBC_NEBB(double* f0, const double dt_, const double c_, 
                double& rho, double ux = 0.0, double uy = 0.0, 
                double fx = 0.0, double fy = 0.0) {
    const double rho_t = (f0[0] + f0[2] + f0[4] 
                + 2 * (f0[3] + f0[6] + f0[7]) - fy * dt_ / 2 / c_) 
                * c_ / (c_ - ux);
    
    double f2_4 = half * (f0[2] - f0[4]);
    

    f0[1] = f0[3] + (twothirds * rho_t * ux - sixth * fx * dt_) / c_;
    f0[5] = f0[7] + rho_t * (half * uy + sixth * ux) / c_ - f2_4 
            - (sixth * fx + quarter * fy) * dt_ / c_;
    f0[8] = f0[6] - rho_t * (half * uy - sixth * ux) / c_ + f2_4 
            - (sixth * fx - quarter * fy) * dt_ / c_;
}

inline void RightBC_NEBB(double* f0, const double dt_, const double c_, 
                double& rho, double ux = 0.0, double uy = 0.0, 
                double fx = 0.0, double fy = 0.0) {
    const double rho_t = (f0[0] + f0[2] + f0[4] 
                + 2 * (f0[1] + f0[5] + f0[8]) + fy * dt_ / 2 / c_) 
                * c_ / (c_ + ux);
    
    double f2_4 = half * (f0[2] - f0[4]);
    

    f0[3] = f0[1] - (twothirds * rho_t * ux - sixth * fx * dt_) / c_;
    f0[7] = f0[5] - rho_t * (half * uy + sixth * ux) / c_ + f2_4 
            + (sixth * fx + quarter * fy) * dt_ / c_;
    f0[6] = f0[8] + rho_t * (half * uy - sixth * ux) / c_ - f2_4 
            + (sixth * fx - quarter * fy) * dt_ / c_;
}
inline void LeftBottomBC_NEBB(double* f0, const double dt_, const double c_, 
                double& rho_, double ux = 0.0, double uy = 0.0, 
                double fx = 0.0, double fy = 0.0, double c_inv = 1.0) {
    f0[1] = f0[3] + twothirds * c_inv * rho_ * ux;
    f0[2] = f0[4] + twothirds * c_inv * rho_ * uy;
    f0[5] = f0[7] +  one_sixth * c_inv * rho_ * (ux + uy);
    f0[6] = one_twelfth * c_inv * (-ux + uy);
    f0[8] = one_twelfth * c_inv * ( ux - uy);
    f0[0] = rho_ - f0[1] - f0[2] - f0[3] - f0[4] 
                 - f0[5] - f0[6] - f0[7] - f0[8];
}
inline void LeftTopBC_NEBB(double* f0, const double dt_, const double c_, 
                double& rho_, double ux = 0.0, double uy = 0.0, 
                double fx = 0.0, double fy = 0.0, double c_inv = 1.0) {
    f0[1] = f0[3] + twothirds * c_inv * rho_ * ux;
    f0[4] = f0[2] - twothirds * c_inv * rho_ * uy;
    f0[8] = f0[6] +  one_sixth * c_inv * rho_ * ( ux - uy);
    f0[5] = one_twelfth * c_inv * ( ux + uy);
    f0[7] = one_twelfth * c_inv * (-ux - uy);
    f0[0] = rho_ - f0[1] - f0[2] - f0[3] - f0[4] 
                 - f0[5] - f0[6] - f0[7] - f0[8];
}
inline void RightBottomBC_NEBB(double* f0, const double dt_, const double c_, 
                double& rho_, double ux = 0.0, double uy = 0.0, 
                double fx = 0.0, double fy = 0.0, double c_inv = 1.0) {
    f0[3] = f0[1] - twothirds * c_inv * rho_ * ux;
    f0[2] = f0[4] + twothirds * c_inv * rho_ * uy;
    f0[6] = f0[8] +  one_sixth * c_inv * rho_ * (-ux + uy);
    f0[5] = one_twelfth * c_inv * ( ux + uy);
    f0[7] = one_twelfth * c_inv * (-ux - uy);
    f0[0] = rho_ - f0[1] - f0[2] - f0[3] - f0[4] 
                 - f0[5] - f0[6] - f0[7] - f0[8];
}
inline void RightTopBC_NEBB(double* f0, const double dt_, const double c_, 
                double& rho_, double ux = 0.0, double uy = 0.0, 
                double fx = 0.0, double fy = 0.0, double c_inv = 1.0) {
    f0[3] = f0[1] - twothirds * c_inv * rho_ * ux;
    f0[4] = f0[2] - twothirds * c_inv * rho_ * uy;
    f0[7] = f0[5] +  one_sixth * c_inv * rho_ * (-ux - uy);
    f0[6] = one_twelfth * c_inv * (-ux + uy);
    f0[8] = one_twelfth * c_inv * ( ux - uy);
    f0[0] = rho_ - f0[1] - f0[2] - f0[3] - f0[4] 
                 - f0[5] - f0[6] - f0[7] - f0[8];
}

// NS equation Boundary BC  ---------HWBB--------
inline void BottomBC_HWBB(double* f0, const double* f1, const double ux, const double uy, const double rho){
    // const double sixth = 1.0 / 6.0;
    f0[2] = f1[4] - 2 * w9[4] * rho * (C_xQ9_nu[4] * ux + C_yQ9_nu[4] * uy) / Cs2_nu;
    f0[5] = f1[7] - 2 * w9[7] * rho * (C_xQ9_nu[7] * ux + C_yQ9_nu[7] * uy) / Cs2_nu;
    f0[6] = f1[8] - 2 * w9[8] * rho * (C_xQ9_nu[8] * ux + C_yQ9_nu[8] * uy) / Cs2_nu;
    // f0[2] = f1[4];
    // f0[5] = f1[7];// + sixth * ux;
    // f0[6] = f1[8];// - sixth * ux;
}
inline void TopBC_HWBB(double* f0, const double* f1, const double ux, const double uy, const double rho){
    const double sixth = 1.0 / 6.0;
    // f0[4] = f1[2];
    // f0[7] = f1[5] - sixth * ux;
    // f0[8] = f1[6] + sixth * ux;
    f0[4] = f1[2] - 2 * w9[2] * rho * (C_xQ9_nu[2] * ux + C_yQ9_nu[2] * uy) / Cs2_nu;
    f0[7] = f1[5] - 2 * w9[5] * rho * (C_xQ9_nu[5] * ux + C_yQ9_nu[5] * uy) / Cs2_nu;
    f0[8] = f1[6] - 2 * w9[6] * rho * (C_xQ9_nu[6] * ux + C_yQ9_nu[6] * uy) / Cs2_nu;
}
inline void LeftBC_HWBB(double* f0, const double* f1, const double ux, const double uy, const double rho){
    // const double sixth = 1.0 / 6.0;
    f0[1] = f1[3] - 2 * w9[3] * rho * (C_xQ9_nu[3] * ux + C_yQ9_nu[3] * uy) / Cs2_nu;
    f0[5] = f1[7] - 2 * w9[7] * rho * (C_xQ9_nu[7] * ux + C_yQ9_nu[7] * uy) / Cs2_nu;
    f0[8] = f1[6] - 2 * w9[6] * rho * (C_xQ9_nu[6] * ux + C_yQ9_nu[6] * uy) / Cs2_nu;
    // f0[1] = f1[3];
    // f0[5] = f1[7] - sixth * ux;
    // f0[8] = f1[6] + sixth * ux;
}
inline void RightBC_HWBB(double* f0, const double* f1, const double ux, const double uy, const double rho){
    // const double sixth = 1.0 / 6.0;
    // f0[3] = f1[1];
    // f0[7] = f1[5] - sixth * ux;
    // f0[6] = f1[8] + sixth * ux;
    f0[3] = f1[1] - 2 * w9[1] * rho * (C_xQ9_nu[1] * ux + C_yQ9_nu[1] * uy) / Cs2_nu;
    f0[7] = f1[5] - 2 * w9[5] * rho * (C_xQ9_nu[5] * ux + C_yQ9_nu[5] * uy) / Cs2_nu;
    f0[6] = f1[8] - 2 * w9[8] * rho * (C_xQ9_nu[8] * ux + C_yQ9_nu[8] * uy) / Cs2_nu;
}


// Convection Diffusion equation Boundary BC  -----Inamuro BC------

inline void BottomBC_Inamuro_D(double* f0, const double TD) {
    double rho_temp = (TD - (f0[0] + f0[1] + f0[3] 
                        + f0[4] + f0[7] + f0[8])) 
                        / (w9[2] + w9[5] + w9[6]);
    
    f0[2] = w9[2] * rho_temp;
    f0[5] = w9[5] * rho_temp;
    f0[6] = w9[6] * rho_temp;
}
inline void TopBC_Inamuro_D(double* f0, const double TD) {
    double rho_temp = (TD - (f0[0] + f0[1] + f0[3] 
                        + f0[2] + f0[5] + f0[6])) 
                        / (w9[4] + w9[7] + w9[8]);
    
    f0[4] = w9[4] * rho_temp;
    f0[7] = w9[7] * rho_temp;
    f0[8] = w9[8] * rho_temp;
}

inline void LeftBC_Inamuro_D(double* f0, double TD) {
    double rho_temp = (TD - (f0[0] + f0[2] + f0[4] 
                        + f0[3] + f0[7] + f0[6])) 
                        / (w9[1] + w9[5] + w9[8]);
    f0[1] = w9[1] * rho_temp;
    f0[5] = w9[5] * rho_temp;
    f0[8] = w9[8] * rho_temp;
}
inline void RightBC_Inamuro_D(double* f0, double TD) {
    double rho_temp = (TD - (f0[0] + f0[2] + f0[4] 
                        + f0[1] + f0[5] + f0[8])) 
                        / (w9[3] + w9[6] + w9[7]);
    f0[3] = w9[3] * rho_temp;
    f0[6] = w9[6] * rho_temp;
    f0[7] = w9[7] * rho_temp;
}
void Initialization(double* rho, double* ux, double* uy, double* Cp, double* Cn, double* Ex, double* Ey, double* V) {
    for (size_t j = 0; j < NY; ++j) {
        for (size_t i = 0; i < NX; ++i) {
            const size_t index = scalar_index(i, j);
            rho[index] = rho0;
            Cp [index] = CP_bulk;
            Cn [index] = CN_bulk;
            ux [index] = 0.0;
            uy [index] = 0.0;
            V [index]  = 0.0;
            Ex [index] = 0.0;
            Ey [index] = 0.0;
            if(j <= 5 * fabs(cos(2 * 3.14159 * 8 / NX * i)) ){
                ux[index] =  5e-5 * sin(i) * cos(j);
                uy[index] = -5e-5 * sin(j) * cos(i);
            }
            // if(j == 0 ){
            //    Cp[index] = C_plus_sur / C_C;
            // }
            // else if(j == NY -1){
            //    Cn[index] = C_minu_sur / C_C;
            // }
        }
    }
}

inline void meqD2Q9_ADE(double* meq, const double phi, const double ux, const double uy){
    const double C2 = phi * (ux * ux + uy * uy);
    const double Cxx = phi * ux * ux;
    const double Cyy = phi * uy * uy;
    const double Cxy = phi * ux * uy;
    meq[0] = phi;
    meq[1] = -2 * phi + 3 * C2;
    meq[2] = phi - 3 * C2;
    meq[3] = phi * ux;
    meq[4] = -phi * ux;
    meq[5] = phi * uy;
    meq[6] = -phi * uy;
    meq[7] = Cxx - Cyy;
    meq[8] = Cxy;
}

inline void compute_Rk(double* Rk, const double& R, const double& ux, const double& uy){
    Rk[0] = w9[0] * R;
    Rk[1] = w9[1] * R * (1.0 + ux * Cs2_D_inv);
    Rk[2] = w9[2] * R * (1.0 + uy * Cs2_D_inv);
    Rk[3] = w9[3] * R * (1.0 - ux * Cs2_D_inv);
    Rk[4] = w9[4] * R * (1.0 - uy * Cs2_D_inv);
    Rk[5] = w9[5] * R * (1.0 + ( ux + uy) * Cs2_D_inv);
    Rk[6] = w9[6] * R * (1.0 + (-ux + uy) * Cs2_D_inv);
    Rk[7] = w9[7] * R * (1.0 + (-ux - uy) * Cs2_D_inv);
    Rk[8] = w9[8] * R * (1.0 + ( ux - uy) * Cs2_D_inv);
}

void Collision_MRT_CDE(double* m1, const double * f0, 
                        const double* s, double& C, 
                        const double ux, const double uy, 
                        const double dt_){
    // double Rk[Q9];
    // compute_Rk(Rk, R, ux, uy);

    const double f0_ = f0[0]; // - dt_ * Rk[0] * 0.5;
    const double f1_ = f0[1]; // - dt_ * Rk[1] * 0.5;
    const double f2_ = f0[2]; // - dt_ * Rk[2] * 0.5;
    const double f3_ = f0[3]; // - dt_ * Rk[3] * 0.5;
    const double f4_ = f0[4]; // - dt_ * Rk[4] * 0.5;
    const double f5_ = f0[5]; // - dt_ * Rk[5] * 0.5;
    const double f6_ = f0[6]; // - dt_ * Rk[6] * 0.5;
    const double f7_ = f0[7]; // - dt_ * Rk[7] * 0.5;
    const double f8_ = f0[8]; // - dt_ * Rk[8] * 0.5;

    C = f0_ + f1_ + f2_ + f3_ + f4_ + f5_ + f6_ + f7_ + f8_ ; //+ dt_ * 0.5 * R;

    double meq[Q9];
    meqD2Q9_ADE(meq, C, ux, uy);

    const double m_0 = f0_ + f1_ + f2_ + f3_ + f4_ + f5_ + f6_ + f7_ + f8_;
    const double m_1 = -4  *  f0_ - f1_ - f2_ - f3_ - f4_ + 2 * f5_ + 2 * f6_ + 2 * f7_ + 2 * f8_;
    const double m_2 = 4 * f0_ - 2 * f1_ - 2 * f2_ - 2 * f3_ - 2 * f4_ + f5_ + f6_ + f7_ + f8_;
    const double m_3 = f1_ - f3_ + f5_ - f6_ - f7_ + f8_;
    const double m_4 = -2 * f1_ + 2 * f3_ + f5_ - f6_ - f7_ + f8_;
    const double m_5 = f2_ - f4_ + f5_ + f6_ - f7_ - f8_;
    const double m_6 = -2 * f2_ + 2 * f4_ + f5_ + f6_ - f7_ - f8_;
    const double m_7 = f1_ - f2_ + f3_ - f4_;
    const double m_8 = f5_ - f6_ + f7_ - f8_;

    // const double p_m_0 = m_0 - (m_0 - meq[0]) * s[0];
    const double p_m_0 = m_0;
    const double p_m_1 = m_1 - (m_1 - meq[1]) * s[1];
    const double p_m_2 = m_2 - (m_2 - meq[2]) * s[2];
    const double p_m_3 = m_3 - (m_3 - meq[3]) * s[3];// + (m_5 - meq[5]) * s[9];
    const double p_m_4 = m_4 - (m_4 - meq[4]) * s[4];
    const double p_m_5 = m_5 - (m_5 - meq[5]) * s[5];// + (m_3 - meq[3]) * s[10];
    const double p_m_6 = m_6 - (m_6 - meq[6]) * s[6];
    const double p_m_7 = m_7 - (m_7 - meq[7]) * s[7];
    const double p_m_8 = m_8 - (m_8 - meq[8]) * s[8];

    // double I_S2MR[Q9];
    // compute_IS2_M_Rk(I_S2MR, Rk, s);

    m1[0] = p_m_0 ; //+ dt_ * I_S2MR[0];
    m1[1] = p_m_1 ; //+ dt_ * I_S2MR[1];
    m1[2] = p_m_2 ; //+ dt_ * I_S2MR[2];
    m1[3] = p_m_3 ; //+ dt_ * I_S2MR[3];
    m1[4] = p_m_4 ; //+ dt_ * I_S2MR[4];
    m1[5] = p_m_5 ; //+ dt_ * I_S2MR[5];
    m1[6] = p_m_6 ; //+ dt_ * I_S2MR[6];
    m1[7] = p_m_7 ; //+ dt_ * I_S2MR[7];
    m1[8] = p_m_8 ; //+ dt_ * I_S2MR[8];
}

void post_collision_m2f(double* f1, const double* m1_, const size_t NX0, const size_t NY0){
    for(size_t j = 0; j < NY0; ++j) {
	    for(size_t i = 0; i < NX0; ++i) {
            const size_t findex = fQ9_index(i, j, 0);
            const double m0 = m1_[findex + 0];
            const double m1 = m1_[findex + 1];
            const double m2 = m1_[findex + 2];
            const double m3 = m1_[findex + 3];
            const double m4 = m1_[findex + 4];
            const double m5 = m1_[findex + 5];
            const double m6 = m1_[findex + 6];
            const double m7 = m1_[findex + 7];
            const double m8 = m1_[findex + 8];

            f1[findex + 0] =  one_ninth * (m0 - m1 + m2);
            f1[findex + 1] = -one_eighteenth * m2 - one_thirtysixth * m1 + one_fourth * m7 + one_sixth * m3 - one_sixth * m4 + one_ninth * m0;
            f1[findex + 2] = -one_eighteenth * m2 - one_thirtysixth * m1 - one_fourth * m7 + one_sixth * m5 - one_sixth * m6 + one_ninth * m0;
            f1[findex + 3] = -one_eighteenth * m2 - one_thirtysixth * m1 + one_fourth * m7 - one_sixth * m3 + one_sixth * m4 + one_ninth * m0;
            f1[findex + 4] = -one_eighteenth * m2 - one_thirtysixth * m1 - one_fourth * m7 - one_sixth * m5 + one_sixth * m6 + one_ninth * m0;
            f1[findex + 5] =  one_twelfth * m4 + one_twelfth * m6 + one_eighteenth * m1 + one_thirtysixth * m2 + one_fourth * m8 + one_sixth * m3 + one_sixth * m5 + one_ninth * m0;
            f1[findex + 6] = -one_twelfth * m4 + one_twelfth * m6 + one_eighteenth * m1 + one_thirtysixth * m2 - one_fourth * m8 - one_sixth * m3 + one_sixth * m5 + one_ninth * m0;
            f1[findex + 7] = -one_twelfth * m4 - one_twelfth * m6 + one_eighteenth * m1 + one_thirtysixth * m2 + one_fourth * m8 - one_sixth * m3 - one_sixth * m5 + one_ninth * m0;
            f1[findex + 8] =  one_twelfth * m4 - one_twelfth * m6 + one_eighteenth * m1 + one_thirtysixth * m2 - one_fourth * m8 + one_sixth * m3 - one_sixth * m5 + one_ninth * m0;
        }
    }
}

inline void meqD2Q9_NSE(double* meq, const double& rho, const double& ux, const double& uy){
    const double ux2 = ux * ux;
    const double uy2 = uy * uy;
    meq[0] = 1.0 * rho;
    meq[1] = rho * ( 3 * ux2 + 3 * uy2 - 2 * c_nu_2) * c_nu_2_inv;
    meq[2] = rho * (-3 * ux2 - 3 * uy2 + 1 * c_nu_2) * c_nu_2_inv;
    meq[3] = 1.0 * ( rho * ux) * c_nu_1_inv;
    meq[4] = 1.0 * (-rho * ux) * c_nu_1_inv;
    meq[5] = 1.0 * ( rho * uy) * c_nu_1_inv;
    meq[6] = 1.0 * (-rho * uy) * c_nu_1_inv;
    meq[7] = 1.0 * rho * (ux2 - uy2) * c_nu_2_inv;
    meq[8] = 1.0 * rho * ux * uy * c_nu_2_inv;
}
inline void m2fQ9(double* f, const double* m){
    const double m0 = m[0];
    const double m1 = m[1];
    const double m2 = m[2];
    const double m3 = m[3];
    const double m4 = m[4];
    const double m5 = m[5];
    const double m6 = m[6];
    const double m7 = m[7];
    const double m8 = m[8];
    f[0] = m0/9 - m1/9   + m2/9;
    f[1] = m0/9 - m1/36 - m2/18 + m3/6 - m4/6  + m7/4;
    f[2] = m0/9 - m1/36 - m2/18 + m5/6 - m6/6  - m7/4;
    f[3] = m0/9 - m1/36 - m2/18 - m3/6 + m4/6  + m7/4;
    f[4] = m0/9 - m1/36 - m2/18 - m5/6 + m6/6  - m7/4;
    f[6] = m0/9 + m1/18 + m2/36 - m3/6 - m4/12 + m5/6 + m6/12 - m8/4;
    f[5] = m0/9 + m1/18 + m2/36 + m3/6 + m4/12 + m5/6 + m6/12 + m8/4;
    f[7] = m0/9 + m1/18 + m2/36 - m3/6 - m4/12 - m5/6 - m6/12 + m8/4;
    f[8] = m0/9 + m1/18 + m2/36 + m3/6 + m4/12 - m5/6 - m6/12 - m8/4;
}

inline void f2mQ9(double* m, const double* f){

    const double f0 = f[0];
    const double f1 = f[1];
    const double f2 = f[2];
    const double f3 = f[3];
    const double f4 = f[4];
    const double f5 = f[5];
    const double f6 = f[6];
    const double f7 = f[7];
    const double f8 = f[8];
    m[0] = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
    m[1] = (-4*f0 - f1 - f2 - f3 - f4 + 2*f5 + 2*f6 + 2*f7 + 2*f8);
    m[2] = (4*f0 - 2*f1 - 2*f2 - 2*f3 - 2*f4 + f5 + f6 + f7 + f8);
    m[3] = (f1 - f3 + f5 - f6 - f7 + f8);
    m[4] = (-2*f1 + 2*f3 + f5 - f6 - f7 + f8);
    m[5] = (f2 - f4 + f5 + f6 - f7 - f8);
    m[6] = (-2*f2 + 2*f4 + f5 + f6 - f7 - f8);
    m[7] = (f1 - f2 + f3 - f4);
    m[8] = (f5 - f6 + f7 - f8);
}
inline void MRT_F(double *R, const double fx, const double fy, 
                  const double ux, const double uy, const double c2_nu) {
    // c_nu_2
    R[0] = 0;
    R[1] =  6.0 * (fx * ux + fy * uy) * c_nu_2_inv;
    R[2] = -6.0 * (fx * ux + fy * uy) * c_nu_2_inv;
    R[3] =  1.0 * fx * c_nu_1_inv;
    R[4] = -1.0 * fx * c_nu_1_inv;
    R[5] =  1.0 * fy * c_nu_1_inv;
    R[6] = -1.0 * fy * c_nu_1_inv;
    R[7] =  2.0 * (fx * ux - fy * uy) * c_nu_2_inv;
    R[8] =  1.0 * (fx * uy + fy * ux) * c_nu_2_inv;
}
void Post_Coll_Add_F(double *m1, const double fx, const double fy,
                     const double ux, const double uy, 
                     const double *s, const double c2_nu)
{
    double MR[Q9];
    MRT_F(MR, fx, fy, ux, uy, c2_nu);

    m1[0] += MR[0] * (2.0 - s[0]) * 0.5;
    m1[1] += MR[1] * (2.0 - s[1]) * 0.5;
    m1[2] += MR[2] * (2.0 - s[2]) * 0.5;
    m1[3] += MR[3] * (2.0 - s[3]) * 0.5;
    m1[4] += MR[4] * (2.0 - s[4]) * 0.5;
    m1[5] += MR[5] * (2.0 - s[5]) * 0.5;
    m1[6] += MR[6] * (2.0 - s[6]) * 0.5;
    m1[7] += MR[7] * (2.0 - s[7]) * 0.5;
    m1[8] += MR[8] * (2.0 - s[8]) * 0.5;
}
void Init_f(double* f0, double* f1, double* ux, double* uy, double* rho) {
    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = 0; j < NY; ++j) {

            const size_t index = scalar_index(i, j);
            const size_t fQ9 = fQ9_index(i, j, 0);

            double feqQ9[Q9], meqQ9[Q9];
            double rho_ = rho[index];
            double ux_ = ux[index];
            double uy_ = uy[index];

            meqD2Q9_NSE(meqQ9, rho_, ux_, uy_);
            m2fQ9(feqQ9, meqQ9);
            for(size_t k = 0; k < Q9; ++k){
                const size_t fQ9_k = fQ9 + k;
                f1[fQ9_k] = meqQ9[k];
                f0[fQ9_k] = feqQ9[k];
            }
        }
    }
}

void Init_D(double* fp0, double* fp1, double* fn0, double* fn1, 
            const double* Cp, const double* Cn, 
            const double* ux, const double* uy) {
    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = 0; j < NY; ++j) {

            const size_t index = scalar_index(i, j);
            const size_t fQ9 = fQ9_index(i, j, 0);

            double feqQ9_P[Q9], meqQ9_P[Q9];
            double feqQ9_N[Q9], meqQ9_N[Q9];

            double Cp_ = Cp[index];
            double Cn_ = Cn[index];
            double ux_ = ux[index];
            double uy_ = uy[index];

            meqD2Q9_ADE(meqQ9_P, Cp_, ux_, uy_);
            meqD2Q9_ADE(meqQ9_N, Cn_, ux_, uy_);
            m2fQ9(feqQ9_P, meqQ9_P);
            m2fQ9(feqQ9_N, meqQ9_N);

            for(size_t k = 0; k < Q9; ++k){
                const size_t fQ9_k = fQ9 + k;
                fp1[fQ9_k] = meqQ9_P[k];
                fp0[fQ9_k] = feqQ9_P[k];

                fn1[fQ9_k] = meqQ9_N[k];
                fn0[fQ9_k] = feqQ9_N[k];
            }
        }
    }
}

inline void CollisionQ9_MRT(double* m1, const double* m0, const double* meq0, const double* s, const double dt_){
    m1[0] = m0[0];
    m1[1] = m0[1] - 1.0* s[1] * (m0[1] - meq0[1]);
    m1[2] = m0[2] - 1.0* s[2] * (m0[2] - meq0[2]);
    m1[3] = m0[3]; 
    m1[4] = m0[4] - 1.0* s[4] * (m0[4] - meq0[4]);
    m1[5] = m0[5]; 
    m1[6] = m0[6] - 1.0* s[6] * (m0[6] - meq0[6]);
    m1[7] = m0[7] - 1.0* s[7] * (m0[7] - meq0[7]);
    m1[8] = m0[8] - 1.0* s[8] * (m0[8] - meq0[8]);
}

inline void StreamingQ9(size_t i, size_t j, double* f0, double* f1) {
    for(size_t k = 0; k < Q9; ++k) {
        int ip = (NX + i - C_xQ9[k]) % NX;
        int jp = (NY + j - C_yQ9[k]) % NY;
        size_t f_p_index = fQ9_index(ip, jp, k);
        f0[k] = f1[f_p_index];
    }
}

inline void Comp_Macro(const double* f0, double& rho, double& ux, double& uy){
    double rho_ = 0;
    double uy_ = 0;
    double ux_ = 0;
    for(size_t k = 0; k < Q9; ++k){
        rho_ += f0[k];
        ux_ += f0[k] * C_xQ9_nu[k];
        uy_ += f0[k] * C_yQ9_nu[k];
    }
    rho = rho_;
    ux = ux_ / rho_;
    uy = uy_ / rho_;
}

void NSE(double* f0, double* f1, 
         double* rho, double* ux, double* uy, 
         const double dt_nu_, const double* S_nu_,
         const double* Cp, const double* Cn, 
         const double* Ex, const double* Ey, 
         const bool err_flag, double& rela_err){

    double temp1{0.0}, temp2{0.0};

    double *f_post  = (double*) malloc(Mesh_Size_Population_Q9);
    for (size_t j = 0; j < NY; ++j) {
        for (size_t i = 0; i < NX; ++i) {

            const size_t index = scalar_index(i, j);
            const size_t fQ9 = fQ9_index(i, j, 0);

            // ------  Compute the conserved quantities delta_rho and u

            double rho_ = rho[index];
            double ux_ = ux[index];
            double uy_ = uy[index];
            // const double ForceX = 0;
            // const double ForceY = 0;
            const double rho_q = (Cp[index] - Cn[index]) * C2q;
            const double ForceY = (Ey[index]) * rho_q;
            const double ForceX = (Ex1 + Ex[index]) * rho_q;

            Comp_Macro(f0 + fQ9, rho_, ux_, uy_);

            ux_ += ForceX * dt_nu_ * 0.5 / rho_;
            uy_ += ForceY * dt_nu_ * 0.5 / rho_;


            // if(err_flag == 1){
            //     temp1 += (ux_ - ux[index]) * (ux_ - ux[index]) + (uy_ - uy[index]) * (uy_ - uy[index]);
            //     temp2 += (ux_) * (ux_) + (uy_) * (uy_);
            // }

            double m0[Q9];
            double meq0[Q9];
            f2mQ9(m0, f0 + fQ9);
            meqD2Q9_NSE(meq0, rho_, ux_, uy_);
            CollisionQ9_MRT(f1 + fQ9, m0, meq0, S_nu_, dt_nu_);
            // Post_Coll_Add_F(f1 + fQ9, ForceX, ForceY, ux_, uy_, S_nu_, c_nu_2);
            f1[fQ9 + 3] += ForceX * dt_nu * c_nu_1_inv ;
            f1[fQ9 + 5] += ForceY * dt_nu * c_nu_1_inv ;

            m2fQ9(f_post + fQ9, f1 + fQ9);
            rho[index] = rho_;
            ux[index]  = ux_ ;
            uy[index]  = uy_ ;
        }
    }
    for (size_t j = 0; j < NY; ++j) {
        for (size_t i = 0; i < NX; ++i) {
            const size_t index = scalar_index(i, j);
            const size_t fQ9 = fQ9_index(i, j, 0);
            StreamingQ9(i, j, f0 + fQ9, f_post);
        }
    }

    for (size_t i = 0; i < NX; ++i){
        
        const size_t i_B = scalar_index(i, 0);
        const size_t f_B = fQ9_index(i, 0, 0);

        const double rho_q_B = (Cp[i_B] - Cn[i_B]) * C2q;
        const double ForceY_B = (Ey[i_B]) * rho_q_B;
        const double ForceX_B = (Ex1 + Ex[i_B]) * rho_q_B;

        BottomBC_NEBB(f0 + f_B, dt_nu, c_nu, rho[i_B], 0.0, 0.00, ForceX_B, ForceY_B);
        // BottomBC_HWBB(f0 + f_B, f_post + f_B, 0.0, 0, rho[i_B]); 

        const size_t i_T = scalar_index(i, NY - 1);
        const size_t f_T = fQ9_index(i, NY - 1, 0);
        const double rho_q_T = (Cp[i_T] - Cn[i_T]) * C2q;
        const double ForceY_T = (Ey[i_T]) * rho_q_T;
        const double ForceX_T = (Ex1 + Ex[i_T]) * rho_q_T;

        TopBC_NEBB(f0 + f_T, dt_nu, c_nu, rho[i_T], 0.0, 0.00, ForceX_T, ForceY_T);
        // TopBC_HWBB(f0 + f_T, f_post + f_T, 0.1, 0, rho[i_T]);
    }
    // for (size_t j = 0; j < NY; ++j){
    //     const size_t i_L = scalar_index(0, j);
    //     const size_t f_L    = fQ9_index(0, j, 0);
    //     // LeftBC_NEBB(f0 + f_L, dt_nu, c_nu, rho[i_L]);
    //     LeftBC_HWBB(f0 + f_L, f_post + f_L, 0.0, 0.0, rho[i_L]);

    //     const size_t i_R = scalar_index(NX - 1, j);
    //     const size_t f_R =    fQ9_index(NX - 1, j, 0);
    //     // RightBC_NEBB(f0 + f_R, dt_nu, c_nu, rho[i_R], 0.0, 0.0);
    //     RightBC_HWBB(f0 + f_R, f_post + f_R, 0.0, 0.0, rho[i_R]);
    // }
    // const size_t fLT = fQ9_index(0     , NY - 1, 0);
    // const size_t fRT = fQ9_index(NX - 1, NY - 1, 0);
    // const size_t fLB = fQ9_index(0     , 0     , 0);
    // const size_t fRB = fQ9_index(NX - 1, 0     , 0);
    // const size_t iLT_neighbor = scalar_index(1     , NY - 2);
    // const size_t iRT_neighbor = scalar_index(NX - 2, NY - 2);
    // const size_t iLB_neighbor = scalar_index(1     , 1     );
    // const size_t iRB_neighbor = scalar_index(NX - 2, 1     );

    // LeftTopBC_NEBB(    f0 + fLT, dt_nu, c_nu, rho[iLT_neighbor], 0.0, 0.0, 0, 0, c_nu_1_inv);
    // LeftBottomBC_NEBB( f0 + fLB, dt_nu, c_nu, rho[iLB_neighbor], 0.0, 0.0, 0, 0, c_nu_1_inv);
    // RightBottomBC_NEBB(f0 + fRB, dt_nu, c_nu, rho[iRB_neighbor], 0.0, 0.0, 0, 0, c_nu_1_inv);
    // RightTopBC_NEBB(   f0 + fRT, dt_nu, c_nu, rho[iRT_neighbor], 0.0, 0.0, 0, 0, c_nu_1_inv);

    free(f_post);
}

void PN_ADE(double* f0_P, double* f1_P, double* f0_N, double* f1_N, 
            double* C_P, double* C_N, 
            const double* ux, const double* uy, 
            const double* phi, const double* Ex, const double* Ey, 
            const double *S_p, const double* S_n,
            const double Kp, const double Kn, const double Dp, const double Dn,
            const size_t NX0, const size_t NY0, const double dt_){


    double *fP_post  = (double*) malloc(Mesh_Size_Population_Q9);
    double *fN_post  = (double*) malloc(Mesh_Size_Population_Q9);
    post_collision_m2f(fN_post, f1_N, NX0, NY0);
    post_collision_m2f(fP_post, f1_P, NX0, NY0);

    #pragma omp for collapse(2) schedule(static,16) 

    for(size_t j = 0; j < NY0; ++j) {
		for(size_t i = 0; i < NX0; ++i) {
            // const size_t index = scalar_indexMB(i, j, Block0.NX);
            // const size_t findex = f_indexQ9MB(i, j, 0, NX0);
            const size_t index = scalar_index(i, j);
            const size_t findex = fQ9_index(i, j, 0);

            StreamingQ9(i, j, f0_N + findex, fN_post );
            StreamingQ9(i, j, f0_P + findex, fP_post );

            //double time_ = exp((1.0 - 2.0 * Pi * Pi * kappa_p) * C_t * time_step);
            //double R = time_ * C_t * ((sin(Pi * C_l * (x[i] + y[j]))) + Pi * (ux + uy) * cos(Pi * C_l * (x[i] + y[j])));
            // double R = 0.0;

            // post_StreamingQ9(f0 + findex, R, ux, uy, dt);
        }
    }

    for (size_t i = 0; i < NX; ++i){
        // --------------- Bottom Boundary -------------------
        const size_t i_B = scalar_index(i, 0);
        const size_t f_B = fQ9_index(i, 0, 0);

        double CP_Bottom = 2 * CP_bulk;
        // double CN_Bottom = CN_bulk * exp(phi_[index]);
        double CN_Bottom = C_Chemical_potential;

        BottomBC_Inamuro_D(f0_P + f_B, CP_Bottom);
        BottomBC_Inamuro_D(f0_N + f_B, CN_Bottom);

        // -------------- Top Boundary --------------------
        const size_t i_T = scalar_index(i, NY - 1);
        const size_t f_T = fQ9_index(i, NY - 1, 0);

        double CP_Top = CP_bulk;
        // double CN_Bottom = CN_bulk * exp(phi_[index]);
        double CN_Top = CN_bulk;

        TopBC_Inamuro_D(f0_P + f_T, CP_Top);
        TopBC_Inamuro_D(f0_N + f_T, CN_Top);
    }
    // for (size_t j = 0; j < NY ; ++j){
    //     // --------------- Left Boundary -------------------
    //     const size_t i_L = scalar_index(0, j);
    //     const size_t f_L    = fQ9_index(0, j, 0);
    //     double CP_Left = CP_bulk;
    //     double CN_Left = CN_bulk;
    //     LeftBC_Inamuro_D(f0_P + f_L, CP_Left);
    //     LeftBC_Inamuro_D(f0_N + f_L, CN_Left);

    //     // --------------- Right Boundary -------------------
    //     const size_t i_R = scalar_index(NX - 1, j);
    //     const size_t f_R =    fQ9_index(NX - 1, j, 0);
    //     // double CP_Right = 2 * CP_bulk;
    //     // double CN_Right = C_Chemical_potential;
    //     double CP_Right = CP_bulk;
    //     double CN_Right = CN_bulk;                                                                                                                                   
        
    //     RightBC_Inamuro_D(f0_P + f_R, CP_Right);
    //     RightBC_Inamuro_D(f0_N + f_R, CN_Right);
    // }

    for (size_t j = 0; j < NY; ++j) {
        for (size_t i = 0; i < NX; ++i) {

            const size_t index = scalar_index(i, j);
            const size_t fQ9 = fQ9_index(i, j, 0);

            // double* C_P_ = C_P + index;
            // double* C_N_ = C_N + index;
            double cp_temp = C_P[index];
            double cn_temp = C_N[index];

            double ux__ = *(ux + index);
            double uy__ = *(uy + index);

            double Ex__ = *(Ex + index);
            double Ey__ = *(Ey + index);

            Collision_MRT_CDE(f1_P + fQ9, f0_P + fQ9, 
                                S_p, cp_temp, 
                                (Kp * Ex__ + ux__), 
                                (Kp * Ey__ + uy__), dt_);

            Collision_MRT_CDE(f1_N + fQ9, f0_N + fQ9, 
                                S_n, cn_temp, 
                                (Kn * Ex__ + ux__), 
                                (Kn * Ey__ + uy__), dt_);

            C_P[index] = cp_temp;
            C_N[index] = cn_temp;
        }
    }
    free(fP_post);
    free(fN_post);
}


