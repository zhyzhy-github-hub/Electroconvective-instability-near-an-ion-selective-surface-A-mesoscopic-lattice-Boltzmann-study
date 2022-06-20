#include "appendix1.h"
void outputTec(int m, double *rho, double *ux, double *uy, double *Cp, double *Cn, double *V, double *Ex, double *Ey)
{
    std::ostringstream name;
    name << "NS_P_NP" << m << ".dat";
    std::ofstream out(name.str().c_str());
    out << "Title= \"NS_PNP Flow\"\n";
    out << "VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"rho\", \"charge_p\", \"charge_n\", \"Voltage\", \"Ex\", \"Ey\", \"UV\", \"EE\", \"rho_E\"\n";
    out << "ZONE T= \"BOX\",I=" << NX << ",J=" << NY << ",F=	POINT" << endl;

    for (size_t j = 0; j < NY; ++j)
    { // if it is changed, should change doubles_saved;
        for (size_t i = 0; i < NX; ++i)
        {
            const size_t index = scalar_index(i, j);
            out << std::fixed << std::scientific //std::setprecision(10)
                << static_cast<double>(i) * C_l << " " << static_cast<double>(j) * C_l
                << " " << ux[index] << " " << uy[index] << " " << rho[index]
                << " " << Cp[index] << " " << Cn[index] << " " << V[index]
                << " " << Ex[index] << " " << Ey[index]
                << " " << sqrt(ux[index] * ux[index] + uy[index] * uy[index])
                << " " << sqrt(Ex[index] * Ex[index] + Ey[index] * Ey[index])
                << " " << (Cp[index] - Cn[index])
                << endl;
            // << static_cast<double>(i) * C_l << " " << static_cast<double>(j) * C_l
            // << " " << ux[index] * C_u << " " << uy[index] * C_u << " " << rho[index] * C_rho / C_rho
            // << " " << Cp[index] * C_C << " " << Cn[index] * C_C << " " << V[index] * C_Phi
            // << " " << Ex[index] * C_E << " " << Ey[index] * C_E
            // << " " << sqrt(ux[index] * ux[index] + uy[index] * uy[index]) * C_u
            // << " " << sqrt(Ex[index] * Ex[index] + Ey[index] * Ey[index]) * C_E << endl;
        }
    }
    out.close();
    cout << "cout*****************************************\n";
}

void WriteBinary_All(const size_t step, const size_t NX_, const size_t NY_,
                     const double *ff0, const double *ff1,
                     const double *gp0, const double *gp1,
                     const double *gn0, const double *gn1,
                     const double *fV0, const double *fV1,
                     const double *rho, const double *ux, const double *uy,
                     const double *phi, const double *Ex, const double *Ey,
                     const double *Cp, const double *Cn)
{
    //std::ostringstream name;
    //name << m * C_t << "_Reatart.dat";
    //std::ofstream w(name.str().c_str( ), std::ios::binary);
    const size_t scalar = NX_ * NY_ * sizeof(double);
    const size_t dis_fun_Q5 = NX_ * NY_ * Q5 * sizeof(double);
    const size_t dis_fun_Q9 = NX_ * NY_ * Q9 * sizeof(double);
    std::ofstream w("Restart.dat1", std::ios::binary);
    if (!w)
    {
        cout << "!!!********Can't write file*****************!!!!!!\n";
        exit(-1);
    }
    w.write((char *)&step, sizeof(size_t));
    // ------------ distribution functions ----------
    w.write((char *)ff0, dis_fun_Q9);
    w.write((char *)ff1, dis_fun_Q9);
    w.write((char *)gp0, dis_fun_Q9);
    w.write((char *)gp1, dis_fun_Q9);
    w.write((char *)gn0, dis_fun_Q9);
    w.write((char *)gn1, dis_fun_Q9);
    w.write((char *)fV1, dis_fun_Q5);
    w.write((char *)fV0, dis_fun_Q5);
    // ----------- scalar -----------------------
    w.write((char *)ux, scalar);
    w.write((char *)uy, scalar);
    w.write((char *)rho, scalar);
    w.write((char *)Ex, scalar);
    w.write((char *)Ey, scalar);
    w.write((char *)phi, scalar);
    w.write((char *)Cp, scalar);
    w.write((char *)Cn, scalar);
    w.close();
}
void ReStartBinary_All(const size_t &step, const size_t NX_, const size_t NY_,
                       const double *ff0, const double *ff1,
                       const double *gp0, const double *gp1,
                       const double *gn0, const double *gn1,
                       const double *fV0, const double *fV1,
                       const double *rho, const double *ux, const double *uy,
                       const double *phi, const double *Ex, const double *Ey,
                       const double *Cp, const double *Cn)
{
    const size_t scalar = NX_ * NY_ * sizeof(double);
    const size_t dis_fun_Q5 = NX_ * NY_ * Q5 * sizeof(double);
    const size_t dis_fun_Q9 = NX_ * NY_ * Q9 * sizeof(double);

    std::ifstream r("Restart.dat1", std::ios::binary);
    if (!r)
    {
        cout << "Can't read file*****************!!!!!!\n";
        exit(-1);
    }
    r.read((char *)&step, sizeof(size_t));
    // ------------ distribution functions ----------
    r.read((char *)ff0, dis_fun_Q9);
    r.read((char *)ff1, dis_fun_Q9);
    r.read((char *)gp0, dis_fun_Q9);
    r.read((char *)gp1, dis_fun_Q9);
    r.read((char *)gn0, dis_fun_Q9);
    r.read((char *)gn1, dis_fun_Q9);
    r.read((char *)fV1, dis_fun_Q5);
    r.read((char *)fV0, dis_fun_Q5);
    // ----------- scalar -----------------------
    r.read((char *)ux, scalar);
    r.read((char *)uy, scalar);
    r.read((char *)rho, scalar);
    r.read((char *)Ex, scalar);
    r.read((char *)Ey, scalar);
    r.read((char *)phi, scalar);
    r.read((char *)Cp, scalar);
    r.read((char *)Cn, scalar);
}
