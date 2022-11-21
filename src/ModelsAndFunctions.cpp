// DEFINING GLOBAL PARAMETERS //
const double R = 1.0; //in fm
const double alpha_em = 1.0/137.036;
const double e = std::sqrt(4.0*M_PI*alpha_em);
const double e_c = 2.0/3.0;
const double e_t = 2.0/3.0;
const double e_b = -1.0/3.0;
const double N_c = 3.0;
const double M_HBARC = 0.197;
const double B_G = 4.0*M_HBARC*M_HBARC;
const double sigma_0 = 2.0*M_PI*B_G;
const double Q_s_0 = 1.0;

const double m_Q_c = 1.275; //in GeV
const double A_c = 0.211; //in GeV^(2/3)

double e_Q;
const double A_Q = std::pow(A_c, 2.0) * 4.0 * M_PI * std::pow(e_c, 2.0) * alpha_em / std::pow(m_Q_c, 2.0); // ONLY FOR CHARM FOR NOW
const double epsilon = 1.0;
//const double Q = 1.0;


void SetQuarkFlavor(char Flavor){
    if(Flavor == 'c'){e_Q = e_c;}
    else if(Flavor == 't'){e_Q = e_t;}
    else if(Flavor == 'b'){e_Q = e_b;}
    //std::cout << "Quark is of flavor " << Flavor << " and with charge " << e_Q << std::endl;
}



// MODEL NAMESPACES //

namespace GBWModel {
    
    double T(double bx, double by) {
        return 1.0 / (2.0*M_PI*B_G) * exp( -(std::pow(bx, 2.0)+std::pow(by, 2.0)) / (2.0*B_G) );
    }

    double Q_s_squared(double bx, double by) {
        return sigma_0 * T(bx, by) * Q_s_0 * Q_s_0;
    }

    double dsigma_dip_d2b(double bx, double by, double rx, double ry) {
        return 2.0 * ( 1 - exp(-0.25* Q_s_squared(bx,by) * std::sqrt(std::pow(rx, 2.0)+std::pow(ry, 2.0))) );
    }
}



namespace PsiPsiSimpleModel {

    double PsiPsi_T_no_delta(double Q, double rx, double ry, double z) {
        if (rx == 0.0) {rx = 1e-30;} // <-- without this GSL produces an error because Suave evaluates at xx=0 and gsl_sf_bessel_Kn(0,x) is only defined for x>0
        if (ry == 0.0) {ry = 1e-30;}
        return -A_Q * std::pow(2.0*m_Q_c*N_c, 1.0/2.0) * e_Q * e *  gsl_sf_bessel_Kn(0, epsilon * std::sqrt(std::pow(rx,2.0)+std::pow(ry,2.0)));
    }

    double PsiPsi_L_no_delta(double Q, double rx, double ry, double z) {
        return 2.0 / m_Q_c * Q * z * (1.0-z) * PsiPsi_T_no_delta(Q,rx,ry,z);
    }
}


namespace A_coherent_Integrand_Function {

    std::complex<double> T(double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z) {
        return 1.0i * PsiPsiSimpleModel::PsiPsi_T_no_delta(Q,rx,ry,z) * exp(-1.0i*(bx*Deltax+by*Deltay + (0.5-z)*(rx*Deltax+ry*Deltay))) * GBWModel::dsigma_dip_d2b(bx,by,rx,ry);
    }

    std::complex<double> L(double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z) {
        return 1.0i * PsiPsiSimpleModel::PsiPsi_L_no_delta(Q,rx,ry,z) * exp(-1.0i*(bx*Deltax+by*Deltay + (0.5-z)*(rx*Deltax+ry*Deltay))) * GBWModel::dsigma_dip_d2b(bx,by,rx,ry);
    }
}

namespace A_incoherent_Integrand_Function {

    std::complex<double> T(double Q, double bx, double by, double bbarx, double bbary, double rx, double ry, double rbarx, double rbary, double Deltax, double Deltay, double z, double zbar) {
        return 1.0i * PsiPsiSimpleModel::PsiPsi_T_no_delta(Q,rx,ry,z) * exp(-1.0i*(bx*Deltax+by*Deltay + (0.5-z)*(rx*Deltax+ry*Deltay))) * GBWModel::dsigma_dip_d2b(bx,by,rx,ry)
        * (-1.0i) * PsiPsiSimpleModel::PsiPsi_T_no_delta(Q,rbarx,rbary,zbar) * exp(+1.0i*(bbarx*Deltax+bbary*Deltay + (0.5-z)*(rbarx*Deltax+rbary*Deltay))) * GBWModel::dsigma_dip_d2b(bbarx,bbary,rbarx,rbary);
    }

    std::complex<double> L(double Q, double bx, double by, double bbarx, double bbary, double rx, double ry, double rbarx, double rbary, double Deltax, double Deltay, double z, double zbar) {
        return 1.0i * PsiPsiSimpleModel::PsiPsi_L_no_delta(Q,rx,ry,z) * exp(-1.0i*(bbarx*Deltax+bbary*Deltay + (0.5-z)*(rbarx*Deltax+rbary*Deltay))) * GBWModel::dsigma_dip_d2b(bx,by,rx,ry)
        * (-1.0i) * PsiPsiSimpleModel::PsiPsi_L_no_delta(Q,rbarx,rbary,zbar) * exp(+1.0i*(bbarx*Deltax+bbary*Deltay + (0.5-z)*(rbarx*Deltax+rbary*Deltay))) * GBWModel::dsigma_dip_d2b(bbarx,bbary,rbarx,rbary);
    }
}