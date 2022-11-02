// DEFINING GLOBAL PARAMETERS //


const double C = 1.0;
const double R = 1.0; //in fm
const double m_Q_c = 1.275; //in GeV
const double A_c = 0.211; //in GeV^(2/3)
const double alpha_em = 1.0/137.036;
const double e = std::sqrt(4.0*M_PI*alpha_em);
const double e_c = 2.0/3.0;
const double e_t = 2.0/3.0;
const double e_b = -1.0/3.0;
const double e_Q = e_c; // <-- DEFINING QUARK FLAVOR
const double A_Q = std::pow(A_c, 2.0) * 4.0 * M_PI * std::pow(e_c, 2.0) * alpha_em / std::pow(m_Q_c, 2.0); // ONLY FOR CHARM
const double N_c = 3.0;
const double epsilon = 1.0;
const double Q = 1.0;
const double Deltax = 0.0;
const double Deltay = 0.0;





// MODEL NAMESPACES //

namespace GBWModel {
    
    double T(double bx, double by) {
        return exp( -(std::pow(bx, 2.0)+std::pow(by, 2.0))/2.0/std::pow(R, 2.0) );
    }

    double Q_s(double bx, double by) {
        return C * T(bx, by);
    }

    double dsigma_qq_d2b(double bx, double by, double rx, double ry) {
        return 2.0 * ( 1 - exp(-0.25*std::pow(Q_s(bx, by), 2.0) * (std::pow(rx, 2.0)+std::pow(ry, 2.0))) );
    }
}



namespace PsiPsiSimpleModel {

    double PsiPsi_T_no_delta(double Q_squared, double rx, double ry, double z) {
        if (rx == 0.0) {rx = 1e-30;} // <-- without this GSL produces an error because Suave evaluates at xx=0 and gsl_sf_bessel_Kn(0,x) is only defined for x>0
        if (ry == 0.0) {ry = 1e-30;}
        return -A_Q * std::pow(2.0*m_Q_c*N_c, 1.0/2.0) * e_Q * e *  gsl_sf_bessel_Kn(0, epsilon * std::pow((std::pow(rx,2.0)+std::pow(ry,2.0)),0.5));
    }

    double PsiPsi_L_no_delta(double Q_squared, double rx, double ry, double z) {
        return 2.0 / m_Q_c * Q_squared * z * (1.0-z) * PsiPsi_T_no_delta(Q_squared,rx,ry,z);
    }
}