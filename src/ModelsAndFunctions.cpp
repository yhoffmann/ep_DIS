// DEFINING GLOBAL PARAMETERS //
const double R = std::sqrt(3.3); //GeVm1
const double alpha_em = 1.0/137.036; //unit 1
const double e = std::sqrt(4.0*M_PI*alpha_em); //unit 1
const double e_c = 2.0/3.0;
const double e_t = 2.0/3.0;
const double e_b = -1.0/3.0;
const double N_c = 3.0;
const double hbarc = 0.1973; //GeV fm
const double BG = 4.0*hbarc*hbarc; //fm2 https://physics.nist.gov/cgi-bin/cuu/Value?rp
const double sigma_0 = 2.0*M_PI*R*R; //GeVm2
const double Q_s_0 = 1.0; //GeV

const double m_Q_c = 1.275; //in GeV
const double A_c = 0.211; //in GeV3/2

double e_Q;
const double A_Q = A_c; // in GeV3/2
const double epsilon = 1.0; // in GeV  // not used, see epsilonFunc
//const double Q = 1.0;

// Unit conversion factors
const double fmToGeVm1 = 1/hbarc;
const double GeVTofmm1 = 1/hbarc;
const double GeVm1Tofm = hbarc;
const double fmm1ToGeV = hbarc;
const double fm2TonB = 1.0e7;


void SetQuarkFlavor (char Flavor)
{
    if (Flavor == 'c'){e_Q = e_c;}
    else if (Flavor == 't'){e_Q = e_t;}
    else if (Flavor == 'b'){e_Q = e_b;}
    //std::cout << "Quark is of flavor " << Flavor << " and with charge " << e_Q << std::endl;
}



// MODEL NAMESPACES //

namespace GBWModel
{
    double T (double bx, double by)
    {
        return 1.0 / (2.0*M_PI*R*R) * exp( -(bx*bx+by*by)*(fmToGeVm1*fmToGeVm1) / (2.0*R*R) ); // in GeV2
    }

    double Q_s_squared (double bx, double by)
    {
        return sigma_0/*GeVm2*/ * Q_s_0/*GeV*/ * Q_s_0/*GeV*/ * T(bx, by)/*GeV2*/; // in GeV2
    }

    double dsigma_dip_d2b (double bx, double by, double rx, double ry)
    {
        return 2.0 * ( 1.0 - exp(-0.25* Q_s_squared(bx,by)/*GeV2*/ * (rx*rx+ry*ry)*(fmToGeVm1*fmToGeVm1)/*GeVm2*/) ); // unit 1, argument of exp also unit 1
    }
}

double epsilonFunc(double Q, double z)
    {
        return std::sqrt(Q*Q*z*(1.0-z)+m_Q_c*m_Q_c); // in GeV
    }

namespace PsiPsi_delta_evaluated
{
    double Trans (double Q, double rx, double ry, double z)
    {
        if (rx == 0.0) {rx = 1e-30;} // <-- without this GSL produces an error because Suave evaluates at xx=0 and gsl_sf_bessel_Kn(0,x) is only defined for x>0
        if (ry == 0.0) {ry = 1e-30;}
        return -A_Q/*GeV3/2*/ * std::sqrt(2.0*m_Q_c*N_c)/*GeV1/2*/ * e_Q/*1*/ * e/*1*/ *  gsl_sf_bessel_K0( epsilonFunc(Q,z)/*GeV*/ * std::sqrt(rx*rx+ry*ry)*fmToGeVm1/*GeVm1*/); // in GeV2
    }

    double Longi (double Q, double rx, double ry, double z)
    {
        return 2.0 / m_Q_c * Q * z * (1.0-z) * Trans(Q,rx,ry,z); // in GeV2
    }

    double Trans (double Q, double rx, double ry)
    {
        if (rx == 0.0) {rx = 1e-30;} // <-- without this GSL produces an error because Suave evaluates at xx=0 and gsl_sf_bessel_Kn(0,x) is only defined for x>0
        if (ry == 0.0) {ry = 1e-30;}
        return -A_Q * std::sqrt(2.0*m_Q_c*N_c) * e_Q * e *  gsl_sf_bessel_K0(epsilonFunc(Q,0.5) * std::sqrt(rx*rx+ry*ry));  // in GeV2 // z = 0.5 for now because of delta function
    }

    double Longi (double Q, double rx, double ry)
    {
        return 2.0 / m_Q_c * Q * Trans(Q,rx,ry); // in GeV2 
    }
}


// in integration: GeV2 -> GeV2 * fm4 = GeV2 * GeVm4 = GeVm2 and then squared after to give GeVm4

namespace A_coherent_Integrand_Function
{
    double Trans_real (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        return -PsiPsi_delta_evaluated::Trans(Q,rx,ry,z) * sin(-1.0*(bx*Deltax+by*Deltay + (0.5-z)*(rx*Deltax+ry*Deltay))*fmToGeVm1) * GBWModel::dsigma_dip_d2b(bx,by,rx,ry); // in GeV2
    }

    double Trans_imag (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        return PsiPsi_delta_evaluated::Trans(Q,rx,ry,z) * cos(-1.0*(bx*Deltax+by*Deltay + (0.5-z)*(rx*Deltax+ry*Deltay))*fmToGeVm1) * GBWModel::dsigma_dip_d2b(bx,by,rx,ry); // in GeV2
    }

    double Longi_real (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        return -PsiPsi_delta_evaluated::Longi(Q,rx,ry,z) * sin(-1.0*(bx*Deltax+by*Deltay + (0.5-z)*(rx*Deltax+ry*Deltay))*fmToGeVm1) * GBWModel::dsigma_dip_d2b(bx,by,rx,ry); // in GeV2
    }

    double Longi_imag (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        return PsiPsi_delta_evaluated::Longi(Q,rx,ry,z) * cos(-1.0*(bx*Deltax+by*Deltay + (0.5-z)*(rx*Deltax+ry*Deltay))*fmToGeVm1) * GBWModel::dsigma_dip_d2b(bx,by,rx,ry); // in GeV2
    }
}





namespace A_coherent_integrand_Function_first_order
{
    double Trans_real (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        return -PsiPsi_delta_evaluated::Trans(Q,rx,ry,z)/*GeV2*/ * sin(-1.0*(bx*Deltax+by*Deltay)*fmToGeVm1) * 2.0 * (-0.25*Q_s_0*Q_s_0)/*GeV2*/ * (rx*rx+ry*ry)*(fmToGeVm1*fmToGeVm1)/*GeVm2*/ * exp( -(bx*bx+by*by)*(fmToGeVm1*fmToGeVm1) / (2.0*R*R) ); // in GeV2
    }

    double Trans_imag (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        return PsiPsi_delta_evaluated::Trans(Q,rx,ry,z) * cos(-1.0*(bx*Deltax+by*Deltay)*fmToGeVm1) * 2.0 * (-0.25*Q_s_0*Q_s_0) * (rx*rx+ry*ry)*(fmToGeVm1*fmToGeVm1)  * exp( -(bx*bx+by*by)*(fmToGeVm1*fmToGeVm1) / (2.0*R*R) );
    }

    double Longi_real (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        return -PsiPsi_delta_evaluated::Longi(Q,rx,ry,z) * sin(-1.0*(bx*Deltax+by*Deltay)*fmToGeVm1) * 2.0 * (-0.25*Q_s_0*Q_s_0) * (rx*rx+ry*ry)*(fmToGeVm1*fmToGeVm1)  * exp( -(bx*bx+by*by)*(fmToGeVm1*fmToGeVm1) / (2.0*R*R) );
    }

    double Longi_imag (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        return PsiPsi_delta_evaluated::Longi(Q,rx,ry,z) * cos(-1.0*(bx*Deltax+by*Deltay)*fmToGeVm1) * 2.0 * (-0.25*Q_s_0*Q_s_0) * (rx*rx+ry*ry)*(fmToGeVm1*fmToGeVm1)  * exp( -(bx*bx+by*by)*(fmToGeVm1*fmToGeVm1) / (2.0*R*R) );
    }
}




namespace A_incoherent_Integrand_Function
{
    double Trans_real (double Q, double bx, double by, double bbarx, double bbary, double rx, double ry, double rbarx, double rbary, double Deltax, double Deltay, double z, double zbar)
    {
        return cos(-1.0*(bx*Deltax+by*Deltay+(0.5-z)*(rx*Deltax+ry*Deltay)-bbarx*Deltax-bbary*Deltay-(0.5-zbar)*(rbarx*Deltax+rbary*Deltay))) * PsiPsi_delta_evaluated::Trans(Q,rx,ry,z) * GBWModel::dsigma_dip_d2b(bx,by,rx,ry)
            * PsiPsi_delta_evaluated::Trans(Q,rbarx,rbary,zbar) * GBWModel::dsigma_dip_d2b(bbarx,bbary,rbarx,rbary);
    }

    double Trans_imag (double Q, double bx, double by, double bbarx, double bbary, double rx, double ry, double rbarx, double rbary, double Deltax, double Deltay, double z, double zbar)
    {
        return sin(-1.0*(bx*Deltax+by*Deltay+(0.5-z)*(rx*Deltax+ry*Deltay)-bbarx*Deltax-bbary*Deltay-(0.5-zbar)*(rbarx*Deltax+rbary*Deltay))) * PsiPsi_delta_evaluated::Trans(Q,rx,ry,z) * GBWModel::dsigma_dip_d2b(bx,by,rx,ry)
            * PsiPsi_delta_evaluated::Trans(Q,rbarx,rbary,zbar) * GBWModel::dsigma_dip_d2b(bbarx,bbary,rbarx,rbary);
    }

    double Longi_real (double Q, double bx, double by, double bbarx, double bbary, double rx, double ry, double rbarx, double rbary, double Deltax, double Deltay, double z, double zbar)
    {
        return cos(-1.0*(bx*Deltax+by*Deltay+(0.5-z)*(rx*Deltax+ry*Deltay)-bbarx*Deltax-bbary*Deltay-(0.5-zbar)*(rbarx*Deltax+rbary*Deltay))) * PsiPsi_delta_evaluated::Longi(Q,rx,ry,z) * GBWModel::dsigma_dip_d2b(bx,by,rx,ry)
            * PsiPsi_delta_evaluated::Longi(Q,rbarx,rbary,zbar) * GBWModel::dsigma_dip_d2b(bbarx,bbary,rbarx,rbary);
    }

    double Longi_imag (double Q, double bx, double by, double bbarx, double bbary, double rx, double ry, double rbarx, double rbary, double Deltax, double Deltay, double z, double zbar)
    {
        return sin(-1.0*(bx*Deltax+by*Deltay+(0.5-z)*(rx*Deltax+ry*Deltay)-bbarx*Deltax-bbary*Deltay-(0.5-zbar)*(rbarx*Deltax+rbary*Deltay))) * PsiPsi_delta_evaluated::Longi(Q,rx,ry,z) * GBWModel::dsigma_dip_d2b(bx,by,rx,ry)
            * PsiPsi_delta_evaluated::Longi(Q,rbarx,rbary,zbar) * GBWModel::dsigma_dip_d2b(bbarx,bbary,rbarx,rbary);
    }
}


namespace Integrand_full
{
    double smallValue = 1.0e-30;
    double Trans_real (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        if (rx==0&&ry==0){rx=smallValue;ry=smallValue;}
        return 1.0/(2.0*M_PI) * A_Q * std::sqrt(2.0*m_Q_c*N_c) * e_Q * e * gsl_sf_bessel_K0(epsilonFunc(Q,z)*std::sqrt(rx*rx+ry*ry)) 
        * sin(bx*Deltax+by*Deltay+(0.5-z)*(rx*Deltax+ry*Deltay)) * (1.0-exp(-0.25*sigma_0*Q_s_0*Q_s_0/(2.0*M_PI*BG)*exp(-(bx*bx+by*by)/(2.0*BG))));
    }

    double Trans_imag (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        if (rx==0&&ry==0){rx=smallValue;ry=smallValue;}
        return 1.0/(2.0*M_PI) * A_Q * std::sqrt(2.0*m_Q_c*N_c) * e_Q * e * gsl_sf_bessel_K0(epsilonFunc(Q,z)*std::sqrt(rx*rx+ry*ry)) 
        * cos(bx*Deltax+by*Deltay+(0.5-z)*(rx*Deltax+ry*Deltay)) * (1.0-exp(-0.25*sigma_0*Q_s_0*Q_s_0/(2.0*M_PI*BG)*exp(-(bx*bx+by*by)/(2.0*BG))));
    }

    double Longi_real (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        if (rx==0&&ry==0){rx=smallValue;ry=smallValue;}
        return 1.0/(1.0*M_PI) * A_Q * std::sqrt(2.0/m_Q_c*N_c) * e_Q * e * Q * z * (z-1.0) * gsl_sf_bessel_K0(epsilonFunc(Q,z)*std::sqrt(rx*rx+ry*ry)) 
        * sin(bx*Deltax+by*Deltay+(0.5-z)*(rx*Deltax+ry*Deltay)) * (1.0-exp(-0.25*sigma_0*Q_s_0*Q_s_0/(2.0*M_PI*BG)*exp(-(bx*bx+by*by)/(2.0*BG))));
    }

    double Longi_imag (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        if (rx==0&&ry==0){rx=smallValue;ry=smallValue;}
        return 1.0/(1.0*M_PI) * A_Q * std::sqrt(2.0/m_Q_c*N_c) * e_Q * e * Q * z * (z-1.0) * gsl_sf_bessel_K0(epsilonFunc(Q,z)*std::sqrt(rx*rx+ry*ry)) 
        * cos(bx*Deltax+by*Deltay+(0.5-z)*(rx*Deltax+ry*Deltay)) * (1.0-exp(-0.25*sigma_0*Q_s_0*Q_s_0/(2.0*M_PI*BG)*exp(-(bx*bx+by*by)/(2.0*BG))));
    }
}