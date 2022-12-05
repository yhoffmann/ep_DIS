// DEFINING GLOBAL PARAMETERS //
const double R = 1.0; // = std::sqrt(3.3); //in fm = GeV-1
const double alpha_em = 1.0/137.036;
const double e = std::sqrt(4.0*M_PI*alpha_em);
const double e_c = 2.0/3.0;
const double e_t = 2.0/3.0;
const double e_b = -1.0/3.0;
const double N_c = 3.0;
const double M_HBARC = 0.197; //GeV fm
const double B_G = 4.0*M_HBARC*M_HBARC;
const double sigma_0 = 2.0*M_PI*B_G;
const double Q_s_0 = 1.0;

const double m_Q_c = 1.275; //in GeV
const double A_c = 0.211; //in GeV^(2/3)

double e_Q;
const double A_Q = A_c*A_c * 4.0 * M_PI * std::pow(e_c, 2.0) * alpha_em / std::pow(m_Q_c, 2.0); // ONLY FOR CHARM FOR NOW
const double epsilon = 1.0;
//const double Q = 1.0;

// STRUCT FOR PASSING PARAMETERS TO INTEGRAND FUNCTION // might need to add functionality for passing epsilon to PsiPsi functions
struct USERDATA
{
    double Q;
    double Deltax;
    double Deltay;
    double z;
};


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
        return 1.0 / (2.0*M_PI*B_G) * exp( -(bx*bx+by*by) / (2.0*B_G) );
        //return 1.0 / (2.0*M_PI*R) * exp( -(bx*bx+by*by) / (2.0*R) );
    }

    double Q_s_squared (double bx, double by)
    {
        return sigma_0 * Q_s_0 * Q_s_0 * T(bx, by);
    }

    double dsigma_dip_d2b (double bx, double by, double rx, double ry)
    {
        return 2.0 * ( 1.0 - exp(-0.25* Q_s_squared(bx,by) * (rx*rx+ry*ry)) );
    }
}



namespace PsiPsi_delta_evaluated
{
    double Trans (double Q, double rx, double ry, double z)
    {
        if (rx == 0.0) {rx = 1e-30;} // <-- without this GSL produces an error because Suave evaluates at xx=0 and gsl_sf_bessel_Kn(0,x) is only defined for x>0
        if (ry == 0.0) {ry = 1e-30;}
        return -A_Q * std::sqrt(2.0*m_Q_c*N_c) * e_Q * e *  gsl_sf_bessel_K0(epsilon * std::sqrt(rx*rx+ry*ry));
    }

    double Longi (double Q, double rx, double ry, double z)
    {
        return 2.0 / m_Q_c * Q * z * (1.0-z) * Trans(Q,rx,ry,z);
    }
}


namespace A_coherent_Integrand_Function
{
    double Trans_real (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        return -1.0 * PsiPsi_delta_evaluated::Trans(Q,rx,ry,z) * sin(-1.0*(bx*Deltax+by*Deltay + (0.5-z)*(rx*Deltax+ry*Deltay))) * GBWModel::dsigma_dip_d2b(bx,by,rx,ry);
    }

    double Trans_imag (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        return PsiPsi_delta_evaluated::Trans(Q,rx,ry,z) * cos(-1.0*(bx*Deltax+by*Deltay + (0.5-z)*(rx*Deltax+ry*Deltay))) * GBWModel::dsigma_dip_d2b(bx,by,rx,ry);
    }

    double Longi_real (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        return -1.0 * PsiPsi_delta_evaluated::Longi(Q,rx,ry,z) * sin(-1.0*(bx*Deltax+by*Deltay + (0.5-z)*(rx*Deltax+ry*Deltay))) * GBWModel::dsigma_dip_d2b(bx,by,rx,ry);
    }

    double Longi_imag (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        return PsiPsi_delta_evaluated::Longi(Q,rx,ry,z) * cos(-1.0*(bx*Deltax+by*Deltay + (0.5-z)*(rx*Deltax+ry*Deltay))) * GBWModel::dsigma_dip_d2b(bx,by,rx,ry);
    }
}


namespace dsigma_dt_coherent_analytical
{
    double Trans (double Q, double Deltax, double Deltay, double z)
    {
        std::complex<double> Amplitude = -1.0i * sigma_0 * Q_s_0 * Q_s_0 / (std::pow(epsilon,4.0)) * exp(-B_G * (Deltax*Deltax+Deltay*Deltay) / 2.0) * A_Q * std::sqrt(2.0*m_Q_c*N_c) * e_Q * e;
        return std::pow(std::abs(Amplitude),2.0) / 16.0 / M_PI * 1.0e7;
    }

    double Longi (double Q, double Deltax, double Deltay, double z)
    {
        std::complex<double> Amplitude = -1.0i * sigma_0 * Q_s_0 * Q_s_0 / (std::pow(epsilon,4.0)) * exp(-B_G * (Deltax*Deltax+Deltay*Deltay) / 2.0) * 2.0 * z * (1.0-z) * A_Q * std::sqrt(2.0*N_c/m_Q_c) * e_Q * e * Q;
        return std::pow(std::abs(Amplitude),2.0) / 16.0 / M_PI * 1.0e7;
    }
}


namespace A_coherent_integrand_Function_first_order
{
    double Trans_real (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        return -PsiPsi_delta_evaluated::Trans(Q,rx,ry,z) * sin(-1.0*(bx*Deltax+by*Deltay + (0.5-z)*(rx*Deltax+ry*Deltay))) * (-0.25*sigma_0*Q_s_0*Q_s_0/2.0/M_PI/B_G) * (rx*rx+ry*ry) * exp(-(bx*bx+by*by)/2.0/B_G);
    }

    double Trans_imag (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        return PsiPsi_delta_evaluated::Trans(Q,rx,ry,z) * cos(-1.0*(bx*Deltax+by*Deltay + (0.5-z)*(rx*Deltax+ry*Deltay))) * (-0.25*sigma_0*Q_s_0*Q_s_0/2.0/M_PI/B_G) * (rx*rx+ry*ry) * exp(-(bx*bx+by*by)/2.0/B_G);
    }

    double Longi_real (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        return -PsiPsi_delta_evaluated::Longi(Q,rx,ry,z) * sin(-1.0*(bx*Deltax+by*Deltay + (0.5-z)*(rx*Deltax+ry*Deltay))) * (-0.25*sigma_0*Q_s_0*Q_s_0/2.0/M_PI/B_G) * (rx*rx+ry*ry) * exp(-(bx*bx+by*by)/2.0/B_G);
    }

    double Longi_imag (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        return PsiPsi_delta_evaluated::Longi(Q,rx,ry,z) * cos(-1.0*(bx*Deltax+by*Deltay + (0.5-z)*(rx*Deltax+ry*Deltay))) * (-0.25*sigma_0*Q_s_0*Q_s_0/2.0/M_PI/B_G) * (rx*rx+ry*ry) * exp(-(bx*bx+by*by)/2.0/B_G);
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

double smallValue = 1.0e-2;
namespace Integrand_full
{
    double Trans_real (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        if (rx==0&&ry==0){rx=smallValue;ry=smallValue;}
        return 1.0/(2.0*M_PI) * A_Q * std::sqrt(2.0*m_Q_c*N_c) * e_Q * e * gsl_sf_bessel_K0(epsilon*std::sqrt(rx*rx+ry*ry)) 
        * sin(bx*Deltax+by*Deltay+(0.5-z)*(rx*Deltax+ry*Deltay)) * (1.0-exp(-0.25*sigma_0*Q_s_0*Q_s_0/(2.0*M_PI*B_G)*exp(-(bx*bx+by*by)/(2.0*B_G))));
    }

    double Trans_imag (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        if (rx==0&&ry==0){rx=smallValue;ry=smallValue;}
        return 1.0/(2.0*M_PI) * A_Q * std::sqrt(2.0*m_Q_c*N_c) * e_Q * e * gsl_sf_bessel_K0(epsilon*std::sqrt(rx*rx+ry*ry)) 
        * cos(bx*Deltax+by*Deltay+(0.5-z)*(rx*Deltax+ry*Deltay)) * (1.0-exp(-0.25*sigma_0*Q_s_0*Q_s_0/(2.0*M_PI*B_G)*exp(-(bx*bx+by*by)/(2.0*B_G))));
    }

    double Longi_real (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        if (rx==0&&ry==0){rx=smallValue;ry=smallValue;}
        return 1.0/(1.0*M_PI) * A_Q * std::sqrt(2.0/m_Q_c*N_c) * e_Q * e * Q * z * (z-1.0) * gsl_sf_bessel_K0(epsilon*std::sqrt(rx*rx+ry*ry)) 
        * sin(bx*Deltax+by*Deltay+(0.5-z)*(rx*Deltax+ry*Deltay)) * (1.0-exp(-0.25*sigma_0*Q_s_0*Q_s_0/(2.0*M_PI*B_G)*exp(-(bx*bx+by*by)/(2.0*B_G))));
    }

    double Longi_imag (double Q, double bx, double by, double rx, double ry, double Deltax, double Deltay, double z)
    {
        if (rx==0&&ry==0){rx=smallValue;ry=smallValue;}
        return 1.0/(1.0*M_PI) * A_Q * std::sqrt(2.0/m_Q_c*N_c) * e_Q * e * Q * z * (z-1.0) * gsl_sf_bessel_K0(epsilon*std::sqrt(rx*rx+ry*ry)) 
        * cos(bx*Deltax+by*Deltay+(0.5-z)*(rx*Deltax+ry*Deltay)) * (1.0-exp(-0.25*sigma_0*Q_s_0*Q_s_0/(2.0*M_PI*B_G)*exp(-(bx*bx+by*by)/(2.0*B_G))));
    }
}