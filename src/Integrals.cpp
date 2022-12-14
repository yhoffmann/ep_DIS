#define MINEVAL 100000
#define MAXEVAL 1000000

#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

#define NNEW 1000
#define NMIN 2
#define FLATNESS 100.

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 0

#include "ModelsAndFunctions.cpp"

int UseSuave; // initialization only, value is read in from command line in function

// Struct for passing parameters and integration ranges to integrand function
struct USERDATA
{
    double Q;
    double Deltax;
    double Deltay;
    double z;
    int WhichIntegrand;
    double bmin;
    double bmax;
};


void FindRoot_Coherent(int n, USERDATA &parameters)
{
    if (parameters.Deltay != 0.0) {std::cerr << "Deltay needs to be 0 for FindRoot to work properly!" << std::endl;}

    switch (parameters.WhichIntegrand)
    {
        case 0: //Trans real, sin
            parameters.bmin = -double(n)*M_PI/(parameters.Deltax*GeVTofmm1);
            parameters.bmax = -double(n-1)*M_PI/(parameters.Deltax*GeVTofmm1);
            break;
        case 1: //Trans imag, cos
            parameters.bmin = (double(n)-0.5)*M_PI/(parameters.Deltax*GeVTofmm1);
            parameters.bmax = (double(n)+0.5)*M_PI/(parameters.Deltax*GeVTofmm1);
            break;
        case 2: //Longi real, sin
            parameters.bmin = -double(n)*M_PI/(parameters.Deltax*GeVTofmm1);
            parameters.bmax = -double(n-1)*M_PI/(parameters.Deltax*GeVTofmm1);
            break;
        case 3: //Longi imag, cos
            parameters.bmin = (double(n)-0.5)*M_PI/(parameters.Deltax*GeVTofmm1);
            parameters.bmax = (double(n)+0.5)*M_PI/(parameters.Deltax*GeVTofmm1);
            break;
        default:
            std::cerr << "Error in FindRoot, WhichFuntion out of range. Exiting!" << std::endl;
            break;
    }
}

// DEFINING MODEL // (done in functions directly)



// INTEGRAND FUNCTIONS //

static int A_coherent_Integrand (const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata)
{
    // GATHERING PARAMETERS FROM STRUCT //
    USERDATA parameters = *(USERDATA *) userdata;

    double Deltax = parameters.Deltax;
    double Deltay = parameters.Deltay;
    double Q = parameters.Q;
    double z = parameters.z;
    int WhichIntegrand = parameters.WhichIntegrand;
    
    // INTEGRAL RANGES //
    double b_RangeFactor = std::sqrt(10.0);
    double b_SignificantRange = GeVm1Tofm*std::sqrt(2.0*BG); // in fm

    double r_RangeFactor = 10.0;
    double r_SignificantRange = 1.0/(epsilonFunc(Q,z)*GeVTofmm1); // in fm

    double bxmin = -b_RangeFactor*b_SignificantRange; // in fm
    double bxmax = b_RangeFactor*b_SignificantRange;

    double bymin = -b_RangeFactor*b_SignificantRange;
    double bymax = b_RangeFactor*b_SignificantRange;

    double rxmin = -r_RangeFactor*r_SignificantRange; // in fm
    double rxmax = r_RangeFactor*r_SignificantRange;

    double rymin = -r_RangeFactor*r_SignificantRange;
    double rymax = r_RangeFactor*r_SignificantRange;

    // DEFINING COORDINATE TRANSFORMS FOR CORRECT INTEGRATION RANGE //
    double bx = bxmin+(bxmax-bxmin)*xx[0];
    double by = bymin+(bymax-bymin)*xx[1];

    double rx = rxmin+(rxmax-rxmin)*xx[2];
    double ry = rymin+(rymax-rymin)*xx[3];
    
    // DEFINING JACOBIAN //
    double Jacobian = (bxmax-bxmin)*(bymax-bymin)*(rxmax-rxmin)*(rymax-rymin)*std::pow(fmToGeVm1,4.0); //in GeVm4 

    // seperating integrand into real and imag parts and defining actual integrand function to be used for integration based on parameter from struct
    switch (WhichIntegrand)
    {
    case 0:
        ff[0] = Jacobian * A_coherent_Integrand_Function::Trans_real(Q,bx,by,rx,ry,Deltax,Deltay,z);
        break;
    case 1:
        ff[0] = Jacobian * A_coherent_Integrand_Function::Trans_imag(Q,bx,by,rx,ry,Deltax,Deltay,z);
        break;
    case 2:
        ff[0] = Jacobian * A_coherent_Integrand_Function::Longi_real(Q,bx,by,rx,ry,Deltax,Deltay,z);
        break;
    case 3:
        ff[0] = Jacobian * A_coherent_Integrand_Function::Longi_imag(Q,bx,by,rx,ry,Deltax,Deltay,z);
        break;
    default:
        std::cout << "No Integrand funciton for this parameter." << std::endl;
        // have to figure out a way to stop integration function
        break;
    }

    return 0;
}





static int A_coherent_Integrand_first_order (const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata)
{
    // GATHERING PARAMETERS FROM STRUCT //
    USERDATA parameters = *(USERDATA *) userdata;

    double Deltax = parameters.Deltax;
    double Deltay = parameters.Deltay;
    double Q = parameters.Q;
    double z = parameters.z;
    int WhichIntegrand = parameters.WhichIntegrand;
    
    // INTEGRAL RANGES //
    double b_RangeFactor = std::sqrt(10.0);
    double b_SignificantRange = GeVm1Tofm*std::sqrt(2.0*BG); // in fm

    double r_RangeFactor = 10.0;
    double r_SignificantRange = 1.0/(epsilonFunc(Q,z)*GeVTofmm1); // in fm

    double bxmin = -b_RangeFactor*b_SignificantRange; // in fm
    double bxmax = b_RangeFactor*b_SignificantRange;

    double bymin = -b_RangeFactor*b_SignificantRange;
    double bymax = b_RangeFactor*b_SignificantRange;

    double rxmin = -r_RangeFactor*r_SignificantRange; // in fm
    double rxmax = r_RangeFactor*r_SignificantRange;

    double rymin = -r_RangeFactor*r_SignificantRange;
    double rymax = r_RangeFactor*r_SignificantRange;

    // DEFINING COORDINATE TRANSFORMS FOR CORRECT INTEGRATION RANGE //
    double bx = bxmin+(bxmax-bxmin)*xx[0];
    double by = bymin+(bymax-bymin)*xx[1];

    double rx = rxmin+(rxmax-rxmin)*xx[2];
    double ry = rymin+(rymax-rymin)*xx[3];
    
    // DEFINING JACOBIAN //
    double Jacobian = (bxmax-bxmin)*(bymax-bymin)*(rxmax-rxmin)*(rymax-rymin)*std::pow(fmToGeVm1,4.0); //in GeVm4 

    // seperating integrand into real and imag parts and defining actual integrand function to be used for integration based on parameter from struct
    switch (WhichIntegrand)
    {
    case 0:
        ff[0] = Jacobian * A_coherent_integrand_Function_first_order::Trans_real(Q,bx,by,rx,ry,Deltax,Deltay,z);
        break;
    case 1:
        ff[0] = Jacobian * A_coherent_integrand_Function_first_order::Trans_imag(Q,bx,by,rx,ry,Deltax,Deltay,z);
        break;
    case 2:
        ff[0] = Jacobian * A_coherent_integrand_Function_first_order::Longi_real(Q,bx,by,rx,ry,Deltax,Deltay,z);
        break;
    case 3:
        ff[0] = Jacobian * A_coherent_integrand_Function_first_order::Longi_imag(Q,bx,by,rx,ry,Deltax,Deltay,z);
        break;
    default:
        std::cout << "No Integrand function for this parameter." << std::endl;
        // have to figure out a way to stop integration function
        break;
    }

    return 0;
}



static int A_incoherent_Integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata)
{
    // GATHERING PARAMETERS FROM STRUCT //
    USERDATA parameters = *(USERDATA *) userdata;

    double Deltax = parameters.Deltax;
    double Deltay = parameters.Deltay;
    double Q = parameters.Q;
    double z = parameters.z;
    double zbar = parameters.z;
    int WhichIntegrand = parameters.WhichIntegrand;

    // INTEGRAL RANGES //
    double b_RangeFactor = std::sqrt(10.0);
    double b_SignificantRange = GeVm1Tofm*std::sqrt(2.0*BG); // in fm

    double r_RangeFactor = 10.0;
    double r_SignificantRange = 1.0/(epsilonFunc(Q,z)*GeVTofmm1); // in fm

    double bxmin = -b_RangeFactor*b_SignificantRange; // in fm
    double bxmax = b_RangeFactor*b_SignificantRange;
    double bymin = -b_RangeFactor*b_SignificantRange;
    double bymax = b_RangeFactor*b_SignificantRange;

    double bbarxmin = -b_RangeFactor*b_SignificantRange; // in fm
    double bbarxmax = b_RangeFactor*b_SignificantRange;
    double bbarymin = -b_RangeFactor*b_SignificantRange;
    double bbarymax = b_RangeFactor*b_SignificantRange;

    double rxmin = -r_RangeFactor*r_SignificantRange; // in fm
    double rxmax = r_RangeFactor*r_SignificantRange;
    double rymin = -r_RangeFactor*r_SignificantRange;
    double rymax = r_RangeFactor*r_SignificantRange;

    double rbarxmin = -r_RangeFactor*r_SignificantRange; // in fm
    double rbarxmax = r_RangeFactor*r_SignificantRange;
    double rbarymin = -r_RangeFactor*r_SignificantRange;
    double rbarymax = r_RangeFactor*r_SignificantRange;

    // DEFINING COORDINATE TRANSFORMS FOR CORRECT INTEGRATION RANGE //
    double bx = bxmin+(bxmax-bxmin)*xx[0];
    double by = bymin+(bymax-bymin)*xx[1];

    double rx = rxmin+(rxmax-rxmin)*xx[2];
    double ry = rymin+(rymax-rymin)*xx[3];

    double bbarx = bbarxmin+(bbarxmax-bbarxmin)*xx[4];
    double bbary = bbarymin+(bbarymax-bbarymin)*xx[5];

    double rbarx = rbarxmin+(rbarxmax-rbarxmin)*xx[6];
    double rbary = rbarymin+(rbarymax-rbarymin)*xx[7];
    
    // DEFINING JACOBIAN //
    double Jacobian = (bxmax-bxmin)*(bymax-bymin)*(rxmax-rxmin)*(rymax-rymin)*(bbarxmax-bbarxmin)*(bbarymax-bbarymin)*(rbarxmax-rbarxmin)*(rbarymax-rbarymin)*std::pow(fmToGeVm1,8.0);

    // seperating integrand into real and imag parts and defining actual integrand function to be used for integration based on parameter from struct
    switch (WhichIntegrand)
    {
    case 0:
        ff[0] = Jacobian * A_incoherent_Integrand_Function::Trans_real(Q,bx,by,bbarx,bbary,rx,ry,rbarx,rbary,Deltax,Deltay,z,zbar);
        break;
    case 1:
        ff[0] = Jacobian * A_incoherent_Integrand_Function::Trans_imag(Q,bx,by,bbarx,bbary,rx,ry,rbarx,rbary,Deltax,Deltay,z,zbar);
        break;
    case 2:
        ff[0] = Jacobian * A_incoherent_Integrand_Function::Longi_real(Q,bx,by,bbarx,bbary,rx,ry,rbarx,rbary,Deltax,Deltay,z,zbar);
        break;
    case 3:
        ff[0] = Jacobian * A_incoherent_Integrand_Function::Longi_imag(Q,bx,by,bbarx,bbary,rx,ry,rbarx,rbary,Deltax,Deltay,z,zbar);
        break;
    default:
        std::cout << "No Integrand function for this parameter." << std::endl;
        // have to figure out a way to stop integration function
        break;
    }

    return 0;
}



std::vector<double> dsigmabydt_coherent (USERDATA &parameters)
{
    std::vector<double> Return(2);

    // PROCEDURE MONITOR //
    int NumberOfRegions,NumberOfEvaluations,ErrorStatus;
        
    // VALUES AND ERRORS //
    int Num_of_Dimensions = 4;
    int Num_of_Integrals = 1;
    int Num_of_Points = 1;
    double epsrel = 1e-3;
    double epsabs = 1e-12;
    int flags1 = 0;
    int flags2 = 4;
    int seed = time(0);

    cubareal Value[Num_of_Integrals], Error[Num_of_Integrals], Probability[Num_of_Integrals];

    double IntegrationResults[4];
    for (parameters.WhichIntegrand=0; parameters.WhichIntegrand<=3; parameters.WhichIntegrand++){
        if (UseSuave)
        {
            Suave(Num_of_Dimensions, Num_of_Integrals,
                A_coherent_Integrand, &parameters, Num_of_Points,
                epsrel, epsabs,
                flags1 | flags2, seed,
                MINEVAL, MAXEVAL,
                NNEW, NMIN,
                FLATNESS, STATEFILE, SPIN,
                &NumberOfRegions,&NumberOfEvaluations, &ErrorStatus, 
                Value, Error, Probability
            );
            //std::cout << "ErrorStatus dsigmabydt_coherent " << ErrorStatus << std::endl;
            // saving integration result of each run in list for later use
            IntegrationResults[parameters.WhichIntegrand] = Value[0];
        }
        else
        {
            Cuhre(Num_of_Dimensions, Num_of_Integrals,
                A_coherent_Integrand, &parameters, Num_of_Points,
                epsrel, epsabs,
                flags1 | flags2,
                MINEVAL, MAXEVAL,
                KEY, STATEFILE, SPIN,
                &NumberOfRegions, &NumberOfEvaluations, &ErrorStatus,
                Value, Error, Probability
            );
            IntegrationResults[parameters.WhichIntegrand] = Value[0];
        }
    }
    
    // Return0 returns Cross section for transverse and Return1 longitudinal and conversion to nb
    Return[0] = 1.0/16.0/M_PI*(IntegrationResults[0]*IntegrationResults[0]+IntegrationResults[1]*IntegrationResults[1]) * (GeVm1Tofm*GeVm1Tofm)*fm2TonB;
    Return[1] = 1.0/16.0/M_PI*(IntegrationResults[2]*IntegrationResults[2]+IntegrationResults[3]*IntegrationResults[3]) * (GeVm1Tofm*GeVm1Tofm)*fm2TonB;

    return Return;
}



std::vector<double> dsigmabydt_coherent_first_order (USERDATA &parameters)
{
    // PROCEDURE MONITOR //
    int NumberOfRegions,NumberOfEvaluations,ErrorStatus;
        
    // VALUES AND ERRORS //
    int Num_of_Dimensions = 4;
    int Num_of_Integrals = 1;
    int Num_of_Points = 1;
    double epsrel = 1e-3;
    double epsabs = 1e-12;
    int flags1 = 0;
    int flags2 = 4;
    int seed = time(0);

    cubareal Value[Num_of_Integrals], Error[Num_of_Integrals], Probability[Num_of_Integrals];

    double IntegrationResults[4];
    for (parameters.WhichIntegrand=0; parameters.WhichIntegrand<=3; parameters.WhichIntegrand++){
        if (UseSuave)
        {
            Suave(Num_of_Dimensions, Num_of_Integrals,
                A_coherent_Integrand_first_order, &parameters, Num_of_Points,
                epsrel, epsabs,
                flags1 | flags2, seed,
                MINEVAL, MAXEVAL,
                NNEW, NMIN,
                FLATNESS, STATEFILE, SPIN,
                &NumberOfRegions,&NumberOfEvaluations, &ErrorStatus, 
                Value, Error, Probability
            );
            //std::cout << "ErrorStatus dsigmabydt_coherent " << ErrorStatus << std::endl;
            // saving integration result of each run in list for later use
            IntegrationResults[parameters.WhichIntegrand] = Value[0];
        }
        else
        {
            Cuhre(Num_of_Dimensions, Num_of_Integrals,
                A_coherent_Integrand_first_order, &parameters, Num_of_Points,
                epsrel, epsabs,
                flags1 | flags2,
                MINEVAL, MAXEVAL,
                KEY, STATEFILE, SPIN,
                &NumberOfRegions, &NumberOfEvaluations, &ErrorStatus,
                Value, Error, Probability
            );
            IntegrationResults[parameters.WhichIntegrand] = Value[0];
        }
    }

    std::vector<double> Return(2);
    // Return0 returns Cross section for transverse and Return1 longitudinal and conversion to nb
    Return[0] = 1.0/16.0/M_PI*(IntegrationResults[0]*IntegrationResults[0]+IntegrationResults[1]*IntegrationResults[1]) * (GeVm1Tofm*GeVm1Tofm)*fm2TonB;
    Return[1] = 1.0/16.0/M_PI*(IntegrationResults[2]*IntegrationResults[2]+IntegrationResults[3]*IntegrationResults[3]) * (GeVm1Tofm*GeVm1Tofm)*fm2TonB;


    return Return;
}


namespace dsigmabydt_coherent_analytical_first_order
{
    double Trans (double Q, double Deltax, double Deltay, double z)
    {
        std::complex<double> Amplitude = -1.0i * (-A_Q) * std::sqrt(2.0*m_Q_c*N_c) * e_Q * e * (-sigma_0*Q_s_0*Q_s_0) * exp(-BG*(Deltax*Deltax+Deltay*Deltay)/2.0) / std::pow(epsilonFunc(Q,z),4.0); // GeVm2
        return std::pow(std::abs(Amplitude),2.0) / 16.0 / M_PI * (GeVm1Tofm*GeVm1Tofm)*fm2TonB; // GeVm2 nb
    }

    double Longi (double Q, double Deltax, double Deltay, double z)
    {
        std::complex<double> Amplitude = -1.0i * 0.5 * (-A_Q) * std::sqrt(2.0*N_c/m_Q_c) * e_Q * e * Q * (-sigma_0*Q_s_0*Q_s_0) * exp(-BG*(Deltax*Deltax+Deltay*Deltay)/2.0) / std::pow(epsilonFunc(Q,z),4.0); //GeVm2
        return std::pow(std::abs(Amplitude),2.0) / 16.0 / M_PI * (GeVm1Tofm*GeVm1Tofm)*fm2TonB; // GeVm2 nb
    }
}



std::vector<double> dsigmabydt_incoherent (USERDATA &parameters)
{
    std::vector<double> Return(2);

    // PROCEDURE MONITOR //
    int NumberOfRegions,NumberOfEvaluations,ErrorStatus;
        
    // VALUES AND ERRORS //
    int Num_of_Dimensions = 8;
    int Num_of_Integrals = 1;
    int Num_of_Points = 1;
    double epsrel = 1e-3;
    double epsabs = 1e-12;
    int flags1 = 0;
    int flags2 = 4;
    int seed = time(0);

    cubareal Value[Num_of_Integrals],Error[Num_of_Integrals],Probability[Num_of_Integrals];

    double IntegrationResults[4];
    for (parameters.WhichIntegrand=0; parameters.WhichIntegrand<=3; parameters.WhichIntegrand++){
        if (UseSuave)
        {
            Suave(Num_of_Dimensions, Num_of_Integrals,
                A_incoherent_Integrand, &parameters, Num_of_Points,
                epsrel, epsabs,
                flags1 | flags2, seed,
                MINEVAL, MAXEVAL,
                NNEW, NMIN,
                FLATNESS, STATEFILE, SPIN,
                &NumberOfRegions,&NumberOfEvaluations, &ErrorStatus, 
                Value, Error, Probability
            );
            //std::cout << "ErrorStatus dsigmabydt_incoherent " << ErrorStatus << std::endl;
            // saving integration result of each run in list for later use
            IntegrationResults[parameters.WhichIntegrand] = Value[0];
        }
        else
        {
            Cuhre(Num_of_Dimensions, Num_of_Integrals,
                A_incoherent_Integrand, &parameters, Num_of_Points,
                epsrel, epsabs,
                flags1 | flags2,
                MINEVAL, MAXEVAL,
                KEY, STATEFILE, SPIN,
                &NumberOfRegions, &NumberOfEvaluations, &ErrorStatus,
                Value, Error, Probability
            );
            IntegrationResults[parameters.WhichIntegrand] = Value[0];
        }
    }
    
    
    // Calculating Transverse ([1]) and Longitudinal ([2]) cross sections tand conversion to nb GeVm2
    Return[0] = 1.0/16.0/M_PI*std::sqrt(IntegrationResults[0]*IntegrationResults[0]+IntegrationResults[1]*IntegrationResults[1]) * (GeVm1Tofm*GeVm1Tofm)*fm2TonB;
    Return[1] = 1.0/16.0/M_PI*std::sqrt(IntegrationResults[2]*IntegrationResults[2]+IntegrationResults[3]*IntegrationResults[3]) * (GeVm1Tofm*GeVm1Tofm)*fm2TonB;
        
    // Calculating Coherent Cross section for calculation of variance
    std::vector<double> CoherentReturn(2);
    CoherentReturn = dsigmabydt_coherent(parameters);

/*
    // Printing out terms of incoherent cross section, to test
    std::cout << "========\nFirst term is:\tT: " << Return[0] << std::endl;
    std::cout << "\t\tL: " << Return[1] << std::endl;
    std::cout << "========\nSecond term is:\tT: " << CoherentReturn[0] << std::endl;
    std::cout << "\t\tL: " << CoherentReturn[1] << std::endl;
*/
    // Calculating incoherent cross section T,L
    Return[0] = Return[0] - CoherentReturn[0];
    Return[1] = Return[1] - CoherentReturn[1];

    return Return;
}