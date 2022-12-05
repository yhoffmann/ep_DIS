#define MINEVAL 100000
#define MAXEVAL 500000

#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

#define NNEW 1000
#define NMIN 2
#define FLATNESS 25.

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
    
    // INTEGRAL RANGES //
    double RangeFactor_b = 20.0;
    double RangeFactor_r = 20.0;

    double bxmin = -RangeFactor_b*R;
    double bxmax = RangeFactor_b*R;

    double bymin = -RangeFactor_b*R;
    double bymax = RangeFactor_b*R;

    double rxmin = -RangeFactor_r*R; // TODO Check r range (maybe depending on model)
    double rxmax = RangeFactor_r*R;

    double rymin = -RangeFactor_r*R;
    double rymax = RangeFactor_r*R;

    // DEFINING COORDINATE TRANSFORMS FOR CORRECT INTEGRATION RANGE //
    double bx = bxmin+(bxmax-bxmin)*xx[0];
    double by = bymin+(bymax-bymin)*xx[1];

    double rx = rxmin+(rxmax-rxmin)*xx[2];
    double ry = rymin+(rymax-rymin)*xx[3];
    
    // DEFINING JACOBIAN //
    double Jacobian = (bxmax-bxmin)*(bymax-bymin)*(rxmax-rxmin)*(rymax-rymin);

    // SEPERATING INTEGRANDS INTO REAL AND IMAG PART //
    ff[0] = Jacobian / (4.0*M_PI) * A_coherent_Integrand_Function::Trans_real(Q,bx,by,rx,ry,Deltax,Deltay,z);
    ff[1] = Jacobian / (4.0*M_PI) * A_coherent_Integrand_Function::Trans_imag(Q,bx,by,rx,ry,Deltax,Deltay,z);

    ff[2] = Jacobian / (4.0*M_PI) * A_coherent_Integrand_Function::Longi_real(Q,bx,by,rx,ry,Deltax,Deltay,z);
    ff[3] = Jacobian / (4.0*M_PI) * A_coherent_Integrand_Function::Longi_imag(Q,bx,by,rx,ry,Deltax,Deltay,z);

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
    
    // INTEGRAL RANGES //
    double RangeFactor_b = 20.0;
    double RangeFactor_r = 20.0;

    double bxmin = -RangeFactor_b*R;
    double bxmax = RangeFactor_b*R;

    double bymin = -RangeFactor_b*R;
    double bymax = RangeFactor_b*R;

    double rxmin = -RangeFactor_r*R; // TODO Check r range (maybe depending on model)
    double rxmax = RangeFactor_r*R;

    double rymin = -RangeFactor_r*R;
    double rymax = RangeFactor_r*R;

    // DEFINING COORDINATE TRANSFORMS FOR CORRECT INTEGRATION RANGE //
    double bx = bxmin+(bxmax-bxmin)*xx[0];
    double by = bymin+(bymax-bymin)*xx[1];

    double rx = rxmin+(rxmax-rxmin)*xx[2];
    double ry = rymin+(rymax-rymin)*xx[3];
    
    // DEFINING JACOBIAN //
    double Jacobian = (bxmax-bxmin)*(bymax-bymin)*(rxmax-rxmin)*(rymax-rymin);

    // SEPERATING INTEGRANDS INTO REAL AND IMAG PART //
    ff[0] = Jacobian / (4.0*M_PI) * A_coherent_integrand_Function_first_order::Trans_real(Q,bx,by,rx,ry,Deltax,Deltay,z);
    ff[1] = Jacobian / (4.0*M_PI) * A_coherent_integrand_Function_first_order::Trans_imag(Q,bx,by,rx,ry,Deltax,Deltay,z);

    ff[2] = Jacobian / (4.0*M_PI) * A_coherent_integrand_Function_first_order::Longi_real(Q,bx,by,rx,ry,Deltax,Deltay,z);
    ff[3] = Jacobian / (4.0*M_PI) * A_coherent_integrand_Function_first_order::Longi_imag(Q,bx,by,rx,ry,Deltax,Deltay,z);

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

    // INTEGRAL RANGES //
    double RangeFactor_b = 15.0;
    double RangeFactor_r = 15.0;

    double bxmin = -RangeFactor_b*R;
    double bxmax = RangeFactor_b*R;
    double bymin = -RangeFactor_b*R;
    double bymax = RangeFactor_b*R;

    double rxmin = -RangeFactor_r*R; // TODO Check r range (maybe depending on model)
    double rxmax = RangeFactor_r*R;
    double rymin = -RangeFactor_r*R;
    double rymax = RangeFactor_r*R;

    double bbarxmin = -RangeFactor_b*R;
    double bbarxmax = RangeFactor_b*R;
    double bbarymin = -RangeFactor_b*R;
    double bbarymax = RangeFactor_b*R;

    double rbarxmin = -RangeFactor_r*R; // TODO Check r range (maybe depending on model)
    double rbarxmax = RangeFactor_r*R;
    double rbarymin = -RangeFactor_r*R;
    double rbarymax = RangeFactor_r*R;

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
    double Jacobian = (bxmax-bxmin)*(bymax-bymin)*(rxmax-rxmin)*(rymax-rymin)*(bbarxmax-bbarxmin)*(bbarymax-bbarymin)*(rbarxmax-rbarxmin)*(rbarymax-rbarymin);

    // SEPERATING INTEGRANDS INTO REAL AND IMAG PART //
    ff[0] = Jacobian / (4.0*M_PI) / (4.0*M_PI) * A_incoherent_Integrand_Function::Trans_real(Q,bx,by,bbarx,bbary,rx,ry,rbarx,rbary,Deltax,Deltay,z,zbar);
    ff[1] = Jacobian / (4.0*M_PI) / (4.0*M_PI) * A_incoherent_Integrand_Function::Trans_imag(Q,bx,by,bbarx,bbary,rx,ry,rbarx,rbary,Deltax,Deltay,z,zbar);

    ff[2] = Jacobian / (4.0*M_PI) / (4.0*M_PI) * A_incoherent_Integrand_Function::Longi_real(Q,bx,by,bbarx,bbary,rx,ry,rbarx,rbary,Deltax,Deltay,z,zbar);
    ff[3] = Jacobian / (4.0*M_PI) / (4.0*M_PI) * A_incoherent_Integrand_Function::Longi_imag(Q,bx,by,bbarx,bbary,rx,ry,rbarx,rbary,Deltax,Deltay,z,zbar);

    return 0;
}



std::vector<double> dsigma_dt_coherent (USERDATA &parameters)
{
    std::vector<double> Return(2);

    // PROCEDURE MONITOR //
    int NumberOfRegions,NumberOfEvaluations,ErrorStatus;
        
    // VALUES AND ERRORS //
    int Num_of_Dimensions = 4;
    int Num_of_Integrals = 4;
    int Num_of_Points = 1;
    double epsrel = 1e-3;
    double epsabs = 1e-12;
    int flags1 = 0;
    int flags2 = 4;
    int seed = 0;

    cubareal Value[Num_of_Integrals], Error[Num_of_Integrals], Probability[Num_of_Integrals];

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

    // Return0 returns Cross section for transverse and Return1 longitudinal and conversion to nb
    Return[0] = 1.0/16.0/M_PI*(Value[0]*Value[0]+Value[1]*Value[1]) * 1.0e7;
    Return[1] = 1.0/16.0/M_PI*(Value[2]*Value[2]+Value[3]*Value[3]) * 1.0e7;

    return Return;
}

std::vector<double> dsigma_dt_coherent_first_oder (USERDATA &parameters)
{
    std::vector<double> Return(2);

    // PROCEDURE MONITOR //
    int NumberOfRegions,NumberOfEvaluations,ErrorStatus;
        
    // VALUES AND ERRORS //
    int Num_of_Dimensions = 4;
    int Num_of_Integrals = 4;
    int Num_of_Points = 1;
    double epsrel = 1e-3;
    double epsabs = 1e-12;
    int flags1 = 0;
    int flags2 = 4;
    int seed = 0;

    cubareal Value[Num_of_Integrals], Error[Num_of_Integrals], Probability[Num_of_Integrals];

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

    // Return0 returns Cross section for transverse and Return1 longitudinal and conversion to nb
    Return[0] = 1.0/16.0/M_PI*(Value[0]*Value[0]+Value[1]*Value[1]) * 1.0e7;
    Return[1] = 1.0/16.0/M_PI*(Value[2]*Value[2]+Value[3]*Value[3]) * 1.0e7;

    return Return;
}


std::vector<double> dsigma_dt_incoherent (USERDATA &parameters)
{
    std::vector<double> Return(2);

    // PROCEDURE MONITOR //
    int NumberOfRegions,NumberOfEvaluations,ErrorStatus;
        
    // VALUES AND ERRORS //
    int Num_of_Dimensions = 8;
    int Num_of_Integrals = 4;
    int Num_of_Points = 1;
    double epsrel = 1e-3;
    double epsabs = 1e-12;
    int flags1 = 0;
    int flags2 = 4;
    int seed = 0;

    cubareal Value[Num_of_Integrals],Error[Num_of_Integrals],Probability[Num_of_Integrals];

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

    // Return0 returns Cross section for transverse and Return1 longitudinal and conversion to nb
    Return[0] = 1.0/16.0/M_PI*std::sqrt(Value[0]*Value[0]+Value[1]*Value[1]) * 1.0e7;
    Return[1] = 1.0/16.0/M_PI*std::sqrt(Value[2]*Value[2]+Value[3]*Value[3]) * 1.0e7;
    
    // Calculating Coherent Cross section for calculation of variance
    std::vector<double> CoherentReturn(2);
    CoherentReturn = dsigma_dt_coherent(parameters);

    // Calculating incoherent cross section T,L
    Return[0] = Return[0] - CoherentReturn[0];
    Return[1] = Return[1] - CoherentReturn[1];

    // Printing out terms of incoherent cross section, to test
    std::cout << "========\nFirst term is:\tT: " << Value[0] << "+i" << Value[1] << std::endl;
    std::cout << "\t\tL: " << Value[2] << "+i" << Value[3] << std::endl;

    std::cout << "========\nSecond term is:\tT: " << CoherentReturn[0] << std::endl;
    std::cout << "\t\tL: " << CoherentReturn[1] << std::endl;

    return Return;
}

/*
struct IntegrationParameters
{
    double Q;
    double Deltax;
    double Deltay;
    double z;
    double bx;
    double by;
};

static int rIntegrand_coherent (const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata)
{
    // GATHERING PARAMETERS FROM STRUCT
    IntegrationParameters parameters = *(IntegrationParameters *) userdata;
    double Q = parameters.Q;
    double Deltax = parameters.Deltax;
    double Deltay = parameters.Deltay;
    double z = parameters.z;
    double bx = parameters.bx;
    double by = parameters.by;
    
    // Integration ranges
    double RangeFactor = 20.0;
    double rxmax = RangeFactor*R;
    double rymax = RangeFactor*R;
    double rxmin = -RangeFactor*R;
    double rymin = -RangeFactor*R;

    double rx = rxmin+(rxmax-rxmin)*xx[0];
    double ry = rymin+(rymax-rymin)*xx[1];

    double Jacobian = (rxmax-rxmin)*(rymax-rymin);

    ff[0] = Jacobian / (4.0*M_PI) * A_coherent_integrand_Function_first_order::Trans_real(Q,bx,by,rx,ry,Deltax,Deltay,z);
    ff[1] = Jacobian / (4.0*M_PI) * A_coherent_integrand_Function_first_order::Trans_imag(Q,bx,by,rx,ry,Deltax,Deltay,z);
    ff[2] = Jacobian / (4.0*M_PI) * A_coherent_integrand_Function_first_order::Longi_real(Q,bx,by,rx,ry,Deltax,Deltay,z);
    ff[3] = Jacobian / (4.0*M_PI) * A_coherent_integrand_Function_first_order::Longi_imag(Q,bx,by,rx,ry,Deltax,Deltay,z);

    return 0;
}


std::vector<double> Make_r_Integration (IntegrationParameters &parameters)
{
    std::vector<double> Return(2);

    // PROCEDURE MONITOR //
    int NumberOfRegions,NumberOfEvaluations,ErrorStatus;
        
    // VALUES AND ERRORS //
    int Num_of_Dimensions = 2;
    int Num_of_Integrals = 4;
    int Num_of_Points = 1;
    double epsrel = 1e-3;
    double epsabs = 1e-12;
    int flags1 = 0;
    int flags2 = 4;
    int seed = 0;

    cubareal Value[Num_of_Integrals],Error[Num_of_Integrals],Probability[Num_of_Integrals];

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

    
}
*/

std::vector<double> Integration_bitbybit (USERDATA &parameters) 

{
    std::vector<double> Return(2);

    // PROCEDURE MONITOR //
    int NumberOfRegions,NumberOfEvaluations,ErrorStatus;
        
    // VALUES AND ERRORS //
    int Num_of_Dimensions = 8;
    int Num_of_Integrals = 4;
    int Num_of_Points = 1;
    double epsrel = 1e-3;
    double epsabs = 1e-12;
    int flags1 = 0;
    int flags2 = 4;
    int seed = 0;

    cubareal Value[Num_of_Integrals],Error[Num_of_Integrals],Probability[Num_of_Integrals];

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

    // Return0 returns Cross section for transverse and Return1 longitudinal and conversion to nb
    Return[0] = 1.0/16.0/M_PI*std::sqrt(Value[0]*Value[0]+Value[1]*Value[1]) * 1.0e7;
    Return[1] = 1.0/16.0/M_PI*std::sqrt(Value[2]*Value[2]+Value[3]*Value[3]) * 1.0e7;
    
    // Calculating Coherent Cross section for calculation of variance
    std::vector<double> CoherentReturn(2);
    CoherentReturn = dsigma_dt_coherent(parameters);

    // Calculating incoherent cross section T,L
    Return[0] = Return[0] - CoherentReturn[0];
    Return[1] = Return[1] - CoherentReturn[1];

    // Printing out terms of incoherent cross section, to test
    std::cout << "========\nFirst term is:\tT: " << Value[0] << "+i" << Value[1] << std::endl;
    std::cout << "\t\tL: " << Value[2] << "+i" << Value[3] << std::endl;

    std::cout << "========\nSecond term is:\tT: " << CoherentReturn[0] << std::endl;
    std::cout << "\t\tL: " << CoherentReturn[1] << std::endl;

    return Return;
}


static int A_full (const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata)
{
    // GATHERING PARAMETERS FROM STRUCT //
    USERDATA parameters = *(USERDATA *) userdata;

    double Deltax = parameters.Deltax;
    double Deltay = parameters.Deltay;
    double Q = parameters.Q;
    double z = parameters.z;
    
    // INTEGRAL RANGES //
    double RangeFactor_b = 20.0;
    double RangeFactor_r = 20.0;

    double bxmin = -RangeFactor_b*R;
    double bxmax = RangeFactor_b*R;

    double bymin = -RangeFactor_b*R;
    double bymax = RangeFactor_b*R;

    double rxmin = -RangeFactor_r*R; // TODO Check r range (maybe depending on model)
    double rxmax = RangeFactor_r*R;

    double rymin = -RangeFactor_r*R;
    double rymax = RangeFactor_r*R;

    // DEFINING COORDINATE TRANSFORMS FOR CORRECT INTEGRATION RANGE //
    double bx = bxmin+(bxmax-bxmin)*xx[0];
    double by = bymin+(bymax-bymin)*xx[1];

    double rx = rxmin+(rxmax-rxmin)*xx[2];
    double ry = rymin+(rymax-rymin)*xx[3];
    
    // DEFINING JACOBIAN //
    double Jacobian = (bxmax-bxmin)*(bymax-bymin)*(rxmax-rxmin)*(rymax-rymin);

    // SEPERATING INTEGRANDS INTO REAL AND IMAG PART //
    ff[0] = Jacobian * Integrand_full::Trans_real(Q,bx,by,rx,ry,Deltax,Deltay,z);
    ff[1] = Jacobian * Integrand_full::Trans_imag(Q,bx,by,rx,ry,Deltax,Deltay,z);

    ff[2] = Jacobian * Integrand_full::Longi_real(Q,bx,by,rx,ry,Deltax,Deltay,z);
    ff[3] = Jacobian * Integrand_full::Longi_imag(Q,bx,by,rx,ry,Deltax,Deltay,z);

    return 0;
}

std::vector<double> dsigma_dt_full (USERDATA &parameters)
{
    std::vector<double> Return(2);

    // PROCEDURE MONITOR //
    int NumberOfRegions,NumberOfEvaluations,ErrorStatus;
        
    // VALUES AND ERRORS //
    int Num_of_Dimensions = 4;
    int Num_of_Integrals = 4;
    int Num_of_Points = 1;
    double epsrel = 1e-3;
    double epsabs = 1e-12;
    int flags1 = 0;
    int flags2 = 4;
    int seed = 1337;

    cubareal Value[Num_of_Integrals], Error[Num_of_Integrals], Probability[Num_of_Integrals];

    Suave(Num_of_Dimensions, Num_of_Integrals,
        A_full, &parameters, Num_of_Points,
        epsrel, epsabs,
        flags1 | flags2, seed,
        MINEVAL, MAXEVAL,
        NNEW, NMIN,
        FLATNESS, STATEFILE, SPIN,
        &NumberOfRegions,&NumberOfEvaluations, &ErrorStatus, 
        Value, Error, Probability
    );

    // Return0 returns Cross section for transverse and Return1 longitudinal and conversion to nb
    Return[0] = 1.0/(16.0*M_PI)*(Value[0]*Value[0]+Value[1]*Value[1]) * 1.0e7;
    Return[1] = 1.0/(16.0*M_PI)*(Value[2]*Value[2]+Value[3]*Value[3]) * 1.0e7;

    return Return;
}