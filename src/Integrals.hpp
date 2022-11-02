#define MINEVAL 0
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

#include "ModelsAndFunctions.hpp"

// DEFINING MODEL // (done in functions directly)



// INTEGRAND FUNCTIONS //
static int A_coherent_Integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) {

    // INTEGRAL RANGES //
    double RangeFactor = 20.0;

    double bxmin = -RangeFactor*R;
    double bxmax = RangeFactor*R;

    double bymin = -RangeFactor*R;
    double bymax = RangeFactor*R;

    double rxmin = -RangeFactor*R; // TODO Check r range (maybe depending on model)
    double rxmax = RangeFactor*R;

    double rymin = -RangeFactor*R;
    double rymax = RangeFactor*R;

    // DEFINING COORDINATE TRANSFORMS FOR CORRECT INTEGRATION RANGE //
    double bx = bxmin+(bxmax-bxmin)*xx[0];
    double by = bymin+(bymax-bymin)*xx[1];

    double rx = rxmin+(rxmax-rxmin)*xx[2];
    double ry = rymin+(rymax-rymin)*xx[3];
    
    // DEFINING JACOBIAN //
    double Jacobian = (bxmax-bxmin)*(bymax-bymin)*(rxmax-rxmin)*(rymax-rymin);

    // DEFINING INTEGRAND FUNCTION //
    std::complex<double> Integrand_T_complex = Jacobian * 1.0i * PsiPsiSimpleModel::PsiPsi_T_no_delta(Q,rx,ry,0.5) * exp(-1.0i*(bx*Deltax+by*Deltay)) * GBWModel::dsigma_qq_d2b(bx,by,rx,ry);
    std::complex<double> Integrand_L_complex = Jacobian * 1.0i * PsiPsiSimpleModel::PsiPsi_L_no_delta(Q,rx,ry,0.5) * exp(-1.0i*(bx*Deltax+by*Deltay)) * GBWModel::dsigma_qq_d2b(bx,by,rx,ry);

    // SEPERATING INTEGRANDS INTO REAL AND IMAG PART //
    ff[0] = Integrand_T_complex.real();
    ff[1] = Integrand_T_complex.imag();

    ff[2] = Integrand_L_complex.real();
    ff[3] = Integrand_L_complex.imag();

    return 0;
}



static int A_incoherent_Integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) {

    // INTEGRAL RANGES //
    double RangeFactor = 20.0;

    double bxmin = -RangeFactor*R;
    double bxmax = RangeFactor*R;

    double bymin = -RangeFactor*R;
    double bymax = RangeFactor*R;

    double rxmin = -RangeFactor*R; // TODO Check r range (maybe depending on model)
    double rxmax = RangeFactor*R;

    double rymin = -RangeFactor*R;
    double rymax = RangeFactor*R;

    double bbarxmin = -RangeFactor*R;
    double bbarxmax = RangeFactor*R;

    double bbarymin = -RangeFactor*R;
    double bbarymax = RangeFactor*R;

    double rbarxmin = -RangeFactor*R; // TODO Check r range (maybe depending on model)
    double rbarxmax = RangeFactor*R;

    double rbarymin = -RangeFactor*R;
    double rbarymax = RangeFactor*R;

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

    // DEFINING INTEGRAND FUNCTION //
    std::complex<double> Integrand_T_complex = Jacobian * 1.0i * PsiPsiSimpleModel::PsiPsi_T_no_delta(Q,rx,ry,0.5) * exp(-1.0i*(bx*Deltax+by*Deltay)) * GBWModel::dsigma_qq_d2b(bx,by,rx,ry) * (-1.0i) * PsiPsiSimpleModel::PsiPsi_T_no_delta(Q,rbarx,rbary,0.5) * exp(+1.0i*(bbarx*Deltax+bbary*Deltay)) * GBWModel::dsigma_qq_d2b(bbarx,bbary,rbarx,rbary);
    std::complex<double> Integrand_L_complex = Jacobian * 1.0i * PsiPsiSimpleModel::PsiPsi_L_no_delta(Q,rx,ry,0.5) * exp(-1.0i*(bx*Deltax+by*Deltay)) * GBWModel::dsigma_qq_d2b(bx,by,rx,ry) * (-1.0i) * PsiPsiSimpleModel::PsiPsi_L_no_delta(Q,rbarx,rbary,0.5) * exp(+1.0i*(bbarx*Deltax+bbary*Deltay)) * GBWModel::dsigma_qq_d2b(bbarx,bbary,rbarx,rbary);

    // SEPERATING INTEGRANDS INTO REAL AND IMAG PART //
    ff[0] = Integrand_T_complex.real();
    ff[1] = Integrand_T_complex.imag();

    ff[2] = Integrand_L_complex.real();
    ff[3] = Integrand_L_complex.imag();

    return 0;
}



std::vector<double> dCoherent_cross_section_dt (USERDATA &parameters){

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

    // Return0 returns Cross section for transverse and Return1 longitudinal
    Return[0] = 1.0/16.0/M_PI*(Value[0]*Value[0]+Value[1]*Value[1]);
    Return[1] = 1.0/16.0/M_PI*(Value[2]*Value[2]+Value[3]*Value[3]);

    return Return;
}


std::vector<double> dIncoherent_cross_section_dt (USERDATA &parameters){

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

    // Return0 returns Cross section for transverse and Return1 longitudinal
    Return[0] = 1.0/16.0/M_PI*(Value[0]*Value[0]+Value[1]*Value[1]);
    Return[1] = 1.0/16.0/M_PI*(Value[2]*Value[2]+Value[3]*Value[3]);

    
    // Calculating Coherent Cross section for calculation of variance
    std::vector<double> CoherentReturn(2);

    CoherentReturn = dCoherent_cross_section_dt(parameters);

    // Calculating incoherent cross section T,L
    Return[0] = Return[0] - CoherentReturn[0];
    Return[1] = Return[1] - CoherentReturn[1];

    return Return;
}