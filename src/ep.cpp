#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cuba.h>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <string>
#include <complex>
#include <iomanip>

/*
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
*/
#include <gsl/gsl_sf_bessel.h>

using namespace std::complex_literals;

#include "Integrals.cpp"


void Check_Commandline_Arguments(int NumberOfArguments)
{
    if (NumberOfArguments != 6)
    {
        std::cerr << "Error. Wrong number of command line arguments. Run './ep.exe [Q] [Deltax] [Deltay] [z] [Suave or Cuhre]'." << std::endl;
        exit(0);
    }
}

void SpecifyIntegrator(char* WhichIntegrator)
{
    switch(*WhichIntegrator)
    {
        case 'S':
            UseSuave = 1;
            std::cout << "(Using Suave for numerical integration)" << std::endl;
            break;
        case 'C':
            UseSuave = 0;
            std::cout << "(Using Cuhre for numerical integration)" << std::endl;
            break;
        default:
            std::cerr << "Last command line argument should specify the integrator routine (S for Suave or C for Cuhre)" << std::endl;
            exit(0);
            break;
    }
}


int main(int argc, char** argv)
{

    Check_Commandline_Arguments(argc);
    SpecifyIntegrator(argv[5]);
    
    USERDATA data;

    data.Q = std::atof(argv[1]);
    data.Deltax = std::atof(argv[2]);
    data.Deltay = std::atof(argv[3]);
    data.z = std::atof(argv[4]);
    data.WhichIntegrand = -1; // -1 as failsafe: will stop program if not changed in integration function

    SetQuarkFlavor('c');
/*
    // defining vectors for return of coherent and incoherent cross sections
    std::vector<double> dsigmabydt_coherent_Result(2);
    std::vector<double> dsigmabydt_incoherent_Result(2);
    std::vector<double> dsigmabydt_coherent_BetweenRoots_Result(2);
    std::vector<double> dsigmabydt_coherent_polar_Result(2);

    // calculating cross sections and saving results in vectors
    dsigmabydt_coherent_Result = dsigmabydt_coherent(data);
    dsigmabydt_incoherent_Result = dsigmabydt_incoherent(data);
    dsigmabydt_coherent_polar_Result = dsigmabydt_coherent_polar(data);

    // printing cross sections results for normal integrals
    std::cout << "###Cross section results\nCoherent Cross Section:\nT: " << dsigmabydt_coherent_Result[0] << " \tL: " << dsigmabydt_coherent_Result[1] << std::endl;
    std::cout << "Incoherent Cross Section:\nT: " << dsigmabydt_incoherent_Result[0] << " \tL: " << dsigmabydt_incoherent_Result[1] << std::endl;
    std::cout << "Coherent Cross Section from between roots\nT: " << dsigmabydt_coherent_BetweenRoots_Result[0] << " \tL: " << dsigmabydt_coherent_BetweenRoots_Result[1] << std::endl;
    std::cout << "Coherent Cross Section with polar coordinate integration\nT: " << dsigmabydt_coherent_polar_Result[0] << " \tL: " << dsigmabydt_coherent_polar_Result[1] << std::endl;
    std::cout << "Karthesian/Polar:\nT: " << dsigmabydt_coherent_Result[0]/dsigmabydt_coherent_polar_Result[0] << " \tL: " << dsigmabydt_coherent_Result[1]/dsigmabydt_coherent_polar_Result[1] << std::endl;
*/
    // Comparing first order analytical result with first order integrated result
    std::vector<double> dsigmabydt_coherent_analytical_first_order_result(2);
    std::vector<double> dsigmabydt_coherent_integrated_first_order_result(2);
    std::vector<double> dsigmabydt_coherent_integrated_first_order_result_UsingBetweenRootCalc(2);

    dsigmabydt_coherent_analytical_first_order_result[0] = dsigmabydt_coherent_analytical_first_order::Trans(data.Q,data.Deltax,data.Deltay,data.z);
    dsigmabydt_coherent_analytical_first_order_result[1] = dsigmabydt_coherent_analytical_first_order::Longi(data.Q,data.Deltax,data.Deltay,data.z);
    dsigmabydt_coherent_integrated_first_order_result = dsigmabydt_coherent_first_order(data);
    dsigmabydt_coherent_integrated_first_order_result_UsingBetweenRootCalc = dsigmabydt_coherent_first_order_UsingBetweenRootCalc(data);

    std::cout << "###Comparing first order analytical with numerical result and numerical reusult with root calculation:\nA:\tT: " << dsigmabydt_coherent_analytical_first_order_result[0] << "\tL: " << dsigmabydt_coherent_analytical_first_order_result[1] << std::endl;
    std::cout << "N:\tT: " << dsigmabydt_coherent_integrated_first_order_result[0] << "\tL: " << dsigmabydt_coherent_integrated_first_order_result[1] << std::endl;
    std::cout << "NRoots:\tT: " << dsigmabydt_coherent_integrated_first_order_result_UsingBetweenRootCalc[0] << "\tL: " << dsigmabydt_coherent_integrated_first_order_result_UsingBetweenRootCalc[1] << std::endl;
    std::cout << "A/N:\tT: " << dsigmabydt_coherent_analytical_first_order_result[0]/dsigmabydt_coherent_integrated_first_order_result[0] << "\tL: " << dsigmabydt_coherent_analytical_first_order_result[1]/dsigmabydt_coherent_integrated_first_order_result[1] << std::endl;
    std::cout << "N/NRs:\tT: " << dsigmabydt_coherent_integrated_first_order_result[0]/dsigmabydt_coherent_integrated_first_order_result_UsingBetweenRootCalc[0] << "\tL: " << dsigmabydt_coherent_integrated_first_order_result[1]/dsigmabydt_coherent_integrated_first_order_result_UsingBetweenRootCalc[1] << std::endl;
    std::cout << "A/NRs:\tT: " << dsigmabydt_coherent_analytical_first_order_result[0]/dsigmabydt_coherent_integrated_first_order_result_UsingBetweenRootCalc[0] << "\tL: " << dsigmabydt_coherent_analytical_first_order_result[1]/dsigmabydt_coherent_integrated_first_order_result_UsingBetweenRootCalc[1] << std::endl;


    // Comparing first order analytical and numerical results plot
    std::ofstream OutStream;
    std::string filename = "Data/FirstOrderResultsAnalyticalVsNumerical.txt";
    OutStream.open(filename);
    // Closing program if it cannot open the file
    if (!OutStream.is_open()) exit(0);

    std::vector<double> QRange = {0.01, 0.02, 0.03, 0.05, 0.1, 0.15, 0.20, 0.3, 0.33, 0.37, 0.4};
    std::vector<double> DeltaRange = {0.001, 0.01, 0.02, 0.05, 0.1, 0.3, 0.5, 0.7, 1.0, 1.2, 1.5, 1.7, 2.0, 2.3, 2.7, 3.0};

    OutStream << "#Q, Delta, Ana T,L, Num T,L;   " << QRange.size() << " values of Q (for Gnuplot)" << std::endl;
    OutStream << "#1, 2,         3,4,     5,6" << std::endl;
    for (double &QValue : QRange)
    {
        OutStream << "'" << "Q=" << QValue << "'" << std::endl;
        for (double &DeltaValue : DeltaRange)
        {
            // Saving Q and Delta values in struct
            data.Q = QValue;
            data.Deltax = DeltaValue;

            // Calculating Cross sections with data from struct
            dsigmabydt_coherent_analytical_first_order_result[0] = dsigmabydt_coherent_analytical_first_order::Trans(data.Q,data.Deltax,data.Deltay,data.z);
            dsigmabydt_coherent_analytical_first_order_result[1] = dsigmabydt_coherent_analytical_first_order::Longi(data.Q,data.Deltax,data.Deltay,data.z);
            dsigmabydt_coherent_integrated_first_order_result = dsigmabydt_coherent_first_order(data);

            // Outputting calculated cross sections and corresponding Q and Delta to file
            OutStream << QValue << " " << DeltaValue << " " << dsigmabydt_coherent_analytical_first_order_result[0] << " " << dsigmabydt_coherent_analytical_first_order_result[1] << " " << dsigmabydt_coherent_integrated_first_order_result[0] << " " << dsigmabydt_coherent_integrated_first_order_result[1] << std::endl;

            // Prodedure monitor
            std::cout << QValue << " " << DeltaValue << std::endl; 
        }
        // Two break lines after each data set for single value of Q
        OutStream << "\n" << std::endl;
    }

    // Cleaning up
    OutStream.close();

    return 0;
}
