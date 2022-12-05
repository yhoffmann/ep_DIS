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


void Check_Commandline_Arguments(int NumberOfArguments)
{
    if (NumberOfArguments != 5)
    {
        std::cerr << "Error. Wrong number of command line arguments" << std::endl;
        exit(0);
    }
}

using namespace std::complex_literals;

#include "Integrals.cpp"


int main(int argc, char** argv)
{
    Check_Commandline_Arguments(argc);

    USERDATA data;

    data.Q = std::atof(argv[1]);
    data.Deltax = std::atof(argv[2]);
    data.Deltay = std::atof(argv[3]);
    data.z = std::atof(argv[4]);

    SetQuarkFlavor('c');

    std::vector<double> Return(2);
    
    Return = dsigma_dt_coherent(data);

    std::cout << "Coherent Cross Sections by t:\n" << "T: " << Return[0] << "  L: " << Return[1] << " [nb]\n======" << std::endl;

    Return = dsigma_dt_incoherent(data);

    std::cout << "Incoherent Cross Sections by t:\n" << "T: " << Return[0] << "  L: " << Return[1] << " [nb]\n======" << std::endl;

    double dsigma_dt_Analytical_Trans = dsigma_dt_coherent_analytical::Trans(data.Q,data.Deltax,data.Deltay,data.z);

    double dsigma_dt_Analytical_Longi = dsigma_dt_coherent_analytical::Longi(data.Q,data.Deltax,data.Deltay,data.z);

    std::cout << "###Analytical Result:\nCoherent Cross Sections by t:\n" << "T: " << dsigma_dt_Analytical_Trans << "  L: " << dsigma_dt_Analytical_Longi << " [nb]\n======" << std::endl;

    Return = dsigma_dt_coherent_first_oder(data);

    std::cout << "###Numerical first order:\nCoherent Cross Sections by t:\n" << "T: " << Return[0] << "  L: " << Return[1] << " [nb]\n======" << std::endl;

    // Test Output for multiple combinations of Delta and Q
    if (true)
    {
        // For version with .real() and .imag()
        std::fstream OutStreamCoherent;
        std::fstream OutStreamIncoherent;
        OutStreamCoherent.open("Data/DeltaAndQ_Coherent_Test.txt");
        OutStreamIncoherent.open("Data/DeltaAndQ_Incoherent_Test.txt");
        OutStreamCoherent << "#Q, Delta squared, T, L, analytical T, analytical L" << std::endl;
        OutStreamIncoherent << "#Q, Delta squared, T, L" << std::endl;

        std::vector<double> DeltaRange {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};

        std::vector<double> QRange {0.01, 0.02, 0.03, 0.05, 0.10, 0.15, 0.20, 0.30};

        std::vector<double> Return2(2);

        for (long unsigned int i=0; i<DeltaRange.size(); i++)
        {
            for (long unsigned int j=0; j<QRange.size(); j++)
            {
                data.Deltax = DeltaRange[i]; // varying Parameters to pass to integration function
                data.Deltay = DeltaRange[i];
                data.Q = QRange[j];

                Return = dsigma_dt_coherent(data);
                dsigma_dt_Analytical_Trans = dsigma_dt_coherent_analytical::Trans(data.Q,data.Deltax,data.Deltay,data.z);
                dsigma_dt_Analytical_Longi = dsigma_dt_coherent_analytical::Longi(data.Q,data.Deltax,data.Deltay,data.z);
                Return2 = dsigma_dt_full(data);
                OutStreamCoherent << data.Q << " " << data.Deltax*data.Deltax << " " << Return[0] << " " << Return[1] << " " << dsigma_dt_Analytical_Trans << " " << dsigma_dt_Analytical_Longi << " " << Return[0] << " " << Return[1] << std::endl;

                std::cout << Return[1]/Return2[1] << " ";

                //Return = dIncoherent_cross_section_dt(data);
                //OutStreamIncoherent << Return[0] << " ";
            }
            OutStreamCoherent << std::endl;
            OutStreamIncoherent << std::endl;
        }
        std::cout << std::endl;
        OutStreamCoherent.close();
        OutStreamIncoherent.close();
    }
    
    return 0;
}
