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

// STRUCT FOR PASSING PARAMETERS TO INTEGRAND FUNCTION // might need to add functionality for passing epsilon to PsiPsi functions
struct USERDATA
{
    double Q;
    double Deltax;
    double Deltay;
    double z;
};


using namespace std::complex_literals;

#include "Integrals.cpp"


int main(int argc, char** argv)
{
    USERDATA data;

    data.Deltax = 1.0;
    data.Deltay = 0.0;
    data.Q = 1.0;
    data.z = 0.5;

    SetQuarkFlavor('c');

    std::vector<double> Return(2);
    
    Return = dCoherent_cross_section_dt(data);

    std::cout << "========\n" << "Coherent Cross Sections by t: \n" << "T: " << Return[0] << "  L: " << Return[1] << " [fm^2]" << std::endl;

    Return = dIncoherent_cross_section_dt(data);

    std::cout << "========\n" << "Incoherent Cross Sections by t: \n" << "T: " << Return[0] << "  L: " << Return[1] << " [fm^2]" << std::endl;

    // Test Output for multiple combinations of Delta and Q
    std::fstream OutStreamCoherent;
    std::fstream OutStreamIncoherent;
    OutStreamCoherent.open("Data/DeltaAndQ_Coherent_Test.txt");
    OutStreamIncoherent.open("Data/DeltaAndQ_Incoherent_Test.txt");

    std::vector<double> DeltaRange {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    std::vector<double> QRange {0.01, 0.02, 0.03, 0.05, 0.10, 0.15, 0.20, 0.30};

    for (int i=0; i<DeltaRange.size(); i++)
    {
        for (int j=0; j<QRange.size(); j++)
        {
            data.Deltax = DeltaRange[i]; // varying Parameters to pass to integration function
            data.Deltay = DeltaRange[i];
            data.Q = QRange[j];

            Return = dCoherent_cross_section_dt(data);
            OutStreamCoherent << data.Q << " " << data.Deltax << " " << Return[0] << " " << Return[1] << std::endl;

            //Return = dIncoherent_cross_section_dt(data);
            //OutStreamIncoherent << Return[0] << " ";
        }
        OutStreamCoherent << std::endl;
        OutStreamIncoherent << std::endl;
    }
    OutStreamCoherent.close();
    OutStreamIncoherent.close();

    return 0;
}