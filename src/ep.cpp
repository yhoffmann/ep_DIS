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
struct USERDATA{
    double Q;
    double Deltax;
    double Deltay;
    double epsilon;
    double e_Q;
};


using namespace std::complex_literals;

#include "Integrals.hpp"



int main(int argc, char** argv) {

    USERDATA data;

    data.Deltax = 1.0;
    data.Deltay = 0.0;
    data.Q = 1.0;

    std::vector<double> Return(2);
    
    Return = dCoherent_cross_section_dt(data);

    std::cout << "========\n" << "Coherent Cross Sections by t: \n" << "T: " << Return[0] << "  L: " << Return[1] << " [fm^2]" << std::endl;

    Return = dIncoherent_cross_section_dt(data);

    std::cout << "========\n" << "Incoherent Cross Sections by t: \n" << "T: " << Return[0] << "  L: " << Return[1] << " [fm^2]" << std::endl;

    return 0;
}