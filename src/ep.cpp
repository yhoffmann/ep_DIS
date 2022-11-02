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

struct USERDATA{ // not needed yet but might be needed later
    double Q;
    double Deltax;
    double Deltay;
};


using namespace std::complex_literals;

#include "Integrals.hpp"



int main(int argc, char** argv) {

    USERDATA data;

    std::vector<double> Return(2);
    
    Return = dCoherent_cross_section_dt(data);

    std::cout << "========\n" << "Coherent Cross Sections by t: \n" << "T: " << Return[0] << "  L: " << Return[1] << " [fm^2]" << std::endl;

    Return = dIncoherent_cross_section_dt(data);

    std::cout << "========\n" << "Incoherent Cross Sections by t: \n" << "T: " << Return[0] << "  L: " << Return[1] << " [fm^2]" << std::endl;

    return 0;
}