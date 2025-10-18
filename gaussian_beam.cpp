#include "gaussian_beam.hpp"
#include <math.h>

using namespace std;

double GaussianBeam::calc_beam_waist(double wavelength, std::complex<double> q)
{
    double w0 = wavelength / M_PI * q.imag();
    return w0;
}

double GaussianBeam::calc_beam_radius(double wavelength, std::complex<double> q)
{
    double w = wavelength / M_PI * (q.imag() + (q.real() * q.real())/q.imag());
    return w;
}