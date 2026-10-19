#include "propagation_material.hpp"
#include <complex>
#include <math.h>

using namespace std;

Material::Material(complex<double> n_0, complex<double> n_2, double lambda, double length, int material, double amp) : TransferMatrix::TransferMatrix()
{
    n0 = n_0;
    n2 = n_2;
    wavelength = lambda;
    wavenumber = 2 * M_PI * n0 / lambda;
    z = length;
    material_type = material;
    amp_0 = amp;
    this->amp = amp_0;

    initialize_matrix();
}

void Material::initialize_matrix()
{
    if(material_type == 0)       // free space
    {
        A = 1;
        B = z/n0;
        C = 0;
        D = 1;
    }
    else if(material_type == 1)  // gaussian duct
    {
        complex<double> gamma = sqrt(n2/n0);    // gamma parameter
        if(abs(gamma)!= 0)
        {
            // A = complex<double> (cos(gamma.real()*z) * cosh(gamma.imag()*z), -sin(gamma.real()*z) * sinh(gamma.imag()*z));
            // B = complex<double> (sin(gamma.real()*z) * cosh(gamma.imag()*z), cos(gamma.real()*z) * sinh(gamma.imag()*z)) / (n0 * gamma);
            // C = -(n0 * gamma) * complex<double> (sin(gamma.real()*z) * cosh(gamma.imag()*z), cos(gamma.real()*z) * sinh(gamma.imag()*z));
            // D = complex<double> (cos(gamma.real()*z) * cosh(gamma.imag()*z), -sin(gamma.real()*z) * sinh(gamma.imag()*z));
            A = cos(gamma * z);
            B = sin(gamma * z) / (n0 * gamma);
            C = -n0 * gamma * sin(gamma*z);
            D = cos(gamma*z);
        }
        else
        {
            A = 1;
            B = z/n0;
            C = 0;
            D = 1;
        }
    }
}

std::complex<double> Material::propagate(std::complex<double> qin)
{
    return calc_q_out(qin, n0);
}

void Material::set_length(double length)
{
    z = length;
    initialize_matrix();
}

void Material::set_n0(double n_0)
{
    n0 = n_0;
    initialize_matrix();
}

void Material::set_n2(double n_2)
{
    n2 = n_2;
    initialize_matrix();
}

void Material::set_material(int material_type)
{
    this->material_type = material_type;
    initialize_matrix();
}

double Material::calc_amplitude(std::complex<double> qin)
{
    amp = amp / (A + B * n0 / qin) * exp(-std::complex<double>(0,1) * wavenumber * z);
    return abs(amp);
}