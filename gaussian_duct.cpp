#include "gaussian_duct.hpp"
#include <math.h>
#include <iostream>

using namespace std;

GaussianDuct::GaussianDuct(double n_0, double n_2, double length) : TransferMatrix::TransferMatrix()
{
    n0 = n_0;
    n2 = n_2;
    z = length;

    initialize_matrix();
}

void GaussianDuct::initialize_matrix()
{
    double gamma = sqrt(n2 / n0);    // gamma parameter (sqrt(n2/n0))

    A = cos(gamma * z);
    B = sin(gamma * z) / (n0 * gamma);
    C = -(n0 * gamma) * sin(gamma * z);
    D = cos(gamma * z);
}

std::complex<double> GaussianDuct::propagate(std::complex<double> qin)
{
    return calc_q_out(n0, qin);
}

void GaussianDuct::set_length(double length)
{
    z = length;
    initialize_matrix();
}

void GaussianDuct::set_n0(double n_0)
{
    n0 = n_0;
    initialize_matrix();
}

void GaussianDuct::set_n2(double n_2)
{
    n2 = n_2;
    initialize_matrix();
}