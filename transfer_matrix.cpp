#include "transfer_matrix.hpp"

#include <complex>
#include <iostream>
#include <math.h>

using namespace std;

TransferMatrix::TransferMatrix(double n0, double n2)
{
    this->n0 = n0;
    this->n2 = n2;
}

TransferMatrix::TransferMatrix(double n0, double n2, complex<double> q_in)
{
    this->n0 = n0;
    this->n2 = n2;
    this->q_input = q_in;
}


complex<double> TransferMatrix::get_q_at(double z, complex<double> q)
{
    double gamma = sqrt(this->n2 / this->n0);    // gamma parameter (sqrt(n2/n0))
    // cout<<gamma*z<<'\n';

    complex<double> A, B, C, D;   // elements of ABCD matrix
    A = cos(gamma * z);
    B = sin(gamma * z) / (this->n0 * gamma);
    C = -(this->n0 * gamma) * sin(gamma * z);
    D = cos(gamma * z);
    // cout<<A<<' '<<B<<' '<<C<<' '<<D<<'\n';
    complex<double> q_out = ((A * q) + (B * this->n0)) / ((C * q / this->n0) + D);
    return q_out;
}


complex<double> TransferMatrix::propagate(double z, complex<double> q)
{
    this->q_output = this->get_q_at(z, q);
    return this->q_output;
}

complex<double> TransferMatrix::propagate(double z)
{
    this->q_output = this->get_q_at(z, this->q_input);
    return this->q_output;
}

double TransferMatrix::calc_beam_radius_at_zero(double wavelength, std::complex<double> q)
{
    double waist = wavelength / M_PI * (q.imag() + pow(q.real(), 2)/q.imag());

    // double waist = sqrt(q.imag() * wavelength / M_PI);
    cout<<waist;
    return waist;
}

double TransferMatrix::beam_waist(double z, double wavelength, std::complex<double> q)
{
    return calc_beam_radius_at_zero(wavelength, q);
}

double TransferMatrix::beam_waist(double z, double wavelength)
{
    return calc_beam_radius_at_zero(wavelength, q_output);
}