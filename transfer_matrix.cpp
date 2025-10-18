#include "transfer_matrix.hpp"

#include <complex>
#include <iostream>
#include <math.h>

using namespace std;


complex<double> TransferMatrix::calc_q_out(double n0, complex<double> q)
{
    complex<double> q_out = ((A * q) + (B * n0)) / ((C * q) + (D * n0));
    return q_out;
}