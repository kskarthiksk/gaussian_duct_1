#include "free_space.hpp"
#include <math.h>
#include <iostream>

using namespace std;

FreeSpace::FreeSpace(double n_0, double length) : TransferMatrix::TransferMatrix()
{
    n0 = n_0;
    z = length;

    initialize_matrix();
}

void FreeSpace::initialize_matrix()
{
    A = 1;
    B = z/n0;
    C = 0;
    D = 1;
}

std::complex<double> FreeSpace::propagate(std::complex<double> qin)
{
    return calc_q_out(qin, 1);
}

void FreeSpace::set_length(double length)
{
    z = length;
    initialize_matrix();
}

void FreeSpace::set_n0(double n_0)
{
    n0 = n_0;
    initialize_matrix();
}