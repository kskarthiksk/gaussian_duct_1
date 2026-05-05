#ifndef TRANSFER_MATRIX_H
#define TRANSFER_MATRIX_H

#include <complex>

/**
 * @brief Class implementing transfer matrix method for a gaussian duct
 * 
 * Stores parameters for the gaussian duct and calculates q_out from q_in
 */
class TransferMatrix
{
protected:
    std::complex<double> A, B, C, D;   // elements of ABCD matrix

public:
    /**
     * @brief Function to calculate q_out
     * 
     * Takes in q and calculates q_out after a distance z through the Gaussian duct
     * 
     * @param  n0 (complex) refractive index at the center of Gaussian duct
     * 
     * @param qin q parameter at input
     * 
     * @return q parameter at output
     */
    std::complex<double> calc_q_out(std::complex<double> qin, std::complex<double> n0);
};

#endif