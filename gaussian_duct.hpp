#ifndef GAUSSIAN_DUCT_H
#define GAUSSIAN_DUCT_H

#include "transfer_matrix.hpp"

/**
 * @brief Class implementing a gaussian duct
 * 
 * Stores n0 and n2 parameters for the gaussian duct
 */
class GaussianDuct : public TransferMatrix
{
    double n0;  // refractive index at the center of Gaussian duct
    double n2;  // n(x) = n0 - 1/2 * n2 * x^2
    double z;   // [m] Length of the gaussian duct

    void initialize_matrix();

    public:
    /**
     * @brief Constructor
     * 
     * Creates new GaussianDuct object taking in n0 and n2
     * 
     * Also initializes the transfer matrix parameters A, B, C and D using parameters of gaussian duct
     * 
     * @param n_0 Refractive index at the center of Gaussian duct
     * 
     * @param n_2 n(x) = n0 - 1/2 * n2 * x^2
     * 
     * @param length [m] Length of the gaussian duct
     */
    GaussianDuct(double n_0, double n_2, double length);

    /**
     * @brief Function to calculate q_out after the gaussian duct
     * 
     * Takes in qin and calculates q_out after a distance z through the Gaussian duct
     * 
     * @param qin q parameter at input
     * 
     * @return q parameter at output
     */
    std::complex<double> propagate(std::complex<double> qin);

    /**
     * @brief Change the length of propagation in the gaussian duct
     * 
     * @param length [m] new distance of propagation or length of Gaussian duct
     */
    void set_length(double length);

    /**
     * @brief Change refractive index at the center of Gaussian duct
     * 
     * @param n_0 new value for refractive index at the center of Gaussian duct
     */
    void set_n0(double n_0);

    /**
     * @brief Change n2 (n(x) = n0 - 1/2 * n2 * x^2)
     * 
     * @param n_2 n(x) = n0 - 1/2 * n2 * x^2
     */
    void set_n2(double n_2);
};

#endif