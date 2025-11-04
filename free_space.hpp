#ifndef FREE_SPACE_H
#define FREE_SPACE_H

#include "transfer_matrix.hpp"

/**
 * @brief Class implementing a length of free space
 * 
 * Stores n0 and z parameters for the free space
 */
class FreeSpace : public TransferMatrix
{
    double n0;  // refractive index of free space
    double z;   // [m] Length of the free space

    void initialize_matrix();

    public:
    /**
     * @brief Constructor
     * 
     * Creates new GaussianDuct object taking in n0 and n2
     * 
     * Also initializes the transfer matrix parameters A, B, C and D using parameters of gaussian duct
     * 
     * @param n_0 Refractive index of free space
     * 
     * @param length [m] Length of the free space
     */
    FreeSpace(double n_0, double length);

    /**
     * @brief Function to calculate q_out after the free space
     * 
     * Takes in qin and calculates q_out after a distance z through the free space
     * 
     * @param qin q parameter at input
     * 
     * @return q parameter at output
     */
    std::complex<double> propagate(std::complex<double> qin);

    /**
     * @brief Change the length of propagation in the free space
     * 
     * @param length [m] new distance of propagation or length of free space
     */
    void set_length(double length);

    /**
     * @brief Change refractive index of free space
     * 
     * @param n_0 new value for refractive index of free space
     */
    void set_n0(double n_0);
};

#endif