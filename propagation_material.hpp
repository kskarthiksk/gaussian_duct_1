#ifndef MATERIAL_H
#define MATERIAL_H

#include "transfer_matrix.hpp"

/**
 * @brief Class implementing a material
 * 
 * Stores variables denoting refractive index (or indices), length and type of material
 * 
 * Material can be free space or gaussian duct
 */
class Material : public TransferMatrix
{
    double n0;  // Central refractive index of material
    double n2;  // n(x) = n0 - 1/2 * n2 * x^2 if material is a gaussian duct (material_type = 1)
    double z;   // [m] Length of the material
    int material_type;   // Integer denoting type of material. 0: Free space 1: Gaussian duct

    /**
     * @brief Initializes the ABCD matrix according to the type of material
     */
    void initialize_matrix();

    public:
    /**
     * @brief Constructor
     * 
     * Creates new Material object taking in n0, n2, z and material_type
     * 
     * Also initializes the transfer matrix parameters A, B, C and D using parameters of material
     * 
     * @param n_0 Central refractive index of material
     * 
     * @param n_2 n(x) = n0 - 1/2 * n2 * x^2 if material is a gaussian duct (material_type = 1)
     * 
     * @param length [m] Length of the material
     * 
     * @param material Integer denoting type of material
     * 0: Free space
     * 1: Gaussian duct
     */
    Material(double n_0, double n_2, double length, int material);

    /**
     * @brief Function to calculate q_out after travelling through the material
     * 
     * Takes in qin and calculates q_out after a distance z through the material
     * 
     * @param qin q parameter at input
     * 
     * @return q parameter at output
     */
    std::complex<double> propagate(std::complex<double> qin);

    /**
     * @brief Change the length of propagation in the material
     * 
     * @param length [m] new distance of propagation or length of material
     */
    void set_length(double length);

    /**
     * @brief Change central refractive index of material
     * 
     * @param n_0 new value for central refractive index of material
     */
    void set_n0(double n_0);

    /**
     * @brief Change n2 (n(x) = n0 - 1/2 * n2 * x^2)
     * 
     * @param n_2 n(x) = n0 - 1/2 * n2 * x^2
     */
    void set_n2(double n_2);

    /**
     * @brief Change the type of material
     * 
     * @param material_type Integer denoting type of material
     * 0: Free space
     * 1: Gaussian duct
     */
    void set_material(int material_type);
};

#endif