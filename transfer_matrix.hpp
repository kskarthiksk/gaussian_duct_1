#ifndef TRANSFER_MATRIX_H
#define TRANSFER_MATRIX_H

#include <complex>

// class GaussianBeamFormulas {
//     public:
//     /**
//      * @brief Calculates the beam waist
//      * 
//      * Calculates beam waist for the given q factor and wavelength at a distance z
//      * 
//      * @param z [m], distance to beam waist
//      * 
//      * @param wavelength [m], wavelength of light
//      * 
//      * @param q q factor of beam
//      * 
//      * @return beam radius at z (metres)
//      */
//     static double calc_beam_radius_at_zero(double wavelength, std::complex<double> q);          // internal function to calculate beam waist
// };



/**
 * @brief Class implementing transfer matrix method for a gaussian duct
 * 
 * Stores parameters for the gaussian duct and calculates q_out from q_in
 */
class TransferMatrix
{
    double n0;           // refractive index at the center of Gaussian duct
    double n2;           // n(x) = n0 - 1/2 * n2 * x^2

   std::complex<double> q_input;    // q parameter at input
  std::complex<double> q_output;   // q parameter at output

    /**
     * @brief Internal function to calculate q_out
     * 
     * Takes in q and calculates q_out after a distance z through the Gaussian duct
     * 
     * @param z distance of propagation or length of Gaussian duct (metres)
     * 
     * @param qin q parameter at input
     * 
     * @return q parameter at output
     */
    std::complex<double> get_q_at(double z, std::complex<double> qin);

    /**
     * @brief Calculates the beam waist
     * 
     * Calculates beam waist for the given q factor and wavelength at a distance z
     * 
     * @param z [m], distance to beam waist
     * 
     * @param wavelength [m], wavelength of light
     * 
     * @param q q factor of beam
     * 
     * @return beam radius at z (metres)
     */
    static double calc_beam_radius_at_zero(double wavelength, std::complex<double> q);          // internal function to calculate beam waist


public:
    /**
     * @brief Constructor
     * 
     * creates new TransferMatrix object taking in n0 and n2
     * 
     * @param n0 Refractive index at the center of Gaussian duct
     * 
     * @param n2 n(x) = n0 - 1/2 * n2 * x^2
     */
    TransferMatrix(double n0, double n2);

    /**
     * @brief Constructor
     * 
     * creates new TransferMatrix object taking in n0, n2 and q_in
     * 
     * @param n0 Refractive index at center
     * 
     * @param n2 n(x) = n0 - 1/2 * n2 * x^2
     * 
     * @param q_in q parameter at input
     */
    TransferMatrix(double n0, double n2, std::complex<double> q_in);

    /**
     * @brief Propagate a distance z through the Gaussian duct
     * 
     * Takes in q_in and calculates q_out after a distance z
     * 
     * @param z distance of propagation or length of Gaussian duct (metres)
     * 
     * @param q q parameter at input
     * 
     * @return q parameter at output
     */
    std::complex<double> propagate(double z, std::complex<double> q);

    /**
     * @brief Propagate a distance z through the Gaussian duct
     * 
     * Uses current q_in to calculates q_out after a distance z
     * 
     * @param z distance of propagation or length of Gaussian duct (metres)
     * 
     * @return q parameter at output
     */
    std::complex<double> propagate(double z);

    /**
     * @brief Calculates the beam waist
     * 
     * Calculates beam waist for the given q factor and wavelength at a distance z
     * 
     * @param z distance to beam waist (metres)
     * 
     * @param wavelength wavelength of light (metres)
     * 
     * @param q q factor of beam
     * 
     * @return beam waist (metres)
     */
    double beam_waist(double z, double wavelength, std::complex<double> q);

    /**
     * @brief Calculates the beam waist
     * 
     * Calculates beam waist for q_out and given wavelength at a distance z
     * 
     * @param z distance to beam waist (metres)
     * 
     * @param wavelength wavelength of light (metres)
     * 
     * @return beam waist (metres)
     */
    double beam_waist(double z, double wavelength);
};

#endif