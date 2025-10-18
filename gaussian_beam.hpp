#ifndef GAUSSIAN_BEAM_H
#define GAUSSIAN_BEAM_H

#include <complex>

class GaussianBeam
{
public:
    /**
     * @brief Calculates the beam waist at zero
     * 
     * Calculates the beam waist at zero for the given q factor and wavelength
     * 
     * @param wavelength [m] wavelength of light
     * 
     * @param q complex beam parameter
     * 
     * @return w0 [m] beam waist at zero
     */
    static double calc_beam_waist(double wavelength, std::complex<double> q);

    /**
     * @brief Calculates beam radius
     * 
     * Calculates the beam radius for the given q factor and wavelength
     * 
     * @param wavelength [m] wavelength of light
     * 
     * @param q complex beam parameter
     * 
     * @return w [m] beam radius
     */
    static double calc_beam_radius(double wavelength, std::complex<double> q);
};

#endif