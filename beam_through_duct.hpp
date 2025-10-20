#ifndef BEAM_THROUGH_DUCT_H
#define BEAM_THROUGH_DUCT_H

#include "gaussian_duct.hpp"
#include "gaussian_beam.hpp"

#include <fstream>

using namespace std;

class BeamThroughDuct
{
    double n0, n2;              // refractive indices
    double duct_length;         // [m] length of gaussian duct
    double wavelength;          // [m] wavelength of light
    int iter;                   // number of iterations
    double w0;                  // [m] beam waist
    double distance_from_waist; // [m] location relative to beam waist
    complex<double> qin;        // q parameter of the beam

public:
    /**
     * @brief Constructor
     * 
     * Takes in a file object and initializes class member variables
     * 
     * @param file File object containing variables
     */
    BeamThroughDuct(ifstream& fin);

    /**
     * @brief Run simulation
     * 
     * Runs simulation for number of iterations set by `iter` variable
     * 
     * @param fout File object to write simulation results into
     */
    void run(ofstream& fout);
};

#endif