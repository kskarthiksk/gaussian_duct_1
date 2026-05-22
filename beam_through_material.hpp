#ifndef BEAM_THROUGH_MATERIAL_H
#define BEAM_THROUGH_MATERIAL_H

#include "propagation_material.hpp"
#include "gaussian_beam.hpp"

#include <fstream>

using namespace std;

class BeamThroughMaterial
{
    complex<double> n0;         // (complex) refractive index of material
    complex<double> n2;         // n(x) = n0 - 1/2 * n2 * x^2 if material is a gaussian duct (material_type = 1)
    double propagation_length;  // [m] length of material
    double wavelength_0;        // [m] wavelength of light
    double wavelength;          // [m] wavelength_0/n0
    int iter;                   // number of iterations
    double w0;                  // [m] beam waist
    double distance_from_waist; // [m] location relative to beam waist
    complex<double> qin;        // q parameter of the beam
    int material_type;          // integer denoting type of material. 0: Free space 1: Gaussian duct
    double angle;               // [radian] half angle of divergence

public:
    /**
     * @brief Constructor
     * 
     * Takes in a file object and initializes class member variables
     * 
     * @param file File object containing variables
     */
    BeamThroughMaterial(ifstream& fin);

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