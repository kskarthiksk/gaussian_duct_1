#include "gaussian_duct.hpp"
#include "gaussian_beam.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

int main()
{
    // filestream for input parameters
    fstream fin;
    
    fin.open("params.txt", ios::in);

    double n0, n2;       // refractive indices
    double duct_length;  // length of gaussian duct (metres)
    double wavelength;   // wavelength of light (metres)
    int iter;           // number of iterations
    double w0;           // beam waist (metres)
    double z;            // location relative to beam waist (metres)

    // string to read key
    string str;

    while(!fin.eof())
    {
        fin>>str;
        if(str == "n0")
            fin>>n0;
        else if(str == "n2")
            fin>>n2;
        else if(str == "duct_length")
            fin>>duct_length;
        else if(str == "wavelength")
            fin>>wavelength;
        else if(str == "iterations")
            fin>>iter;
        else if(str == "w0")
            fin>>w0;
        else if(str == "distance_from_waist")
            fin>>z;
    }

    fin.close();    // closing input file

    // Rayleigh length (m)
    double zR = M_PI * w0 * w0 / wavelength;

    // q parameter of the beam
    complex<double> q = z + 1i * zR;

    // TransferMatrix object
    GaussianDuct duct = GaussianDuct(n0, n2, z);

    // filestream for storing output (q parameter and beam waist at output)
    fstream fout;
    fout.open("output.csv", ios::out);

    complex<double> q_out = q;   // output q parameter after each iteration
    double w;                // beam waist after each iteration (metres)

    // adding column titles
    fout<<"q param real,q param imag,beam waist\n";


    int limit {numeric_limits<double>::digits10};
    fout<<fixed<<setprecision(limit);

    // adding initial values to file
    fout<<q_out.real()<<','<<q_out.imag()<<","<<w0<<'\n';
    
    for(int i = 0; i < iter; ++i)
    {
        q_out = duct.propagate(q_out);
        w = GaussianBeam::calc_beam_waist(wavelength, q_out);
        fout<<q_out.real()<<','<<q_out.imag()<<","<<w<<'\n';
    }

    // closing output file
    fout.close();

    return 0;
}