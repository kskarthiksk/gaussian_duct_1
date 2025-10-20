#include "beam_through_duct.hpp"

#include <iomanip>

using namespace std;

BeamThroughDuct::BeamThroughDuct(ifstream& fin)
{
    string str; // string to read key

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
            fin>>distance_from_waist;
    }

    // Rayleigh length (m)
    double zR = M_PI * w0 * w0 / wavelength;

    qin = distance_from_waist + 1i * zR;
}

void BeamThroughDuct::run(ofstream& fout)
{
    // TransferMatrix object
    GaussianDuct duct = GaussianDuct(n0, n2, duct_length);

    complex<double> q_out = qin;    // output q parameter after each iteration
    double w;                       // [m] beam waist after each iteration

    // adding column titles
    fout<<"q param real,q param imag,beam waist\n";


    int limit {numeric_limits<double>::digits10};
    fout<<fixed<<setprecision(limit);

    // adding initial values to file
    fout<<q_out.real()<<','<<q_out.imag()<<","<<w0*w0<<'\n';
    
    for(int i = 0; i < iter; ++i)
    {
        q_out = duct.propagate(q_out);
        w = GaussianBeam::calc_beam_waist(wavelength, q_out);
        fout<<q_out.real()<<','<<q_out.imag()<<","<<w<<'\n';
    }
}