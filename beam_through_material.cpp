#include "beam_through_material.hpp"

#include <iomanip>

using namespace std;

BeamThroughMaterial::BeamThroughMaterial(ifstream& fin)
{
    string str; // string to read key

    while(!fin.eof())
    {
        fin>>str;
        if(str == "material")
            fin>>material_type;
        else if(str == "n0")
            fin>>n0;
        else if(str == "n2")
            fin>>n2;
        else if(str == "propagation_length")
            fin>>propagation_length;
        else if(str == "wavelength")
            fin>>wavelength;
        else if(str == "iterations")
            fin>>iter;
        else if(str == "w0")
            fin>>w0;
        else if(str == "distance_from_waist")
            fin>>distance_from_waist;
    }

    if(material_type == 0)
        n2 = 0;

    // Rayleigh length (m)
    double zR = M_PI * w0 * w0 / wavelength;

    qin = distance_from_waist + 1i * zR;
}

void BeamThroughMaterial::run(ofstream& fout)
{
    // TransferMatrix object
    Material material = Material(n0, n2, propagation_length, material_type);

    complex<double> q_out = qin;    // output q parameter after each iteration
    double w;                       // [m] beam waist after each iteration
    double r;                       // [m] beam radius after each iteration

    // adding column titles
    fout<<"q param real,q param imag,beam waist,beam radius\n";


    int limit {numeric_limits<double>::digits10};
    fout<<fixed<<setprecision(limit);

    r = GaussianBeam::calc_beam_radius(wavelength, q_out);

    // adding initial values to file
    fout<<q_out.real()<<','<<q_out.imag()<<','<<w0*w0<<','<<r<<'\n';
    
    for(int i = 0; i < iter; ++i)
    {
        q_out = material.propagate(q_out);
        w = GaussianBeam::calc_beam_waist(wavelength, q_out);
        r = GaussianBeam::calc_beam_radius(wavelength, q_out);
        fout<<q_out.real()<<','<<q_out.imag()<<','<<w<<','<<r<<'\n';
    }
}