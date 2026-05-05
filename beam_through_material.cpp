#include "beam_through_material.hpp"

#include <iomanip>
#include <map>
#include <sstream>

using namespace std;

BeamThroughMaterial::BeamThroughMaterial(ifstream& fin)
{

    string line;                // variable to store a line in the param file
    map<string, double> params; // map to store parameters as key-value pairs

    // reading file line by line and storing parameters into the map
    while(getline(fin, line))
    {
        istringstream parameter(line);
        string key;
        parameter>>key;
        if(key[0] == '#' || key[0] == '\0')
            continue;
        double value;
        parameter>>value;
        params[key] = value;
    }

    // saving parameters into member variables
    material_type = params["material"];
    n0.real(params["n0_real"]);
    n0.imag(params["n0_imag"]);
    n2.real(params["n2_real"]);
    n2.imag(params["n2_imag"]);
    propagation_length = params["propagation_length"];
    wavelength = params["wavelength"];
    iter = params["iterations"];
    w0 = params["w0"];
    distance_from_waist = params["distance_from_waist"];
    wavelength /= n0.real();

    if(material_type == 0)
        n2 = 0;

    // Rayleigh length (m)
    double zR = M_PI * w0 * w0 / wavelength;

    qin = distance_from_waist + 1i * zR;
}

void BeamThroughMaterial::run(ofstream& fout)
{
    // TransferMatrix object
    Material material = Material(n0, n2, wavelength, propagation_length, material_type);

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