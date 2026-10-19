#include "beam_through_material.hpp"
#include "gaussian_beam.hpp"

#include <complex>
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
    
    n2.real(params["n2_real"]);
    
    propagation_length = params["zScreen"];
    wavelength_0 = params["wavelength"];
    n0.imag((params["g0"] - params["a0"]) * wavelength_0 / (2 * M_PI));
    n2.imag((params["g2"] - params["a2"]) * wavelength_0 / (2 * M_PI));
    if (n2 != std::complex<double>(0,0))
        material_type = 1;
    steps = params["steps"];
    angle = params["angle"];
    amp_0 = params["amplitude"];
    // w0 = params["w0"];
    distance_from_waist = -params["zw"];
    step_size = propagation_length/steps;
    if(step_size > 1)
    {
        step_size = 1;
        steps = propagation_length;
    }
    wavelength = wavelength_0/n0.real();

    angle = asin(sin(angle)/n0.real());
    w0 = wavelength / (M_PI * tan(angle));
    if(material_type == 0)
        n2 = 0;

    // Rayleigh length (m)
    double zR = M_PI * w0 * w0 / wavelength;

    qin = distance_from_waist + 1i * zR;
}

void BeamThroughMaterial::run(ofstream& fout)
{
    // TransferMatrix object
    Material material = Material(n0, n2, wavelength_0, step_size, material_type, amp_0);

    complex<double> q_out = qin;    // output q parameter after each iteration
    double w;                       // [m] beam waist after each iteration
    double r;                       // [m] beam radius after each iteration
    double amp = amp_0;             // amplitude after each iteration
    double fwhm;                    // fwhm after each iteration

    // adding column titles
    fout<<"q param real,q param imag,beam waist,beam radius,fwhm,amplitude\n";


    int limit {numeric_limits<double>::digits10};
    fout<<fixed<<setprecision(limit);

    r = GaussianBeam::calc_beam_radius(wavelength, q_out);
    fwhm = GaussianBeam::calc_fwhm(r);

    // adding initial values to file
    fout<<q_out.real()<<','<<q_out.imag()<<','<<w0<<','<<r<<','<<fwhm<<','<<amp<<'\n';
    
    for(int i = 0; i < steps; ++i)
    {
        amp = material.calc_amplitude(q_out);
        q_out = material.propagate(q_out);
        w = GaussianBeam::calc_beam_waist(wavelength, q_out);
        r = GaussianBeam::calc_beam_radius(wavelength, q_out);
        fwhm = GaussianBeam::calc_fwhm(r);
        fout<<q_out.real()<<','<<q_out.imag()<<','<<w<<','<<r<<','<<fwhm<<','<<amp<<'\n';
    }
}