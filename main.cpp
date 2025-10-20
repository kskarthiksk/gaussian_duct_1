#include "gaussian_duct.hpp"
#include "gaussian_beam.hpp"
#include "beam_through_duct.hpp"

#include <fstream>
#include <iostream>


using namespace std;

int main()
{
    ifstream fin;   // filestream for input parameters
    fin.open("params.txt", ios::in);
    BeamThroughDuct simulation = BeamThroughDuct(fin);
    fin.close();    // closing input file

    // filestream for storing output (q parameter and beam waist at output)
    ofstream fout;
    fout.open("output.csv", ios::out);
    simulation.run(fout);   // run simulation
    fout.close();

    return 0;
}