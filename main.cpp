#include "beam_through_material.hpp"

#include <fstream>
#include <iostream>


using namespace std;

int main()
{
    ifstream fin;   // filestream for input parameters
    fin.open("params1.txt", ios::in);
    BeamThroughMaterial simulation = BeamThroughMaterial(fin);
    fin.close();    // closing input file

    // filestream for storing output (q parameter and beam waist at output)
    ofstream fout;
    fout.open("output.csv", ios::out);
    simulation.run(fout);   // run simulation
    fout.close();

    return 0;
}