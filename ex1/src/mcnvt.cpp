#include <iostream>
#include <sstream>
#include <vector>
// #include "routines.hh"  

//SP4E built-in types
//UInt = unsigned int;
//Real = double;
//complex = std::complex<Real>;

// These functioins should later be moved to routines.cpp
#include <fstream>
#include <Eigen/Dense>
#include <cmath>
struct Atom {
    std::string element;
    //double x, y, z;
    Eigen::Vector3d coordinate;
};
void readdata(std::string dataxyz, std::vector<Atom> &atoms, int &natoms){


    std::ifstream input_stream(dataxyz.c_str());
    std::string input_line;
    if (input_stream.is_open() == false){
        std::cerr << "Cannot open file" << dataxyz << std::endl;
        throw;
    }

    while (input_stream.good()){
        getline(input_stream, input_line);
        //TODO: figure out a way to read the natoms separately 
        if (input_line[0] == '#'|| input_line.size() == 0){
            continue;
        }

        std::stringstream inputline_stream(input_line);

        Atom this_atom;
        inputline_stream >> this_atom.element 
        >> this_atom.coordinate[0] >> this_atom.coordinate[1] >> this_atom.coordinate[2];
        atoms.push_back(this_atom);
    }
    input_stream.close();
    natoms = atoms.size();
}
//Compute the potential energy
inline void get_U(int natoms, std::vector<Atom> &atoms, double &potential_e){

    potential_e = 0.0;
    Eigen::Vector3d distance;
    double distance_norm;
    for (int i=0; i<natoms; i++){
        for (int j=i+1; j<natoms; j++){
            distance = atoms[i].coordinate - atoms[j].coordinate;
            distance_norm = distance.norm();
            potential_e = potential_e + 
                          (1.0/pow(distance_norm,12) -
                           1.0/pow(distance_norm, 6));
        }
    }
    potential_e = potential_e*4.0;
}



int main(int argc, char **argv)
{
    // Dummy arguments for now
    int seed = 1357;
    double temp = 0.23;
    std::string dataxyz = "input_dummy.xyz";
    //int nstep = 10000000;
    int nstep = 5000;
    int stridetrj = 50000;
    int stridelog = 100;
    double mcstep = 0.01;
    std::string outputf = "output_dummy.xyz";

    //TODO: enable command-line arguments later
    //if (argc != 9)
    //{
    //    std::cout << "Usage: " << argv[0] << " seed temp dataxyz nstep stridetrj stridelog mcstep outputf " << std::endl;
    //    std::exit(EXIT_FAILURE);
    //}

    //// Read input arguments
    //int seed;
    //std::stringstream(argv[1]) >> seed;
    //double temp;
    //std::stringstream(argv[2]) >> temp;
    //std::string dataxyz = argv[3];
    //int nstep;
    //std::stringstream(argv[4]) >> nstep;
    //int stridetrj;
    //std::stringstream(argv[5]) >> stridetrj;
    //std::string stridelog = argv[6];
    //double mcstep;
    //std::stringstream(argv[7]) >> mcstep;
    //std::string outputf = argv[8];
    //// beta is derived from temp
    //double beta;
    //beta = 1.0/temp;
    ////TODO: initialize random number generator

    // Read the box and the number of atoms
    std::vector<Atom> atoms;
    int natoms;
    readdata(dataxyz, atoms, natoms);
    for (int i = 0; i < natoms; i++){
        Atom atom = atoms[i];
        std::cout << atom.element << " "
        << atom.coordinate[0] << " " 
        << atom.coordinate[1] << " " 
        << atom.coordinate[2] << std::endl;
    }

    // Start MC code
    int attempts = 0;
    int accepted = 0;
    double potential_e;
    double old_potential_e;
    double new_potential_e;
    // Compute potential energy of the initial configuration
    get_U(natoms, atoms, potential_e);
    std::cout << "potential_e: " << potential_e << std::endl;
    for (int istep=0; istep<nstep; istep++){
        attempts = attempts + 1;
        // Save the potential energy before the move
        old_potential_e = potential_e;
        // Select a particle randomly
        //TODO
        // Save the original position
        //TODO
        // Randomly displace the selected particles
        //TODO
        // Compute potential energy after the move
        std::cout << "istep: " << istep << std::endl;
        get_U(natoms, atoms, new_potential_e);

    }
    //std::cout << std::log(0.0) << std::endl;
    std::cout << errno << std::endl;
}
