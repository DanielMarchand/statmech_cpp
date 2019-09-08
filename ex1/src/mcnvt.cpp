#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <random>
#include "routines.h"

int main(int argc, char **argv)
{
    // Dummy arguments for now

    //TODO: enable command-line arguments later
    if (argc != 9)
    {
        std::cout << "Usage: " << argv[0] << " seed temp dataxyz nstep stridetrj stridelog mcstep outputxyz " << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Read input arguments
    int seed;
    std::stringstream(argv[2]) >> seed;
    double temp;
    std::stringstream(argv[2]) >> temp;
    std::string dataxyz = argv[3];
    int nstep;
    std::stringstream(argv[4]) >> nstep;
    int stridetrj;
    std::stringstream(argv[5]) >> stridetrj;
    int stridelog;
    std::stringstream(argv[6]) >> stridelog;
    double mcstep;
    std::stringstream(argv[7]) >> mcstep;
    std::string outputxyz = argv[8];
    // beta is derived from temp
    double beta;
    beta = 1.0/temp;
    // initialize random number generator
    std::random_device rd;
    //std::mt19937 rng(rd());
    std::mt19937 rng(seed);
    // generate the gaussian function to be used later
    std::normal_distribution<double> random_gaussian(0, mcstep);

    // Read the box and the number of atoms
    int n_atoms;
    Eigen::MatrixX3d atoms = readdata(dataxyz, n_atoms);
    std::uniform_int_distribution<int> random_uniform_int(0, n_atoms-1);
    std::uniform_real_distribution<double> random_uniform_real(0, 1);

    // Start MC code
    int attempts = 0;
    int accepted = 0;
    int random_index;
    double potential_e;
    double old_potential_e;
    double new_potential_e;
    Eigen::MatrixX3d oldposition;
    // Compute potential energy of the initial configuration
    get_U(n_atoms, atoms, potential_e);
    for (int istep=0; istep<nstep; istep++){
        attempts = attempts + 1;
        // Save the potential energy before the move
        old_potential_e = potential_e;
        // Randomly select a particle for displacement
        random_index = random_uniform_int(rng);
        // Save the original position
        oldposition = atoms.row(random_index);
        // Randomly displace the selected particles
        for (int i=0; i<3; i++){
            atoms.row(random_index)[i] += random_gaussian(rng);
        }
        // Compute potential energy after the move
        get_U(n_atoms, atoms, new_potential_e);

        // Acceptance test
        double del_E;
        del_E  = new_potential_e - old_potential_e;
        if (random_uniform_real(rng) < exp(-beta*(new_potential_e-old_potential_e))){
            accepted += 1;
            potential_e = new_potential_e;
        }
        else{
            atoms.row(random_index) = oldposition;
        }
        if (istep % stridelog == 0){
            std::cout << "step: " << istep 
              << " accepted: " << accepted 
              << " E: "<< potential_e << std::endl;
        }
        if (istep % stridetrj == 0){
           //put the origin of the Cartesian coordinate system in the CM
           //call adjustqcm(natoms,q)
           write_xyz(outputxyz, istep, n_atoms, atoms);
        }

    }
}
