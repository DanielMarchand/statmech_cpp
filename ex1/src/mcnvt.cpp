#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <random>
// #include "routines.hh"  

// These functions should later be moved to routines.cpp
#include <fstream>
#include <Eigen/Dense>
#include <cmath>
Eigen::MatrixX3d readdata(std::string dataxyz,  int &n_atoms){
    std::ifstream input_stream(dataxyz.c_str());
    std::string input_line;
    if (input_stream.is_open() == false){
        std::cerr << "Cannot open file" << dataxyz << std::endl;
        throw;
    }

    std::vector<double> tmp_storage;
    std::string foo_element;
    while (input_stream.good()){
        getline(input_stream, input_line);
        //TODO: figure out a way to read the n_atoms separately 
        if (input_line[0] == '#'|| input_line.size() == 0){
            continue;
        }

        std::stringstream inputline_stream(input_line);

        double x,y,z;
        inputline_stream >> foo_element
          >> x >> y >> z;
        tmp_storage.push_back(x);
        tmp_storage.push_back(y);
        tmp_storage.push_back(z);
    }
    input_stream.close();

    // store the result into an Eigen matrix
    int num_columns = 3;
    int num_rows = tmp_storage.size()/num_columns;
    Eigen::MatrixX3d atoms(num_rows, num_columns);
    int insert_row;
    int insert_col;
    for (unsigned int i=0; i<tmp_storage.size(); i++){
        insert_col = i%3;
        insert_row = i/3;
        atoms.row(insert_row).col(insert_col) << tmp_storage.at(i);
    }

    n_atoms = atoms.rows();
    return atoms;
}

//Compute the potential energy
template <typename Derived>
void get_U(int n_atoms, const Eigen::MatrixBase<Derived> &atoms, double &potential_e){
    potential_e = 0.0;
    double d_2, d_6;
    for (int i=0; i<n_atoms; i++){
        for (int j=i+1; j<n_atoms; j++){
            d_2 = (atoms.row(i)-atoms.row(j)).squaredNorm();
            d_6 = d_2*d_2*d_2;
            potential_e = potential_e + 1.0/(d_6*d_6) - 1.0/(d_6);
        }
    }
    potential_e = potential_e*4.0;
}



int main(int argc, char **argv)
{
    // Dummy arguments for now
    int seed = 1237;
    double temp = 0.23;
    std::string dataxyz = "input_dummy.xyz";
    int nstep =  10000000; // 5000
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

    }
}
