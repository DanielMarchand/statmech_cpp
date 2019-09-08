// These functions should later be moved to routines.cpp
#include <fstream>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <random>
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

template <typename Derived>
void write_xyz(std::string outputxyz, int &istep, 
               int &n_atoms, const Eigen::MatrixBase<Derived> &atoms){
    std::ofstream outfile (outputxyz, std::ios_base::app);
    if (outfile.is_open()){
        outfile << n_atoms << std::endl; 
        outfile << "# Step" << istep << std::endl;
        for (int i=0; i<n_atoms; i++){
            outfile << "Ar   " 
                    << atoms.row(i).col(0) << "  "
                    << atoms.row(i).col(1) << "  "
                    << atoms.row(i).col(2) << "  "
                    << std::endl;
        }
        outfile.close();
    }
    else std::cout << "Unable to open file";
}

//  subroutine adjustqcm(natoms,q)
//    ! put the origin in the CM
//    implicit none
//    integer, intent(in)  :: natoms
//    double precision, intent(inout)  :: q(3,natoms)
//    double precision :: qcm(3)
//    integer i
//    qcm=0.0d0
//    ! compute the position of CM
//    do i=1,natoms
//      qcm=qcm+q(:,i)
//    enddo
//    qcm=qcm/natoms
//    ! traslate
//    do i=1,natoms
//      q(:,i)=q(:,i)-qcm
//    enddo   
//  end subroutine adjustqcm

//Compute the potential energy
template <typename Derived>
void get_U(int n_atoms, const Eigen::MatrixBase<Derived> &atoms, double &potential_e){
    potential_e = 0.0;
    double d_2, d_6;
    for (int i=0; i<n_atoms; i++){
        for (int j=i+1; j<n_atoms; j++){
            // The pow function is very very slow, so we avoid it!
            d_2 = (atoms.row(i)-atoms.row(j)).squaredNorm();
            d_6 = d_2*d_2*d_2;
            potential_e = potential_e + 1.0/(d_6*d_6) - 1.0/(d_6);
        }
    }
    potential_e = potential_e*4.0;
}
