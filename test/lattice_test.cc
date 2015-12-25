
#include "lattice.h"

int main()
{
    Square_Lattice_Torus square_torus({3,3});
    square_torus.print_lattice_inf();

    //Square_Lattice_Cylinder square_cylinder({3,3}); 
    //square_cylinder.print_lattice_inf();

    //Square_Lattice_Open square_open({3,3});
    //square_open.print_lattice_inf();

    //Kagome_Cirac_Lattice_Torus kagome_torus({3,3});
    //kagome_torus.print_lattice_inf();

    std::stringstream ss;
    ss << "/home/jiangsb/code/peps_itensor/result/test/lattice_test"; 
    std::string file_name=ss.str();
    writeToFile(file_name,square_torus);
    //std::ofstream s(file_name.c_str());
    //square_torus.write(s);
    //s.close();

    
    Lattice_Base lattice_from_file;
    readFromFile(file_name,lattice_from_file);
    //std::ifstream si(file_name.c_str());
    //lattice_from_file.read(si);
    //s.close();

    lattice_from_file.print_lattice_inf();

    return 0;
}
