
#include "lattice.h"

int main()
{
    //Square_Lattice_Torus square_torus({3,3});
    //square_torus.print_lattice_inf();

    //Square_Lattice_Cylinder square_cylinder({3,3}); 
    //square_cylinder.print_lattice_inf();

    //Square_Lattice_Open square_open({3,3});
    //square_open.print_lattice_inf();

    Kagome_Lattice_Cirac_Torus kagome_torus({3,3});
    kagome_torus.print_lattice_inf();
    return 0;
}
