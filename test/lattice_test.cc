
#include "lattice.h"

int main()
{
    //Square_Lattice_Cylinder square{std::array<int,2>{4,4}};
    //square.print_lattice_inf();

    Honeycomb_Lattice_Torus honeycomb(std::array<int,2>{2,2});
    honeycomb.print_lattice_inf();

    return 0;
}
