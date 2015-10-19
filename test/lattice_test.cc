
#include "lattice.h"

int main()
{
    //Square_Lattice_Cylinder square{std::array<int,2>{4,4}};
    //square.print_lattice_inf();

    Square_Lattice_Open square(std::array<int,2>{4,3});
    square.print_lattice_inf();

    return 0;
}
