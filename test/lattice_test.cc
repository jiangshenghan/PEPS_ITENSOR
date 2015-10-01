
#include "lattice.h"

int main()
{
    Square_Lattice_Cylinder square{std::array<int,2>{4,4}};
    square.print_lattice_inf();
    return 0;
}
