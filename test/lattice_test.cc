
#include "lattice.h"

int main()
{
    Square_Lattice_Torus square{std::array<int,2>{3,3}};
    square.print_lattice_inf();
    return 0;
}
