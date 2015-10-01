
#include "peps.h"

int main()
{
    //Square_Lattice_Torus square_lattice{std::array<int,2>{2,2}};
    Square_Lattice_Cylinder square_lattice{std::array<int,2>{6,2}};
   
    square_lattice.print_lattice_inf();


    //PEPS_IndexSet index_set(2,2,square_lattice);
    //PEPS peps_test(square_lattice,index_set);

    IQPEPS_IndexSet_SpinHalf index_set(6,square_lattice);
    //IQPEPS_IndexSet_SpinHalf index_set(17,std::vector<int>{1,3,0,0,2},square_lattice);

    cout << index_set.name() << endl << endl;
    cout << "Physical Legs: " << endl;
    for (const auto &leg : index_set.phys_legs())
    {
        cout << leg << endl;
    }
    cout << "Virtual Legs: " << endl;
    for (const auto &leg : index_set.virt_legs())
    {
        cout << leg << endl;
    }

    IQPEPS peps_test(square_lattice,index_set);

    cout << peps_test.name() << endl << endl;

    cout << "Site Tensors: " << endl;
    for (const auto &tensor : peps_test.site_tensors())
    {
        cout << tensor << endl;
    }

    cout << "Bond Tensors: " << endl;
    for (const auto &tensor : peps_test.bond_tensors())
    {
        cout << tensor << endl;
    }

    cout << "Boundary Tensors: " << endl;
    for (const auto &one_side_boundary_tensors : peps_test.boundary_tensors())
    {
        for (const auto &tensor : one_side_boundary_tensors)
        {
            cout << tensor << endl;
        }
    }

    return 0;
}
