
#include "peps.h"

int main()
{
    Square_Lattice_Torus square_lattice{std::array<int,2>{2,2}};
   
    square_lattice.print_lattice_inf();


    //PEPS_IndexSet index_set(2,3,square_lattice);
    //PEPS_Torus peps_test(square_lattice,index_set);

    //for (const auto &tensor : peps_test.site_tensors())
    //{
    //    //cout << tensor << endl;
    //    PrintDat(tensor);
    //}


    //IQPEPS_IndexSet_SpinHalf index_set(15,square_lattice);
    IQPEPS_IndexSet_SpinHalf index_set(17,std::vector<int>{1,3,0,0,2},square_lattice);
    IQPEPS_Torus peps_test(square_lattice,index_set);

    for (const auto &tensor : peps_test.site_tensors())
    {
        cout << tensor << endl;
    }

    return 0;
}
