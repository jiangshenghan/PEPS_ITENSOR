
#ifndef _TROTTER_H_
#define _TROTTER_H_

//
//class for trotter gate for NN Heisenberg Hamiltonian exp(-tH)
//
class NN_Heisenberg_Trotter_Gate
{
    public:
        //
        //Constructors
        //
        NN_Heisenberg_Trotter_Gate() {}

        NN_Heisnberg_Trotter_Gate(std::array<IQIndex,2> nn_sites, int t=0.1);


    private:
        //
        //Data Member
        //
        std::array<std::array<IQIndex,2>,2> phys_legs_;
        std::array<IQIndex,2> virt_legs_;
        std::array<IQTensor,2> site_tensors_;
        IQTensor bond_tensor_;
        //t_ labels imaginary time
        int t_;
}

#endif
