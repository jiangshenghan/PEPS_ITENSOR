
#ifndef _IPEPS_H_
#define _IPEPS_H_

//
//ipeps class for simple update, which has diagonal tensors on links
//
class iPEPS_Simple_Update
{
    public:
        //
        //Constructors
        //
        iPEPS_Simple_Update(int& d, int& D, const iLattice_Base& ilattice);

        //
        //Destructors
        //
        ~iPEPS_Simple_Update();

        //
        //Access Methods
        //

    private:
        //
        //Data Member
        //
        //d_ is physical index dimension while D_ is virtual leg dim.
        int d_, D_;
        const iLattice_Base& ilattice_;
        vector<Index> site_physical_legs_;
        vector< vector<Index> > site_virtual_legs_;
        vector< array<Index,2> > link_legs_;
        vector<ITensor> site_tensors_, link_tensors_;
}

#endif
