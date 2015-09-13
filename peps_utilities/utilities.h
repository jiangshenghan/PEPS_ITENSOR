
#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include "predef.h"

//Spin stores total spin quantum number (times two) as well as the direction
struct Spin
{
    public:
        Spin() {}

        Spin(int spin_qn, Arrow dir=Out) : spin_qn_(spin_qn), dir_(dir) { }

        int spin_qn() const { return spin_qn_; }
        const Arrow &dir() const { return dir_; }
        
        void dag() { dir_=-dir_; }

    private:
        int spin_qn_;
        Arrow dir_;
};

//class IndexSpin is an IQIndex which stores indices with fixed spin quantum number
class IndexSpin 
{
    public:
        //
        //Constructor
        //
        IndexSpin() {}

        explicit IndexSpin(const IQIndex &iqind) : spin_(iqind.nindex()-1,iqind.dir()), sz_legs_(iqind) { assert(valid()); }

        IndexSpin(const IQIndex &iqind, const Spin &spin) : spin_(spin), sz_legs_(iqind) { assert(valid()); }

        //construct IndexSpin from spin quantum number, no extra deg
        explicit IndexSpin(const Spin &spin): spin_(spin)
        {
            std::vector<IndexQN> sz_indqns;
            for (int sz=spin_.spin_qn(); sz>=-(spin_.spin_qn()); sz-=2)
            {
                std::stringstream ss;
                ss << "S=" << spin_.spin_qn()/2.0 << ", Sz=" << sz/2.0 << " leg";
                std::string sz_leg_name=ss.str();
                sz_indqns.push_back(IndexQN(Index(sz_leg_name,1),QN(sz)));
            }
            sz_legs_=IQIndex(nameint("leg with 2S=", spin_.spin_qn()),sz_indqns,spin_.dir());
            
            assert(valid());
        }

        //
        //Access method
        //
        Arrow dir() const { return spin_.dir(); }

        int spin_qn() const { return spin_.spin_qn(); }
        
        const Spin &spin() const { return spin_; }

        const IQIndex &leg() const { return sz_legs_; }

        //return the first found IQIndexVal with qn=sz
        IQIndexVal indval (int sz) const
        {
            int n=0;
            for (const auto &indqn : sz_legs_.indices())
            {
                if (indqn.qn.sz()==sz) return sz_legs_.operator()(n+1);
                n+=indqn.m();
            }
            return IQIndexVal();
        }


        IQIndexVal operator()(int n) const
        {
            return sz_legs_.operator()(n);
        }
        
        int m() const { return sz_legs_.m(); }

        //
        //Method
        //
        void dag()
        {
            spin_.dag();
            sz_legs_.dag();
        }

        //
        //Constructor helpers
        //
        //To determine if the input IQIndex and spin are valid for a SU(2) representation \mathbb{D}\otime\mathbb{V}_S
        bool valid()
        {
            
            if (spin_.dir()!=sz_legs_.dir()) //check direction
                return false;

            int dim=sz_legs_[0].m();
            for (const auto &indqn : sz_legs_.indices())
            {
                int sz=indqn.qn.sz();
                if (dim!=indqn.m() || //consistent extra deg
                    std::abs(sz)>spin_.spin_qn() ||
                    std::abs(sz)%2!=spin_.spin_qn()%2) //check sz qn
                {
                    return false;
                }
            }
            return true;   //every sz has same extra deg, all int/half-int, and go from -spin,...,spin
        }

    private:
        Spin spin_;
        IQIndex sz_legs_;
};


inline std::ostream &operator<<(std::ostream &s, const Coordinate &coord)
{
    return s << "(" << coord[0] << "," << coord[1] << "," << coord[2] << ")";
}

inline std::ostream &operator<<(std::ostream &s, const Spin &spin)
{
    return s << "S=" << spin.spin_qn()/2.0 << ", dir=" << spin.dir() << endl;
}

inline std::ostream &operator<<(std::ostream &s, const IndexSpin &indspin)
{
    return s << indspin.spin() << indspin.leg();
}

#endif
