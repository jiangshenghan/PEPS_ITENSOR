
#ifndef _SINGLE_LAYER_ENV_H_
#define _SINGLE_LAYER_ENV_H_

#include "peps.h"

//
//obtain single layer environement tensors of some symmetric PEPS
//
class SL_Env_Tensors
{
    public:
        SL_Env_Tensors(const std::vector<IQTensor> &site_tensors, double cut_off=1e-10);

        void update_env_tensor_one_step();

    private:
        IQIndex comm_leg_;
        //site_tensors_ formed by two tensors with one common leg
        //env_tensors_ are formed by tensors with two legs: dag(link) and prime(link), where link is out virtual leg of one site tensor. These tensors approx single layer env
        std::vector<IQTensor> site_tensors_, env_tensors_;
        //spectrum_sqrt_ stores sqaure root of entanglement between two single tensors
        std::vector<double> spectrum_sqrt_;
        //cut_off_ control the env precision
        double cut_off_;
};
//



#endif
