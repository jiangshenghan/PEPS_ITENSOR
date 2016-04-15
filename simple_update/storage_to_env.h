#ifndef _STORAGE_TO_ENV_H
#define _STORAGE_TO_ENV_H

#include "tnetwork_storage.h"

template <class Tensor> class Storage_to_Env{
public:
  Storage_to_Env(const Tnetwork_Storage<Tensor> & input, itpp::Vec<int> options);
  virtual void update_env(const Tnetwork_Storage<Tensor> & new_input, itpp::Vec<int> options);
  virtual Tensor output_whole_env() const;
  virtual itpp::Array<Tensor> output_env_MPOs() const;
};

#endif