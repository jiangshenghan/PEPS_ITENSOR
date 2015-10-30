#ifndef _TNETWORK_STORAGE_H
#define _TNETWORK_STORAGE_H

#include "itpp/itbase.h"
#include "core.h"

template <class Tensor>  class Tnetwork_Storage{
public:
  //label the type of this tensor network
  //let us define: 
  //_tnetwork_type=1 for square square, 
  //_tnetwork_type=2 for honeycomb lattice,
  //_tnetwork_type=3 for triangular lattice
  //_tnetwork_type=4 for kagome lattice,
  //_tnetwork_type=5 for kagome Cirac lattice,...
  int _tnetwork_type;
  //label the boundary condition of the tensor network,
  //let us define:
  //_boundary_condition=1 for torus,
  //_boundary_condition=2 for cylinder periodic along y direction(i.e., (x,y,n)=(x,y+Ly,n))
  //_boundary_condition=3 for open boundary in both directions
  int _boundary_condition;
  //the number of sublattices in one unit cell
  //e.g.: _n_subl=1 for square lattice, _n_subl=2 for honeycomb lattice
  int _n_subl;
  //the number of unit cells along x-direction
  int _Lx;
  //the number of unit cells along y-direction
  int _Ly;
  //store the list of all site tensors (single-layer).
  itpp::Array<Tensor> _tensor_list;
  //convert the site-coordinate: (x,y,n) to the
  //index of that site in _tensor_list=_coor_to_siteind(x,y)(n)
  itpp::Mat< itpp::Vec<int> > _coor_to_siteind; 
  
  //read to disk file:
  void read(std::istream& s){
    s.read((char *) &_tnetwork_type, sizeof(_tnetwork_type));
    s.read((char *) &_boundary_condition, sizeof(_boundary_condition));
    s.read((char *) &_n_subl, sizeof(_n_subl));
    s.read((char *) &_Lx, sizeof(_Lx));
    s.read((char *) &_Ly, sizeof(_Ly));
    int n_sites=_n_subl*_Lx*_Ly;
    _tensor_list.set_size(n_sites);
    for(int i=0;i<_tensor_list.size();i++){
      _tensor_list(i).read(s);
    }
    _coor_to_siteind.set_size(_Lx,_Ly);
    for(int i=0;i<_coor_to_siteind.rows();i++){
      for(int j=0;j<_coor_to_siteind.cols();j++){
	_coor_to_siteind(i,j).set_size(_n_subl);
	for(int k=0;k<_coor_to_siteind(i,j).size();k++){
	  s.read((char *) &(_coor_to_siteind(i,j)(k)),sizeof((_coor_to_siteind(i,j)(k))));
	}
      }
    }
    return;
  }
  
  //write from disk file.
  void write(std::ostream& s) const{
    s.write((char *) &_tnetwork_type, sizeof(_tnetwork_type));
    s.write((char *) &_boundary_condition, sizeof(_boundary_condition));
    s.write((char *) &_n_subl, sizeof(_n_subl));
    s.write((char *) &_Lx, sizeof(_Lx));
    s.write((char *) &_Ly, sizeof(_Ly));
    for(int i=0;i<_tensor_list.size();i++){
      _tensor_list(i).write(s);
    }
    for(int i=0;i<_coor_to_siteind.rows();i++){
      for(int j=0;j<_coor_to_siteind.cols();j++){
	for(int k=0;k<_coor_to_siteind(i,j).size();k++){
	  s.write((char *) &(_coor_to_siteind(i,j)(k)),sizeof((_coor_to_siteind(i,j)(k))));
	}
      }
    }
    return;
  }
  
};

#endif