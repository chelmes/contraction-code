/***************************************************************************/
//
// Short demonstration of gaugefield class
//
// Author: Christopher Helmes
//
/***************************************************************************/
#include <array>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <typeinfo>
#include <vector>

#include "CorrelatorIo2pt.h"
//#include "BasicOperator.h"
//#include "GlobalData.h"
#include "ranlxs.h"
//#include "typedefs.h"

int main(int ac, char* av[]) {

// dummy data
  // initialize one OpInfo struct
  pdg op_info;
  // one three momentum
  std::array<int,3> mom {{0,1,2}};
  // one displacement vector
  std::array<int,3> der {{3,4,5}};
  // gamma configuration 
  std::array<int,4> dirac {{5,4,4,4}};
  // initialize operator info
  op_info.p = mom;
  op_info.dis = der;
  op_info.gamma = dirac;
  Tag so_si;
  for (int i = 0; i < 2; ++i){
    so_si.mom[i] = i+1;
    so_si.dis[i] = i+2;
    so_si.gam[i] = i+3;
  }

  // print everything as a check
  std::cout << "p:"; 
  for (auto mom_comp : op_info.p) 
    std::cout << " " << mom_comp;
    std::cout << std::endl;
  std::cout << "dis:"; 
  for (auto der : op_info.dis) 
    std::cout << " " << der;
    std::cout << std::endl;
  std::cout << "gamma:"; 
  for (auto dir : op_info.gamma) 
    std::cout << " " << dir;
    std::cout << std::endl;
    const size_t T = 96;
    std::vector<std::complex<double> > dummy_corr (T);
  //generate random numbers for test
  float re[96];
  float im[96];
  for (int i = 0; i < 3; i++){
    int seed=1337*1;
    rlxs_init(2,1337);
    ranlxs(re, 96);
    ranlxs(im, 96);

    for (int t = 0; t < T; ++t){
     dummy_corr[t]=( std::complex<double>(re[t], im[t]) );
    }
//  for (auto& el : dummy_corr) std::cout << el << std::endl;
  //read from testfile
    write_2pt_lime("big_test", op_info, op_info, dummy_corr);
  }
  for (auto& el : dummy_corr) std::cout << el << std::endl;
  return 0;
}

