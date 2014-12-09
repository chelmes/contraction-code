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

#include "boost/multi_array.hpp"
#include "boost/crc.hpp"
#include "CorrelatorIo2pt.h"
//#include "BasicOperator.h"
//#include "GlobalData.h"
#include "ranlxs.h"
#include "typedefs.h"

void fill_corr_rand(std::vector<cmplx>& one_corr, const int mult){
  const int T = one_corr.size();
  float re[T];
  float im[T];
    ranlxs(re, T);
    ranlxs(im, T);

    for (int t = 0; t < T; ++t){
     one_corr[t]=( std::complex<double>(re[t], im[t]) );
    }
}

int main(int ac, char* av[]) {
  
  // Global Message for the computation containing number of random vectors,
  // parameters and global checksum
  GlobalDat run_id;
  rlxs_init(2,1337);
  // 100 Correlators
  std::vector<vec> correlators(100);
  for (auto& el : correlators) el.resize(96);
  // attributes of the correlators
  std::vector<Tag> attributes(100);

  // Fill correlators with random numbers
  for (auto& el : correlators) fill_corr_rand(el, &el-&correlators[0]);
  write_2pt_lime("final_write", run_id, attributes, correlators);
  //swap_correlators(correlators);
  //swap_correlators(correlators);
  //for (auto& el : correlators.at(0)) std::cout << el << std::endl;
  //concatenate all correlation functions in one vector
  std::vector<cmplx> collect(correlators.size()*correlators[0].size());
  size_t length = correlators.size()*correlators[0].size()*2*sizeof(double);
  for(auto& c : correlators)
    for (auto& el : c) collect.push_back(el);
  boost::crc_32_type chk_agent;
  size_t bytes = correlators.size()*(correlators[0]).size()*2*sizeof(double);
  std::cout << bytes << std::endl;
  chk_agent.process_bytes(collect.data(), bytes); 
  std::cout << "Checksum for Correlators is:" << chk_agent.checksum() << std::endl;

  std::vector<Tag> tags_in(100);
  std::vector<vec> correlators_in(100);
  std::cout << "read_in from file: final_write " << std::endl;
  read_2pt_lime("final_write", tags_in, correlators_in);
//  for (auto& el : correlators_in.at(0)) std::cout << el << std::endl;
 

  return 0;
}

