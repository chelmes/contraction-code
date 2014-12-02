#ifndef CORRELATOR_IO_2PT_H_
#define CORRELATOR_IO_2PT_H_

#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <complex>
#include <iostream>
#include <vector>
//TODO: After pull request stick together (collect structs in single header?)
//#include "OperatorStructure.h"
#include "lime.h"

struct pdg{
  size_t id;
  std::array<int,3> p;
  std::array<int,3> dis;
  std::array<int,4> gamma;
};

struct Tag{
  int mom[2];
  int dis[2];
  int gam[2];
};
void write_2pt_lime(const char* filename, const Tag& tag
                    std::vector<std::complex<double> >& corr);

void read_2pt_lime(const char* filename, const Tag& tag,
                   std::vector<std::complex<double> >& corr);


#endif //IO_2PT_H_
