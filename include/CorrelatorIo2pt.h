#ifndef CORRELATOR_IO_2PT_H_
#define CORRELATOR_IO_2PT_H_

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <complex>
#include <iostream>
#include <vector>
//TODO: After pull request stick together (collect structs in single header?)
//#include "OperatorStructure.h"
#include "boost/crc.hpp"
#include "boost/integer.hpp"
#include "lime.h"
#include "typedefs.h"

struct pdg{
  size_t id;
  std::array<int,3> p;
  std::array<int,3> dis;
  std::array<int,4> gamma;
};

struct Tag {
  int mom[2];
  int dis[2];
  int gam[2];
};

struct GlobalDat {
  std::vector<size_t> rnd_seeds;
  size_t nb_rnd_vecs;
  size_t nb_perambs;
};

void swap_correlators(std::vector<vec>& correlators);
void swap_tag_vector(std::vector<Tag>& tags);
void write_2pt_lime(const char* filename, GlobalDat& dat, std::vector<Tag>& tags,
                    std::vector<vec>& corr);

void read_2pt_lime(const char* filename, const Tag& tag,
                   std::vector<std::complex<double> >& corr);


#endif //IO_2PT_H_
