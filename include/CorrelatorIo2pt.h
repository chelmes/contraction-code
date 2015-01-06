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
#include "IoHelpers.h"
#include "boost/crc.hpp"
#include "boost/integer.hpp"
#include "lime.h"
#include "typedefs.h"

// Writes a vector of 2pt correlation functions and their tags to a file with
// filename
void write_2pt_lime(const char* filename, GlobalDat& dat, std::vector<Tag>& tags,
                    std::vector<vec>& corr);

// Reads in a file of correlation functions into a vector of correlation
// functions including their tags. Global Checksum is checked implicitly
void read_2pt_lime(const char* filename, std::vector<Tag>& tag,
                   std::vector<vec>& corr);

// Gets a single correlation function from a file utilizing read_2pt_lime
void get_2pt_lime(const char* filename, const size_t num_corrs,
                  const size_t corr_length, const Tag& tag,
                  std::vector<cmplx >& corr);

// Dump ASCII version of correlator on screen
void ASCII_dump_2pt(const char* filename, size_t g_so, size_t g_si, size_t p_so,
                  size_t p_si, size_t dis_so, size_t dis_si );
#endif //IO_2PT_H_
