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

void write_2pt_lime(const char* filename, GlobalDat& dat, std::vector<Tag>& tags,
                    std::vector<vec>& corr);

void read_2pt_lime(const char* filename, const Tag& tag,
                   std::vector<std::complex<double> >& corr);


#endif //IO_2PT_H_
