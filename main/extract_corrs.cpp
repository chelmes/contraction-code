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
#include "TagHandling.h"
#include "typedefs.h"


int main(int ac, char* av[]) {

if (ac < 4){
  std::cout << "Not enough options given!" << std::endl;
  std::cout << "Usage: extract_corrs <CORR-FILE> <NUM-CORRS> <LT>" << std::endl;
  std::cout << "CORR-FILE: Path to lime file containing correlation functions\n"
            << "NUM-CORRS: Number of correlators in file\n"
            << "LT: Time extent of lattice" << std::endl;
  return 1;
}

const size_t Lt = atoi(av[3]);
ASCII_dump_corr(av[1], "tag_pars", Lt, atoi(av[2])); 
 

  return 0;
}

