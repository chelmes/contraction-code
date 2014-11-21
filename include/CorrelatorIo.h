#ifndef CORRELATOR_IO_
#define CORRELATOR_IO_

#include <fstream>
#include <utility>
#include <cstdio>

namespace LapH {

  typedef std::vector<std::pair> corr_info;

  void write_2pt_lime(const char* name, const corr_info* meta, const array_cd_d6* C2_mes);
  void write_4pt_lime(const char* name, const corr_info* meta, const array_cd_d5* C4_mes);



}// end of namespace



#endif //CORRELATOR_IO_
