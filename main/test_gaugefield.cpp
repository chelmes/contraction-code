#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <typeinfo>

#include "Gaugefield.h"
#include "GlobalData.h"
#include "BasicOperator.h"
#include "ReadWrite.h"
#include "typedefs.h"

int main(int ac, char* av[]) {

  GlobalData* global_data = GlobalData::Instance();
  global_data->read_parameters(ac, av);
  GaugeField gf(0,3,64,3);
  gf.init(4,4,4);
  gf.read_gauge_field(0,3);
  std::cout << gf(3,4,2) << std::endl;
  gf.smearing_hyp(3,0.62,0.62,2);
  std::cout << gf(3,4,2) << std::endl;
  std::cout << gf.get_up(4,1) << std::endl;
  return 0;

}
