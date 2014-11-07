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
  std::cout << gf(0,0,0) << std::endl;
  gf.read_gauge_field(0,3);
  std::cout << gf(0,0,0) << std::endl;
  return 0;

}
