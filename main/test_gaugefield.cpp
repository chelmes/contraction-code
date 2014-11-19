/***************************************************************************/
//
// Short demonstration of gaugefield class
//
// Author: Christopher Helmes
//
/***************************************************************************/
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <typeinfo>

#include "GaugeField.h"
#include "GlobalData.h"
#include "BasicOperator.h"
#include "typedefs.h"

int main(int ac, char* av[]) {
  //initialize parameters, gaugeconfiguration and eigenvectors
  GlobalData* global_data = GlobalData::Instance();
  global_data->read_parameters(ac, av);
  GaugeField gf(0,3,64,3);
  LapH::EigenVector evec(4,4*4*4*3, 4);
  //read in eigenvectors
  evec.read_eigen_vector("/hiskp2/eigensystems/test4x4x4x4/eigenvectors.1000.",1);
  //initialize lookup table
  std::cout << "\ninitialize lookup tables:" << std::endl;
  gf.init(4,4,4);
  //read in gaugefield
  std::cout << "\nRead in complete 4^4 configuration:" << std::endl;
  gf.read_gauge_field(0,3);
  std::cout << "\n U at t = 3, vol = 4, dir = z:" << std::endl;
  std::cout << gf(3,4,2) << std::endl;
  //smear gauge field
  //gf.smearing_stout(3,0.62,3);
  std::cout << "\nPlaquette value at timeslice 0: " << 
    gf.plaque_ts(0) << std::endl;
  // trafo all timeslices
  std::cout << "\nGauge transformation of whole confiugration" << std::endl;
  gf.trafo(0,3);
  std::cout << "\nPlaquette value at timeslice 0 after gauge_config:" <<
    gf.plaque_ts(0) << std::endl;
  Eigen::MatrixXcd evec_t = gf.trafo_ev(evec[0]);
  std::cout << (evec[0] * evec[0].adjoint()).trace() << std::endl;
  //displace eigensystem
  Eigen::MatrixXcd W = Eigen::MatrixXcd::Zero(4*4*4*3,4);
    W.col(1) = gf.disp(evec[0].col(1),0,0,true);
  //std::cout << (W * W.adjoint()).trace() << std::endl;
  //std::cout << gf(3,4,2) << std::endl;
  return 0;

}
