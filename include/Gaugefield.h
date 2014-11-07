#ifndef GAUGEFIELD_H_
#define GAUGEFIELD_H_

#include <cmath>
#include <complex>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>

//#include "Eigen/Dense"
//#include "Eigen/Core"
//#include "Eigen/SparseCore"
#include "boost/multi_array.hpp" 

//#include "lime.h"
#include "config_utils.h"
#include "EigenVector.h"
#include "GlobalData.h"
#include "io_utils.h"
#include "propagator_io.h"
#include "quark.h"
#include "RandomVector.h"
#include "typedefs.h"

//2dim array for lookup tables
typedef boost::multi_array<int, 2> look;
typedef boost::multi_array<Eigen::Matrix3cd, 2> array_3cd_d2_eigen;
class GaugeField {
  public:
    //!
    //!standard constructor and destructor,
    //!t0 is first, tf is last timeslice to be read in, v3 denotes the spatial
    //!lattice sites LX * LY * LZ
    //!
    GaugeField(const size_t t0, const size_t tf, const size_t v3,
               const size_t ndir);
    ~GaugeField() {};
    //!
    //![]-Operator gets overloaded to return one SU3-gauge matrix
    //!
    inline const Eigen::Matrix3cd& operator()(const size_t t, const size_t v,
                                              const size_t dir) const {
        return (tslices.at(t))[v][dir];
    };

    void read_gauge_field(const size_t slice_i, const size_t slice_f);
    void map_timeslice_to_eigen(const size_t t, const double* timeslice);
    //void smearing_hyp(const double a_1, const double a_2, const int it);
    //void smearing_stout(double rho, int iter);
    //void smearing_ape( double alpha, int iter);
  private:
  
void read_lime_gauge_field_doubleprec_timeslices(double* gaugefield, const char* filename, const size_t slice_i, const size_t slice_f);
  //vector holding multiarrays for range of timeslices
  std::vector<array_3cd_d2_eigen> tslices;
  //!
  //!Memory for the whole gauge configuration in lime layout
  //!
  look iup;
  look idown;
};

#endif //GAUGEFIELD_H_
