/*! \file Gaugefield.h
    \brief Contains the class for gaugefield handling

    The Gauge fields of a range of timeslices are read into a vector of
    Multiarrays. Functions: - 
 */

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

//! \typedef 2dim array for lookup tables
typedef boost::multi_array<int, 2> look;
//! \typedef 2dim array for gauge field matrices of one timeslice.
//! 
typedef boost::multi_array<Eigen::Matrix3cd, 2> array_3cd_d2_eigen;
class GaugeField {
  public:
    //! \brief Constructor

    //!t0 is first, tf is last timeslice to be read in, v3 denotes the spatial
    //!lattice sites LX * LY * LZ
    //!
    GaugeField(const size_t t0, const size_t tf, const size_t v3,
               const size_t ndir);
    ~GaugeField() {};
    //! \brief overloaded operator() returning gauge field matrix
    //!
    //!
    inline const Eigen::Matrix3cd& operator()(const size_t t, const size_t v,
                                              const size_t dir) const {
        return (tslices.at(t))[v][dir];
    };

    void read_gauge_field(const size_t slice_i, const size_t slice_f);
    void map_timeslice_to_eigen(const size_t t, const double* timeslice);
    void smearing_hyp(const size_t t, const double alpha_1,
                      const double alpha_2, const size_t iter);
    void smearing_stout(const size_t t, const double rho, const size_t iter);
    void smearing_ape(const size_t t, const double alpha_1, const size_t iter);
    void init(const size_t LX, const size_t LY, const size_t LZ);
    //Test functions for navigation
    int get_up(const int pos, const int dir);
    int get_dn(const int pos, const int dir);
    //Displacement as Vector and Matrix form
    Eigen::MatrixXcd disp(const Eigen::MatrixXcd& v, const size_t t,
                             const size_t dir, bool sym);
  private:
  
  void read_lime_gauge_field_doubleprec_timeslices(double* gaugefield,
                                                   const char* filename,
                                                   const size_t slice_i,
                                                   const size_t slice_f);
  //vector holding multiarrays for range of timeslices
  std::vector<array_3cd_d2_eigen> tslices;
  //One gaugefield of size V3
  Eigen::MatrixXcd omega;
  //!
  //!Memory for the whole gauge configuration in lime layout
  //!
  look iup;
  look idown;
};

#endif //GAUGEFIELD_H_
