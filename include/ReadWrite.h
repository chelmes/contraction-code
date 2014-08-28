/*
 * ReadWrite.h
 *
 *  Created on: Apr 13, 2014
 *      Author: werner
 */

#ifndef _READ_WRITE_H_
#define _READ_WRITE_H_

// TODO: check if they are all necessary. Doesn't matter much though
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <typeinfo>
#include <vector>

//#include "Eigen/Dense"
//#include "Eigen/Core"
//#include "Eigen/SparseCore"
//#include "boost/multi_array.hpp" 

#include "lime.h"
#include "GlobalData.h"
#include "propagator_io.h"
#include "quark.h"
#include "io_utils.h"
#include "config_utils.h"
#include "typedefs.h"

/***************************Input from files**********************************/


//typedef boost::multi_array<Eigen::MatrixXcd, 3> dim3_eigen_array;
//typedef std::vector<Eigen::MatrixXcd> vec_eigen;
//typedef std::complex<double> cmplx;
//typedef std::vector<cmplx> vec;
//typedef boost::multi_array<cmplx, 2> dim2_array;

class ReadWrite {

public:
  //EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  ReadWrite ();
//  virtual ~ReadWrite (); 
  void build_source_matrix (const int config_i);
  void read_perambulators_from_file (const int config_i);
  void read_rnd_vectors_from_file (const int config_i);
  void read_lime_gauge_field_doubleprec_timeslices(const int config_i);

  inline vec_Xcd_eigen get_perambulator() {
    return perambulator;
  }
  inline array_Xcd_d3_eigen get_basicoperator() {
    return basicoperator;
  }
  inline vec_Xcd_eigen get_random_vector() {
    return rnd_vec;
  }
//  int* mom_squared;
//  int number_of_momenta;

protected:
  vec_Xcd_eigen perambulator;
  vec_Xcd_eigen rnd_vec;
  array_Xcd_d3_eigen basicoperator;
  array_cd_d2 momentum;
  Eigen::Matrix3cd** eigen_timeslice;
  double* gaugefield;
  int** iup;
  int** idown;
  void read_eigenvectors_from_file (Eigen::MatrixXcd& V, const int config_i, 
      const int t);
};

#endif // _READ_WRITE_H__
