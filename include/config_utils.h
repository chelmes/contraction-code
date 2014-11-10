/*
 * configs.h
 * Declares: - lookup-tables in 3d (hopping3d)
 *					 - mapping from ildg to Eigen Array (map_timeslice_to_eigen)
 *					 - building of Laplacian in E4- and C-space (BuildLaplacian)
 *					 - gauge transform (transform_ts)
 *					 - check for gauge invariance (check_gauge)
 * Created on: Aug 26, 2013
 * Author: christopher helmes
 */

#ifndef CONFIGS_H_
#define CONFIGS_H_

#include <complex>
#include <iomanip>
#include <iostream>
#include <Eigen/Eigen>
#include <unsupported/Eigen/MatrixFunctions>
#include <array>
#include "GlobalData.h"

// navigation through timeslice
void hopping3d(int** iup, int** idown);

// write timeslice of ildg-message to array of eigen-matrices
//void map_timeslice_to_eigen(Eigen::Matrix3cd **eigen, double *timeslice);

// Right displacement of one eigensystem
void right_displacement_one_dir(Eigen::Matrix3cd** config, int** iup, 
    int** idown, const int dir, Eigen::MatrixXcd& V, Eigen::MatrixXcd& W );

// Hyp-Smearing of one timeslice eigeen timeslice is overwritten smearing takes 
// place in 3 dimensions using alpha_1, alpha_2 as staple weights and iter as 
// iteration number
void smearing_hyp(Eigen::Matrix3cd **eigen_timeslice, double alpha_1, double alpha_2, int iter);

// empty function header for compiling! No single precision support
//void read_lime_gauge_field_singleprec(double *config, const char * filename,
//    const int T, const int LX, const int LY, const int LZ);

#endif /* CONFIGS_H_ */
