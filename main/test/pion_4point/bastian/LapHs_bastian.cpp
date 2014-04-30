//============================================================================
// Name        : LapHs.cpp
// Author      : BK
// Version     :
// Copyright   : Copies are prohibited so far
// Description : stochastic LapH code
//============================================================================

#include <iostream>
#include <complex>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <typeinfo>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "GlobalData.h"
#include "BasicOperator.h"
#include "ReadWrite.h"

int main (int ac, char* av[]) {

	Eigen::initParallel();

	// Reading in global parameters from input file
	GlobalData* global_data = GlobalData::Instance();
	global_data->read_parameters(ac, av);

	Eigen::setNbThreads(4);

	// global variables from input file needed in main function
	const int Lt = global_data->get_Lt();
	const int end_config = global_data->get_end_config();
	const int delta_config = global_data->get_delta_config();
	const int start_config = global_data->get_start_config();
	const int number_of_max_mom = global_data->get_number_of_max_mom();
	const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();

	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
	//const int number_of_inversions = quarks[0].number_of_dilution_T
  //			* quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D;

	clock_t time;

	const std::complex<double> I(0.0, 1.0);

	char outfile[400];
	FILE *fp = NULL;

	// ***************************************************************************
	// memory allocation *********************************************************
	// memory for the operator in Dirac and eigenvector space and in time
	// ***************************************************************************
  // ***************************************************************************

//  Eigen::MatrixXcd* op_tsource_tsink = new Eigen::MatrixXcd[number_of_rnd_vec];
//  Eigen::MatrixXcd* op_tsink_tsource = new Eigen::MatrixXcd[number_of_rnd_vec];

//  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){
//    (op_tsource_tsink[rnd_i]).resize(
//        4 * number_of_eigen_vec, 4 * number_of_eigen_vec);
//    (op_tsink_tsource[rnd_i]).resize(
//        4 * number_of_eigen_vec, 4 * number_of_eigen_vec);
//  }


	Eigen::MatrixXcd op_D_tsource = Eigen::MatrixXcd::Zero(
			quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D,
			4 * number_of_eigen_vec);
	Eigen::MatrixXcd op_D_tsource_1 = Eigen::MatrixXcd::Zero(
			quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D,
			4 * number_of_eigen_vec);
	Eigen::MatrixXcd op_D_tsink = Eigen::MatrixXcd::Zero(
			quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D,
			4 * number_of_eigen_vec);
	Eigen::MatrixXcd op_D_tsink_1 = Eigen::MatrixXcd::Zero(
			quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D,
			4 * number_of_eigen_vec);
	Eigen::MatrixXcd** X = new Eigen::MatrixXcd*[number_of_rnd_vec];
	Eigen::MatrixXcd** Y = new Eigen::MatrixXcd*[number_of_rnd_vec];
	for(int rnd1 = 0; rnd1 < number_of_rnd_vec; rnd1++){
		X[rnd1] = new Eigen::MatrixXcd[number_of_rnd_vec];
		Y[rnd1] = new Eigen::MatrixXcd[number_of_rnd_vec];
		for(int rnd2 = 0; rnd2 < number_of_rnd_vec; rnd2++){
			X[rnd1][rnd2] = Eigen::MatrixXcd::Zero(
					4 * number_of_eigen_vec, 4 * number_of_eigen_vec);
			Y[rnd1][rnd2] = Eigen::MatrixXcd::Zero(
					4 * number_of_eigen_vec, 4 * number_of_eigen_vec);
		}
	}

	// memory for the correlation function
	std::complex<double>*** C2_mes = new std::complex<double>**[number_of_max_mom];
	for(int p = 0; p < number_of_max_mom; ++p){
		C2_mes[p] = new std::complex<double>*[16]; // 16 Dirac matrices
		for(int dirac = 0; dirac < 16; ++dirac){
			C2_mes[p][dirac] = new std::complex<double>[Lt];
		}
	}

	// everything for operator handling
	BasicOperator* basic = new BasicOperator;

  // everything for Input files
  ReadWrite* rewr = new ReadWrite;

	// intermidiate memory for traces
	double part1[number_of_rnd_vec][number_of_rnd_vec];
	double part2[number_of_rnd_vec][number_of_rnd_vec];

	int norm = 0;
	for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
	  for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){
	    for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){
	      if((rnd2 != rnd1) && (rnd2 != rnd3)){
		      for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
			      if((rnd4 != rnd1) && (rnd4 != rnd3)){
				      norm++;
				      //std::cout << "\n\nnorm: " << norm << rnd1 << rnd3 << rnd2 << rnd4 << std::endl;
			      }
		      }
	      }
      }
    }
  }

	std::cout << "\n\tNumber of contraction combinations: " << norm << std::endl;
	const double norm1 = Lt * norm;

  // Memory for propagation matrices (is that a word?) from t_source to t_sink
  // (op_1) and vice versa (op_2)
  // additional t_source to t_sink (op_3) and t_sink to t_source (op_4) for
  // 4-point functions
  Eigen::MatrixXcd* op_1 = new Eigen::MatrixXcd[number_of_rnd_vec];
  Eigen::MatrixXcd* op_2 = new Eigen::MatrixXcd[number_of_rnd_vec];
  Eigen::MatrixXcd* op_3 = new Eigen::MatrixXcd[number_of_rnd_vec];
  Eigen::MatrixXcd* op_4 = new Eigen::MatrixXcd[number_of_rnd_vec];

  for(int rnd_i = 0; rnd_i < number_of_rnd_vec; ++rnd_i){
    op_1[rnd_i] = Eigen::MatrixXcd(4 * number_of_eigen_vec, 
        4 * number_of_eigen_vec);
    op_2[rnd_i] = Eigen::MatrixXcd(4 * number_of_eigen_vec, 
        4 * number_of_eigen_vec);
    op_3[rnd_i] = Eigen::MatrixXcd(4 * number_of_eigen_vec, 
        4 * number_of_eigen_vec);
    op_4[rnd_i] = Eigen::MatrixXcd(4 * number_of_eigen_vec, 
        4 * number_of_eigen_vec);
  }

	// ***************************************************************************
	// ***************************************************************************
	// Loop over all configurations **********************************************
	// ***************************************************************************
	// ***************************************************************************

  for(int config_i = start_config; config_i <= end_config; config_i +=
			delta_config){

		std::cout << "\nprozessing configuration: " << config_i << "\n\n";

		rewr->read_perambulators_from_file(config_i);
		rewr->read_eigenvectors_from_file(config_i);
		rewr->read_rnd_vectors_from_file(config_i);
		rewr->build_source_matrix();

		// *************************************************************************
		// CONNECTED CONTRACTION 1 *************************************************
		// *************************************************************************

		// setting the correlation function to zero
		std::cout << "\n\tcomputing the connected contribution of pi_+/-:\n";
		time = clock();
		for(int p = 0; p < number_of_max_mom; ++p)
			for(int dirac = 0; dirac < 16; ++dirac)
				for(int t1 = 0; t1 < Lt; ++t1)
					C2_mes[p][dirac][t1] = std::complex<double>(0.0, 0.0);

    int dirac_min = 5;
    int dirac_max = 6;

		for(int t_source = 0; t_source < Lt; ++t_source){
			for(int t_sink = 0; t_sink < Lt; ++t_sink){

        // initialize contraction[rnd_i] = perambulator * basicoperator
        // = D_u^-1
        // choose 'i' for interlace or 'b' for block time dilution scheme
        basic->init_operator(t_source, t_sink, rewr, 'b');

        for(int dirac = dirac_min; dirac < dirac_max; ++dirac){

          // "multiply contraction[rnd_i] with gamma structure"
          // contraction[rnd_i] are the columns of D_u^-1 which get
          // reordered by gamma multiplication. No actuall multiplication
          // is carried out
          basic->get_operator(op_1, dirac);

          // same as get_operator but with gamma_5 trick. D_u^-1 is
          // daggered and multipied with gamma_5 from left and right
          basic->get_operator_g5(op_2, dirac);


          for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
            for(int rnd2 = rnd1 + 1; rnd2 < number_of_rnd_vec; ++rnd2){
              for(int p = 0; p < number_of_max_mom; ++p){
                // building Correlation function get quite intuitive
                // C2 = tr(D_d^-1 Gamma D_u^-1 Gamma)
                // TODO: find signflip of imaginary part
                C2_mes[p][dirac][abs((t_sink - t_source - Lt) % Lt)] += 
                    (op_2[rnd2] * op_1[rnd1]).trace();
              }
            }
          }   
        }

			}
		}

		double norm3 = Lt * number_of_rnd_vec * (number_of_rnd_vec - 1) * 0.5;
		for(int t = 0; t < Lt; ++t)
			for(int p = 0; p < number_of_max_mom; ++p)
				for(int dirac = dirac_min; dirac < dirac_max; ++dirac)
					C2_mes[p][dirac][t] /= norm3;

    // output to binary file
		sprintf(outfile, "./C2_pi+-_conf%04d.dat", config_i);
		if((fp = fopen(outfile, "wb")) == NULL)
			std::cout << "fail to open outputfile" << std::endl;
		for(int p = 0; p < number_of_max_mom; ++p)
			for(int dirac = dirac_min; dirac < dirac_max; ++dirac)
				fwrite((double*) C2_mes[p][dirac], sizeof(double), 2 * Lt, fp);
		fclose(fp);

    // output to terminal
		printf("\n");
		for(int dirac = dirac_min; dirac < dirac_max; ++dirac){
			printf("\tdirac = %02d\n", dirac);
			printf(
					"\t t\tRe(C2_con)\tIm(C2_con)\n\t----------------------------------\n");
			for(int t1 = 0; t1 < Lt; ++t1){
				printf("\t%02d\t%.5e\t%.5e\n", t1, real(C2_mes[0][dirac][t1]),
						imag(C2_mes[0][dirac][t1]));
			}
			printf("\n");
		}
		printf("\n");

    time = clock() - time;
		printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);

		// *************************************************************************
		// FOUR PT CONTRACTION 1 ***************************************************
		// *************************************************************************

		// setting the correlation function to zero
		std::cout << "\n\tcomputing the connected contribution of C4_1:\n";
		time = clock();
		for(int p = 0; p < number_of_max_mom; ++p)
			for(int dirac = 0; dirac < 16; ++dirac)
				for(int t1 = 0; t1 < Lt; ++t1)
					C2_mes[p][dirac][t1] = std::complex<double>(0.0, 0.0);

		for(int t_source = 0; t_source < Lt; ++t_source){
			for(int t_sink = 0; t_sink < Lt; ++t_sink){

				int t_source_1 = (t_source + 1) % Lt;
				int t_sink_1 = (t_sink + 1) % Lt;

        basic->init_operator(t_source_1, t_sink_1, rewr, 'b');

        basic->get_operator(op_1, 5);
        basic->get_operator_g5(op_2, 5);

        basic->init_operator(t_source, t_sink, rewr, 'b');

        basic->get_operator(op_3, 5);
        basic->get_operator_g5(op_4, 5);

        // first trace
        for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
          // first u quark: t_source_1 -> t_sink_1
          for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){
            // first d quark: t_sink_1 -> t_source_1
            part1[rnd1][rnd3] = std::real((op_2[rnd3] * op_1[rnd1]).trace());
          }
        }

        // second trace
        for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){			
          // second u quark: t_source -> t_sink
          for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
            // second d quark: t_sink -> t_source
            part2[rnd2][rnd4] = std::real((op_4[rnd4] * op_3[rnd2]).trace());
          }
        }

        // complete diagramm
        // every quark line must have its own random vec
				for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
				  for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){
				    for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){			
      				if((rnd2 != rnd1) && (rnd2 != rnd3)){
					      for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
						      if((rnd4 != rnd1) && (rnd4 != rnd3)){
							      C2_mes[0][5][abs((t_sink - t_source - Lt) % Lt)] += 
                        part1[rnd1][rnd3] * part2[rnd2][rnd4];
						      }
					      }
              }
            }
          }
        }

      }
    }

    for(int t = 0; t < Lt; ++t)
			C2_mes[0][5][t] /= norm1;

    // output to binary file
		sprintf(outfile, "./C4_1_conf%04d.dat", config_i);
		if((fp = fopen(outfile, "wb")) == NULL)
			std::cout << "fail to open outputfile" << std::endl;
		fwrite((double*) C2_mes[0][5], sizeof(double), 2 * Lt, fp);
		fclose(fp);

    // output to terminal
		printf("\n");
		for(int dirac = 5; dirac < 6; ++dirac){
			printf("\tdirac = %02d\n", dirac);
			printf(
					"\t t\tRe(C2_con)\tIm(C2_con)\n\t----------------------------------\n");
			for(int t1 = 0; t1 < Lt; ++t1){
				printf("\t%02d\t%.5e\t%.5e\n", t1, real(C2_mes[0][dirac][t1]),
						imag(C2_mes[0][dirac][t1]));
			}
			printf("\n");
		}
		printf("\n");

    time = clock() - time;
		printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);


		// *************************************************************************
		// FOUR PT CONTRACTION 2 ***************************************************
		// *************************************************************************

		// setting the correlation function to zero
		std::cout << "\n\tcomputing the connected contribution of C4_2:\n";
		time = clock();
		for(int p = 0; p < number_of_max_mom; ++p)
			for(int dirac = 0; dirac < 16; ++dirac)
				for(int t1 = 0; t1 < Lt; ++t1)
					C2_mes[p][dirac][t1] = std::complex<double>(0.0, 0.0);

		for(int t_source = 0; t_source < Lt; ++t_source){
			for(int t_sink = 0; t_sink < Lt-1; ++t_sink){

				int t_source_1 = (t_source + 1) % Lt;
				int t_sink_1 = (t_sink + 1) % Lt;

        basic->init_operator(t_source_1, t_sink, rewr, 'b');

        basic->get_operator(op_1, 5);
        basic->get_operator_g5(op_2, 5);

        basic->init_operator(t_source, t_sink_1, rewr, 'b');

        basic->get_operator(op_3, 5);
        basic->get_operator_g5(op_4, 5);

        // first trace
        for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
          // first u quark: t_source_1 -> t_sink
          for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){
            // first d quark: t_sink -> t_source_1
            part1[rnd1][rnd3] = std::real((op_2[rnd3] * op_1[rnd1]).trace());
          }
        }

        // second trace
        for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){			
          // second u quark: t_source -> t_sink_1
          for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
            // second d quark: t_sink_1 -> t_source
            part2[rnd2][rnd4] = std::real((op_4[rnd4] * op_3[rnd2]).trace());
          }
        }

        // complete diagramm
        // every quark line must have its own random vec
				for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
				  for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){
				    for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){			
      				if((rnd2 != rnd1) && (rnd2 != rnd3)){
					      for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
						      if((rnd4 != rnd1) && (rnd4 != rnd3)){
							      C2_mes[0][5][abs((t_sink - t_source - Lt) % Lt)] += 
                        part1[rnd1][rnd3] * part2[rnd2][rnd4];
						      }
					      }
              }
            }
          }
        }

      }
    }

		for(int t = 0; t < Lt; ++t)
			C2_mes[0][5][t] /= norm1;

    // output to binary file
		sprintf(outfile, "./C4_2_conf%04d.dat", config_i);
		if((fp = fopen(outfile, "wb")) == NULL)
			std::cout << "fail to open outputfile" << std::endl;
		fwrite((double*) C2_mes[0][5], sizeof(double), 2 * Lt, fp);
		fclose(fp);

    // output to terminal
		printf("\n");
		for(int dirac = 5; dirac < 6; ++dirac){
			printf("\tdirac = %02d\n", dirac);
			printf(
					"\t t\tRe(C2_con)\tIm(C2_con)\n\t----------------------------------\n");
			for(int t1 = 0; t1 < Lt; ++t1){
				printf("\t%02d\t%.5e\t%.5e\n", t1, real(C2_mes[0][dirac][t1]),
						imag(C2_mes[0][dirac][t1]));
			}
			printf("\n");
		}
		printf("\n");

    time = clock() - time;
		printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);

#endif

		// *************************************************************************
		// FOUR PT CONTRACTION 3 ***************************************************
		// *************************************************************************

		// setting the correlation function to zero
		std::cout << "\n\tcomputing the connected contribution of C4_3:\n";
		time = clock();
		for(int p = 0; p < number_of_max_mom; ++p)
			for(int dirac = 0; dirac < 16; ++dirac)
				for(int t1 = 0; t1 < Lt; ++t1)
					C2_mes[p][dirac][t1] = std::complex<double>(0.0, 0.0);

		for(int t_source = 0; t_source < Lt; ++t_source){
			for(int t_sink = 0; t_sink < Lt; ++t_sink){

				int t_source_1 = (t_source + 1) % Lt;
				int t_sink_1 = (t_sink + 1) % Lt;

#if 0
        basic->init_operator(t_source, t_sink, rewr, 'b');
        basic->get_operator(op_1, 5);

        basic->init_operator(t_source, t_sink_1, rewr, 'b');
        basic->get_operator_g5(op_2, 5);

        basic->init_operator(t_source_1, t_sink_1, rewr, 'b');
        basic->get_operator(op_3, 5);

        basic->init_operator(t_source_1, t_sink, rewr, 'b');
        basic->get_operator_g5(op_4, 5);

        // first part
        for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
          // first u quark: t_source -> t_sink
          for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){
            // first d quark: t_sink_1 -> t_source
            X[rnd1][rnd3] = op_1[rnd3] * op_2[rnd1];
          }
        }

        // second part
        for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){			
          // second u quark: t_source_1 -> t_sink_1
          for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
            // second d quark: t_sink -> t_source_1
            Y[rnd2][rnd4] = op_3[rnd4] * op_4[rnd2];
          }
        }

        // complete diagramm
        // every quark line must have its own random vec
				for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
				  for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){
				    for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){			
      				if((rnd2 != rnd1) && (rnd2 != rnd3)){
					      for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
						      if((rnd4 != rnd1) && (rnd4 != rnd3)){
                    // -= accounts for - sign from wick contraction
							      C2_mes[0][5][abs((t_sink - t_source - Lt) % Lt)] += 
                        (X[rnd1][rnd3] * Y[rnd2][rnd4]).trace();
						      }
					      }
              }
            }
          }
        }
#endif
        //test other indices
        basic->init_operator(t_source_1, t_sink_1, rewr, 'b');
        basic->get_operator(op_1, 5);

        basic->init_operator(t_source, t_sink_1, rewr, 'b');
        basic->get_operator_g5(op_2, 5);

        basic->init_operator(t_source, t_sink, rewr, 'b');
        basic->get_operator(op_3, 5);

        basic->init_operator(t_source_1, t_sink, rewr, 'b');
        basic->get_operator_g5(op_4, 5);

        // first part
        for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
          // first u quark: t_source -> t_sink
          for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){
            // first d quark: t_sink_1 -> t_source
            X[rnd1][rnd3] = op_1[rnd1] * op_2[rnd3];
          }
        }

        // second part
        for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){			
          // second u quark: t_source_1 -> t_sink_1
          for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
            // second d quark: t_sink -> t_source_1
            Y[rnd2][rnd4] = op_3[rnd2] * op_4[rnd4];
          }
        }

        // complete diagramm
        // every quark line must have its own random vec
				for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
				  for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){
				    for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){			
      				if((rnd2 != rnd1) && (rnd2 != rnd3)){
					      for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
						      if((rnd4 != rnd1) && (rnd4 != rnd3)){
                    // -= accounts for - sign from wick contraction
							      C2_mes[0][5][abs((t_sink - t_source - Lt) % Lt)] += 
                        (X[rnd1][rnd3] * Y[rnd2][rnd4]).trace();
						      }
					      }
              }
            }
          }
        }

      }
    }


#if 0

				for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
					// first u quark: t_source_1 -> t_sink_1
					peram_1 = (basic->perambulator[rnd1]).block(
							4 * number_of_eigen_vec * t_source_1,
							quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D
									* ( (t_sink_1/2 + t_sink_1%2) % (quarks[0].number_of_dilution_T)),  // changed for block interlace dilution
							4 * number_of_eigen_vec,
							quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D);
					op_D_tsink_1 = peram_1
							* basic->basicoperator[rnd1][t_sink_1][5];
					for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){
						// first d quark: t_sink_1 -> t_source_1
						peram_3 = (basic->perambulator[rnd3]).block(
								4 * number_of_eigen_vec * t_source,
								quarks[0].number_of_dilution_D * quarks[0].number_of_dilution_E
										* ( (t_sink_1/2 + t_sink_1%2) % quarks[0].number_of_dilution_T),		// changed for block interlace dilution
								4 * number_of_eigen_vec,
								quarks[0].number_of_dilution_D
										* quarks[0].number_of_dilution_E);
						op_D_tsource_1 =
								(peram_3 * basic->basicoperator[rnd3][t_sink_1][5]).adjoint();
						op_D_tsource_1 = basic->mul_r_gamma(op_D_tsource_1,
								op_D_tsource_1.rows(), op_D_tsource_1.cols(), 5, 3);

						X[rnd1][rnd3] = op_D_tsink_1 * op_D_tsource_1;

					}
				}
				for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){
					// second u quark: t_source -> t_sink
					peram_2 = (basic->perambulator[rnd2]).block(
							4 * number_of_eigen_vec * t_source,
							quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D
									* ( (t_sink/2 + t_sink%2) % quarks[0].number_of_dilution_T),    // changed for block interlace dilution
							4 * number_of_eigen_vec,
							quarks[0].number_of_dilution_E
									* quarks[0].number_of_dilution_D);
					op_D_tsink = peram_2
							* basic->basicoperator[rnd2][t_sink][5];
					for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
						// second d quark: t_sink -> t_source
						peram_4 = (basic->perambulator[rnd4]).block(
								4 * number_of_eigen_vec * t_source_1,
								quarks[0].number_of_dilution_D * quarks[0].number_of_dilution_E
										* ( (t_sink/2 + t_sink%2) % quarks[0].number_of_dilution_T),		// changed for block interlace dilution
								4 * number_of_eigen_vec,
								quarks[0].number_of_dilution_D
										* quarks[0].number_of_dilution_E);
						op_D_tsource =
								(peram_4 * basic->basicoperator[rnd4][t_sink][5]).adjoint();
						op_D_tsource = basic->mul_r_gamma(op_D_tsource,
								op_D_tsource.rows(), op_D_tsource.cols(), 5, 3);

						Y[rnd2][rnd4] = op_D_tsink * op_D_tsource;

					}
				}


				// building correlation function
				for(int rnd1 = 0; rnd1 < number_of_rnd_vec; ++rnd1){
				for(int rnd3 = rnd1 + 1; rnd3 < number_of_rnd_vec; ++rnd3){
				for(int rnd2 = 0; rnd2 < number_of_rnd_vec; ++rnd2){			
				if((rnd2 != rnd1) && (rnd2 != rnd3)){
					for(int rnd4 = rnd2 + 1; rnd4 < number_of_rnd_vec; ++rnd4){
						if((rnd4 != rnd1) && (rnd4 != rnd3)){
							C2_mes[0][5][abs((t_sink - t_source - Lt) % Lt)]  += (X[rnd1][rnd3] * Y[rnd2][rnd4]).trace();
						}
					}
				}}}}


			}
		}

		for(int t = 0; t < Lt; ++t)
			C2_mes[0][5][t] /= norm1;

		// output to binary file
		sprintf(outfile, "./C4_3_conf%04d.dat", config_i);
		if((fp = fopen(outfile, "wb")) == NULL)
			std::cout << "fail to open outputfile" << std::endl;
		fwrite((double*) C2_mes[0][5], sizeof(double), 2 * Lt, fp);
		fclose(fp);

    // output to terminal
		printf("\n");
		for(int dirac = 5; dirac < 6; ++dirac){
			printf("\tdirac = %02d\n", dirac);
			printf(
					"\t t\tRe(C2_con)\tIm(C2_con)\n\t----------------------------------\n");
			for(int t1 = 0; t1 < Lt; ++t1){
				printf("\t%02d\t%.5e\t%.5e\n", t1, real(C2_mes[0][dirac][t1]),
						imag(C2_mes[0][dirac][t1]));
			}
			printf("\n");
		}
		printf("\n");

    time = clock() - time;
		printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);
		
		// *************************************************************************
		// FOUR PT CONTRACTION 4 ***************************************************
		// *************************************************************************

    // identical to FOUR PT CONTRACTION 3

	} // loop over configs ends here

	// TODO: freeing all memory!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	//delete basic;
	return 0;
}

