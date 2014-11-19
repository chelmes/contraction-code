#include "VdaggerV.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
LapH::VdaggerV::VdaggerV() : vdaggerv(), rvdaggervr(), momentum(), nb_mom(1),
                             is_vdaggerv_set(false) {

  // is needed in the whole class
  nb_mom = global_data->get_number_of_momenta();

  const size_t nb_ev = global_data->get_number_of_eigen_vec();
  const size_t Lt = global_data->get_Lt();
  const size_t Vs = global_data->get_Lx() * global_data->get_Ly() * 
                 global_data->get_Lz();         
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t dilE = quarks[0].number_of_dilution_E;

  // TODO: just a workaround, must be read in via infile
  const size_t nb_dis = 1;

  // only half of the array is stored to save memory. But be carefull, it 
  // must be mapped correctly from outside by addressing the memomentum
  // correctly and daggering
  vdaggerv.resize(boost::extents[nb_mom/2+1][Lt][4]);

  //rvdaggervr.resize(boost::extents[nb_mom][Lt][4][nb_rnd][nb_rnd]);
  rvdaggervr.resize(boost::extents[nb_mom/2+1][Lt][nb_dis][nb_rnd][nb_rnd]);

  momentum.resize(boost::extents[nb_mom][Vs]);
  create_momenta();

}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::VdaggerV::create_momenta () {

  const int Lx = global_data->get_Lx();
  const int Ly = global_data->get_Ly();
  const int Lz = global_data->get_Lz();

  static const std::complex<double> I(0.0, 1.0);

  //const int number_of_max_mom = global_data->get_number_of_max_mom();
  const int max_mom_in_one_dir = global_data->get_max_mom_in_one_dir();

  // helper variables for momenta
  const double px = 2. * M_PI / (double) Lx;
  const double py = 2. * M_PI / (double) Ly;
  const double pz = 2. * M_PI / (double) Lz;
  int p = 0;
  int max_mom_squared = global_data->get_number_of_max_mom();

  // running over all momentum components
  for(int ipx = -max_mom_in_one_dir; ipx <= max_mom_in_one_dir; ++ipx){
    for(int ipy = -max_mom_in_one_dir; ipy <= max_mom_in_one_dir; ++ipy){
      for(int ipz = -max_mom_in_one_dir; ipz <= max_mom_in_one_dir; ++ipz){
        if((ipx * ipx + ipy * ipy + ipz * ipz) > max_mom_squared) {
          continue;
        }
        //TODO: for Lx == Ly == Lz ipxH and ipxHipyH may be integers and px, 
        //py get multiplied in the exponential
        // running over all lattice points
        for(int x = 0; x < Lx; ++x){
          const int xH = x * Ly * Lz; // helper variable
          const double ipxH = ipx * px * x; // helper variable
          for(int y = 0; y < Ly; ++y){
            const int xHyH = xH + y * Lz; // helper variable
            const double ipxHipyH = ipxH + ipy * py * y; // helper variable
            for(int z = 0; z < Lz; ++z){
              momentum[p][xHyH + z] = exp(-I * (ipxHipyH + ipz * pz * z));
            }
          }
        }
        ++p;
      }
    }
  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
static void read_eigenvectors_from_file (LapH::EigenVector& V,
                                         const int config_i, const int t) {
    char name[200];
    std::string filename = global_data->get_path_eigenvectors() + "/"
        + global_data->get_name_eigenvectors();
    sprintf(name, "%s.%04d.%03d", filename.c_str(), config_i, t);

    V.read_eigen_vector(name, 0, 0);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::VdaggerV::build_vdaggerv (const int config_i) {

  clock_t t2 = clock();
  std::cout << "\tbuild vdaggerv:";

  const size_t Lt = global_data->get_Lt();
  const size_t dim_row = global_data->get_dim_row();
  const size_t displ_min = global_data->get_displ_min();
  const size_t displ_max = global_data->get_displ_max();
  const size_t nb_ev = global_data->get_number_of_eigen_vec();

  std::fill(vdaggerv.origin(), vdaggerv.origin() + vdaggerv.num_elements(), 
            Eigen::MatrixXcd::Zero(nb_ev, nb_ev));

#pragma omp parallel
{
  Eigen::MatrixXcd W_t = Eigen::MatrixXcd::Zero(dim_row, nb_ev);
  Eigen::VectorXcd mom = Eigen::VectorXcd::Zero(dim_row);
  LapH::EigenVector V_t(1, dim_row, nb_ev);// each thread needs its own copy
  #pragma omp for schedule(dynamic)
  for(size_t t = 0; t < Lt; ++t){

    read_eigenvectors_from_file(V_t, config_i, t);
    // zero momentum 
    (vdaggerv[nb_mom/2][t][0]) = Eigen::MatrixXcd::Identity(nb_ev, nb_ev);

    for(size_t p = 0; p < nb_mom/2; p++){
      // momentum vector contains exp(-i p x). Divisor 3 for colour index. 
      // All three colours on same lattice site get the same momentum.
      for(size_t x = 0; x < dim_row; ++x) {
        mom(x) = momentum[p][x/3];
      }
      vdaggerv[p][t][0] = V_t[0].adjoint() * mom.asDiagonal() * V_t[0];
    } // loop over momentum
  } // loop over time
}// pragma omp parallel ends here

  // set flag that vdaggerv is set
  is_vdaggerv_set = true;

  t2 = clock() - t2;
  std::cout << std::setprecision(1) << "\t\t\tSUCCESS - " << std::fixed 
    << ((float) t2)/CLOCKS_PER_SEC << " seconds" << std::endl;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::VdaggerV::build_rvdaggervr(const int config_i,
                             const std::vector<LapH::RandomVector>& rnd_vec) {
  // check of vdaggerv is already build
  if(not is_vdaggerv_set){
    std::cout << "\n\n\tCaution: vdaggerv is not set and rvdaggervr cannot be" 
              << " computed\n\n" << std::endl;
    exit(0);
  }

  clock_t t2 = clock();
  std::cout << "\tbuild rvdaggervr:";

  const size_t Lt = global_data->get_Lt();
  const size_t dim_row = global_data->get_dim_row();
  const size_t displ_min = global_data->get_displ_min();
  const size_t displ_max = global_data->get_displ_max();
  const size_t nb_ev = global_data->get_number_of_eigen_vec();
  const size_t nb_mom = global_data->get_number_of_momenta();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t dilE = quarks[0].number_of_dilution_E;
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;

  std::fill(rvdaggervr.origin(), 
            rvdaggervr.origin() + rvdaggervr.num_elements(), 
            Eigen::MatrixXcd::Zero(dilE, 4*dilE));

  // TODO: just a workaround
  const size_t nb_dis = 1;

  #pragma omp parallel for schedule(dynamic)
  for(size_t t = 0; t < Lt; t++){
  for(size_t p = 0; p <= nb_mom/2; p++){
    for(size_t dis = 0; dis < nb_dis; ++dis) {
      for(size_t rnd_i = 0; rnd_i < nb_rnd; ++rnd_i) {

        Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(nb_ev, 4*dilE);
        for(size_t dirac = 0; dirac < 4; dirac++){
          for(size_t vec_i = 0; vec_i < nb_ev; ++vec_i) {
            size_t blk_i =  dirac + vec_i * 4 + 4 * nb_ev * t;
          
              M.block(0, vec_i%dilE + dilE*dirac, nb_ev, 1) += 
                                vdaggerv[p][t][dis].col(vec_i) * 
                      rnd_vec[rnd_i][blk_i];
          }
        }
        for(size_t rnd_j = 0; rnd_j < nb_rnd; ++rnd_j){
          if(rnd_i != rnd_j){
            for(size_t dirac = 0; dirac < 4; dirac++){
              for(size_t vec_j = 0; vec_j < nb_ev; ++vec_j) {
                size_t blk_j =  dirac + vec_j * 4 + 4 * nb_ev * t;

                rvdaggervr[p][t][dis][rnd_j][rnd_i]
                              .block(vec_j%dilE, dilE*dirac , 1, dilE) +=
                                     M.block(vec_j, dilE*dirac, 1, dilE) * 
                                             std::conj(rnd_vec[rnd_j][blk_j]);
              }
            }
            vdaggerv[p][t][0] = V_t[0].adjoint() * mom.asDiagonal() * V_t[0];
            vdaggerv[nb_mom - p - 1][t][0] =
                (vdaggerv[p][t][0]).adjoint();
             
           } // end if momentum
         // case displacement
         } 
        else {

          (W_t).setZero();

          //Time Slice of Configuration
          //Factor 3 for second color index, 2 for complex numbers
          double* timeslice = gaugefield + (t * dim_row * 3 * 2);

          //Write Timeslice in Eigen Array
          map_timeslice_to_eigen(eigen_timeslice, timeslice);

          //displacement in one direction i acting to the right
//          right_displacement_one_dir(eigen_timeslice, iup, idown, dir - 1, 
//              V_t[0], W_t);
          // dir = 3 for z-displacement (1 x 2 y)
          // W holds DV

          if(p == nb_mom/2) {
            for(size_t ev = 0; ev < nb_ev; ev++){
              W_t.col(ev) = field -> lr_disp(this, p, t, dir, ev);
            }
            vdaggerv[nb_mom/2][t][dir] = V_t[0].adjoint() * W_t
                - W_t.adjoint() * V_t[0];

            std::cout << vdaggerv[nb_mom/2][t][dir].block(0,0,3,3) << std::endl; 
          } 
//          else {
//
//            // momentum vector contains exp(-i p x)
//            // Divisor 3 for colour index. All three colours on same lattice site get
//            // the same momentum
//            for(int x = 0; x < dim_row; ++x) {
//              mom(x) = momentum[p][x/3];
//            }
//
//            //TODO: check if saving V[t].adjoint in own matrix is faster
//            //TODO: is that efficient?
//            // build vdaggerv = V^dagger exp(-ipx) V 
//            // opposite momentum is just vdaggerv daggered
//            vdaggerv[p][t][dir] = V_t[0].adjoint() * mom.asDiagonal() *  W_t -
//                W_t.adjoint() * mom.asDiagonal() * V_t[0];
//            vdaggerv[nb_mom - p - 1][t][dir] =
//                (-1) *  (vdaggerv[p][t][dir]).adjoint();
//
//          }
//        } // end if displacement

      } // end for momentum
    } // end for displacement

  } // loop over time ends here
          }
        }
      }
    }
  }}// time and momemtum loop ends here

  t2 = clock() - t2;
  std::cout << std::setprecision(1) << "\t\tSUCCESS - " << std::fixed 
    << ((float) t2)/CLOCKS_PER_SEC << " seconds" << std::endl;

}





