#include "Gaugefield.h"

static GlobalData * const global_data = GlobalData::Instance();

GaugeField::GaugeField(const size_t t0, const size_t tf, const size_t v3,
                      const size_t ndir) : tslices(),
                                           iup(boost::extents[v3][ndir]),
                                           idown(boost::extents[v3][ndir]) {
  tslices.resize(tf-t0+1);
  for(auto& t: tslices) t.resize(boost::extents[v3][ndir]);

}

///////////////////////////////////////////////////////////////////////////////
//Initialize the lookup tables/////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void GaugeField::init(const size_t L1, const size_t L2, const size_t L3){
  int* x0_h = new int[3];
  int* x1_h = new int[3];
  int* x2_h = new int[3];

  int L0_h = L1 * L2;

  for ( int x0 = 0; x0 < L1; ++x0 ) {//loop x0
    x0_h[2] = x0 * L0_h;
    //negative direction (index at lower boundary)
    if ((x0_h[0] = x0 - 1) < 0) x0_h[0] = L0_h * (L1 - 1);
    else x0_h[0] *= L0_h;
    //positive direction (index at upper boundary)
    if ((x0_h[1] = x0 + 1) == L1) x0_h[1] = 0;
    else x0_h[1] *= L0_h;

    for ( int x1 = 0; x1 < L1; ++x1 ) {//loop x1
      x1_h[2] = x1 * L2;
      //neg. dir.
      if ((x1_h[0] = x1 - 1) < 0) x1_h[0] = L2 * (L1 - 1);
      else x1_h[0] *= L2;
      //pos. dir.
      if ((x1_h[1] = x1 + 1) == L1) x1_h[1] = 0;
      else x1_h[1] *= L2;

      for ( int x2 = 0; x2 < L2; ++x2 ) {//loop x2
        x2_h[2] = x2;
        //neg. dir.
        if ((x2_h[0] = x2 - 1) < 0) x2_h[0] = L2 -1;
        //pos. dir.
        if ((x2_h[1] = x2 +1) == L2) x2_h[1] = 0;
        //overall volume index
        int i = x0_h[2] + x1_h[2] + x2_h[2];
        //std::cout << x0 << " " << x1 << " " << x2 << " " << i << std::endl;
        //upwards
        iup[i][0] = x0_h[1] + x1_h[2] + x2_h[2];
        iup[i][1] = x0_h[2] + x1_h[1] + x2_h[2];
        iup[i][2] = x0_h[2] + x1_h[2] + x2_h[1];
        //downwards
        idown[i][0] = x0_h[0] + x1_h[2] + x2_h[2];
        idown[i][1] = x0_h[2] + x1_h[0] + x2_h[2];
        idown[i][2] = x0_h[2] + x1_h[2] + x2_h[0];
      }//end loop x2
    }//end loop x1
  }//end loop x0
  delete x0_h;
  delete x1_h;
  delete x2_h;
 
}

//get index in positive direction
int GaugeField::get_up(const int pos, const int dir){
  return iup[pos][dir];
}
//get index in negative direction
int GaugeField::get_dn(const int pos, const int dir){
  return idown[pos][dir];
}
///////////////////////////////////////////////////////////////////////////////
//Transform one timeslice from lime array to array_3cd_d2_eigen////////////////
///////////////////////////////////////////////////////////////////////////////
//mapping from gauge config to Eigen 3x3 complex matrix arrays
void GaugeField::map_timeslice_to_eigen(const size_t t, const double* timeslice) {
  int L1 = global_data->get_Lx();
  int L2 = global_data->get_Ly();
  int L3 = global_data->get_Lz();
  
  const int V3 = L1 * L2 * L3;
  const int V_TS = global_data->get_V_TS();
  
  //Number of directions
  const int NDIR = 4;
  //Number of colors
  const int NCOL = 3;

  //read in elements
  int el_input = 0;
  for (int z = 0; z < L3; ++z) {//spatial loops
    for (int y = 0; y < L2; ++y) {
      for (int x = 0; x < L1; ++x) {
        for (int mu = 1; mu < 4; ++mu) {//direction loop
          std::complex< double > array[9];
          for (int a = 0; a < 3; ++a) {//colour loops
            for (int b = 0; b < 3; ++b) {
              //timeslice index of real part
              int ind_r = z*V_TS/L3+y*V_TS/(L3*L2)+x*V_TS/(V3)+
                mu*V_TS/(V3*NDIR)+a*V_TS/(V3*NDIR*NCOL)
                +b*V_TS/(V3*NDIR*NCOL*NCOL)+0;
              //timeslice index of imaginary part
              int ind_i = z*V_TS/L3+y*V_TS/(L3*L2)+x*V_TS/(V3)+
                mu*V_TS/(V3*NDIR)+a*V_TS/(V3*NDIR*NCOL)
                +b*V_TS/(V3*NDIR*NCOL*NCOL)+1;
		
              std::complex<double> pair(timeslice[ind_r], timeslice[ind_i]);
              //array to be mapped to Eigen Array
              array[3*b+a] = pair;
              ++el_input;
            }
          }
          Eigen::Map<Eigen::Matrix3cd> dummy(array);
          //spatial index
          int ind = z*L2*L1+y*L1+x;
           tslices.at(t)[ind][mu-1] = dummy;
        }
      }
    }
  }
  std::cout << el_input << " doubles read in from ildg timeslice " << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
//Useful helpers for smearing and gaugetrafos//////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
std::array< int, 2 > get_dirs(int mu) {
  std::array<int,2> dirs;  
  if (mu == 0) {
    dirs.at(0) = 1;
    dirs.at(1) = 2;
  }
  else if (mu == 1) {
    dirs.at(0) = 0;
    dirs.at(1) = 2; 
  }
  else if (mu == 2) {
    dirs.at(0) = 0;
    dirs.at(1) = 1; 
  }
  
  return(dirs);
}

//storage map for dec_timeslice
static int decor (int dir, int smear_plane) {
  int ret;
  if (dir == 0) {
    if (smear_plane == 1) ret = 0;
    else ret = 3; 
  }
  else if (dir == 1) {
    if (smear_plane == 0) ret = 1;
    else ret = 4;
  }
  else {
    if (smear_plane == 0) ret = 2;
    else ret = 5;
  }
  return ret;
}

static Eigen::Matrix3cd proj_to_su3_imp(Eigen::Matrix3cd& in){
  //avoid possible aliasing issues:
  Eigen::Matrix3cd out = ( (in.adjoint() * in).sqrt() ).inverse();
  in = in*out;
  std::complex<double> det = 1./pow(in.determinant(),1./3.);
  in *= det;
  return in;
  
}
///////////////////////////////////////////////////////////////////////////////
//Smearing methods/////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//Stout-Smearing
void GaugeField::smearing_stout(const size_t t, const double rho, const size_t iter) {

  //parameter passing still to be improved
  const int LX = global_data -> get_Lx();
  const int LY = global_data -> get_Ly();
  const int LZ = global_data -> get_Lz();
  const int V3 = LX * LY * LZ;

  std::complex<double> im_half(0,0.5);
  array_3cd_d2_eigen eigen_timeslice_ts(boost::extents[V3][3]);
 // Eigen::Matrix3cd **eigen_timeslice_ts = new Eigen::Matrix3cd *[V3]; 
 // for ( auto i = 0; i < V3; ++i ) {
 //   eigen_timeslice_ts[i] = new Eigen::Matrix3cd[3];
 // }
  for (int j = 0; j < iter; ++j) {
    for (int i = 0; i < V3; ++i) {
      for (int dir = 0; dir < 3; ++dir) {
        Eigen::Matrix3cd staple = Eigen::Matrix3cd::Zero(); //Holding all smearing matrices for one link
        //filling each element of smearer using 3 links
        //For each link calculate staples summing them up
        for (int not_dir = 0; not_dir < 3; ++not_dir) {
          if (dir != not_dir) {
            int mu = iup[i][not_dir];
            int nu = iup[mu][dir];
            int eta = idown[nu][not_dir];
            //Staples in positive direction
            staple += (tslices.at(t))[i][not_dir]*
                      ( (tslices.at(t))[mu][dir]*
                      ((tslices.at(t))[eta][not_dir].adjoint()));

            mu = idown[i][not_dir];
            nu = iup[mu][dir];
            //Staples in negative direction
            staple += ( (tslices.at(t))[mu][not_dir].adjoint() )*
                      ( (tslices.at(t))[mu][dir] * (tslices.at(t))[nu][dir] );
          }
        }
        Eigen::Matrix3cd omega = ( staple * ( rho/4. ) ) *
                                 ( (tslices.at(t))[i][dir].adjoint() );
        Eigen::Matrix3cd q = ( ( omega.adjoint() - omega ) * 0.5 ) -
                             ( ( ( ( omega.adjoint() - omega ).trace() ) *
                              Eigen::Matrix3cd::Identity() ) * (1./6.) );
        eigen_timeslice_ts[i][dir] = ( ( q*(-1) ).exp() ) *
                                     (tslices.at(t))[i][dir];
      }
    }
    for ( auto i = 0; i < V3; ++i ) {
      for ( auto mu = 0; mu < 3; ++mu) {
        tslices.at(t)[i][mu] = eigen_timeslice_ts[i][mu];
      }
    }
  }
}

//APE-Smearing
void GaugeField::smearing_ape(const size_t t, const double alpha_1, const size_t iter){

  //parameter passing still to be improved
  const int LX = global_data -> get_Lx();
  const int LY = global_data -> get_Ly();
  const int LZ = global_data -> get_Lz();
  const int V3 = LX * LY * LZ;
  //Eigen::Matrix3cd **eigen_timeslice_ts = new Eigen::Matrix3cd *[V3]; 
  //for ( auto i = 0; i < V3; ++i ) {
  //  eigen_timeslice_ts[i] = new Eigen::Matrix3cd[3];
  //}
  //temporal timeslice from decorated links
  array_3cd_d2_eigen eigen_timeslice_ts(boost::extents[V3][3]);  
  for (int j = 0; j < iter; ++j) {
    for (int i = 0; i < V3; ++i) {
      for (int dir = 0; dir < 3; ++dir) {
        Eigen::Matrix3cd staple = Eigen::Matrix3cd::Zero(); //Holding all smearing matrices for one link
        //filling each element of smearer using 3 links
        //Position indices mu nu eta
        for (int not_dir = 0; not_dir < 3; ++not_dir) {
          if (dir != not_dir) {
            int mu = iup[i][not_dir];
            int nu = iup[mu][dir];
            int eta = idown[nu][not_dir];
            //Staples in positive direction
            staple += (tslices.at(t))[i][not_dir]*
              ( (tslices.at(t))[mu][dir] * ((tslices.at(t))[eta][not_dir].adjoint()) );

            mu = idown[i][not_dir];
            nu = iup[mu][dir];
            //Staples in negative direction
            staple += (tslices.at(t)[mu][not_dir].adjoint()) *
              ((tslices.at(t))[mu][dir] * (tslices.at(t))[nu][dir]);
          }
        }
        eigen_timeslice_ts[i][dir] = (tslices.at(t)[i][dir] * (1.-alpha_1)) + (staple * alpha_1/4.);
      }
    }
    for ( auto i = 0; i < V3; ++i ) {
      for ( auto mu = 0; mu < 3; ++mu) {
        tslices.at(t)[i][mu] = proj_to_su3_imp(eigen_timeslice_ts[i][mu]);
      }
    }
  }
}

//HYP-Smearing

void GaugeField::smearing_hyp( const size_t t, const double alpha_1, const double alpha_2,
                               const size_t iter) {
  //parameter passing still to be improved
  const int LX = global_data -> get_Lx();
  const int LY = global_data -> get_Ly();
  const int LZ = global_data -> get_Lz();
  const int V3 = LX * LY * LZ;
  //temporal timeslice twice the size for decorated links
  array_3cd_d2_eigen dec_timeslice(boost::extents[V3][6]);

  //temporal timeslice from decorated links
  array_3cd_d2_eigen eigen_timeslice_ts(boost::extents[V3][3]);  
  //temporal integers holding directions for decorated smearing
  int mu, nu, eta;
  for (int run = 0; run < iter; ++run) {
    //calculate inner staple from original timeslice, store in dec_timeslice each link can get smeared in two planes
    for (int vol = 0; vol < V3; ++vol) {
      for (int dir = 0; dir < 3; ++dir) {
        //inner staple
        Eigen::Matrix3cd inner_staple = Eigen::Matrix3cd::Zero();
        std::array< Eigen::Matrix3cd, 2 > tmp_staples;
        std::array<int, 2> perpendics = get_dirs( dir );
        for (auto it_perp_dir = perpendics.begin(); it_perp_dir != perpendics.end(); ++it_perp_dir ) {
          int perp_dir = *it_perp_dir;
          //up-type smearing_indices
          mu = iup[vol][perp_dir];
          nu = iup[mu][dir];
          eta = iup[vol][dir];
          //product of up matrices
          inner_staple = tslices.at(t)[vol][perp_dir] * ( tslices.at(t)[mu][dir] *
                        ( tslices.at(t)[eta][perp_dir].adjoint() ) ); 
          //down-type smearing indices
          mu = idown[vol][perp_dir];
          nu = iup[mu][dir];
          //eta is same endpoint no adjoint necessary here
          //product of down matrices
          inner_staple += ( tslices.at(t)[mu][perp_dir].adjoint() ) *
                          ( tslices.at(t)[mu][dir] * tslices.at(t)[nu][perp_dir] );
          //Careful placement of decorated links in dec_timeslices:
          //dir=0 has placement in dec_dir = 0 (smeared in 1 plane)
          //                       dec_dir = 3 (smeared in 2 plane)
          //dir=1 has placement in dec_dir = 1 (smeared in 0 plane)
          //                       dec_dir = 4 (smeared in 2 plane)
          //dir=2 has placement in dec_dir = 2 (smeared in 0 plane)
          //                       dec_dir = 5 (smeared in 1 plane)
          Eigen::Matrix3cd stac = ( tslices.at(t)[vol][dir] *
                              (1-alpha_2) ) +  ( inner_staple * alpha_2/2.);  
          int n_el = it_perp_dir - perpendics.begin();
          tmp_staples.at(n_el) = proj_to_su3_imp(stac);
        }

        //staple link in direction dir in non participating and negative directions
        dec_timeslice[vol][dir] = tmp_staples.at(0);
        dec_timeslice[vol][dir+3] = tmp_staples.at(1);
      }
    }
    //calculate outer staple from dec_timeslice as modified ape-smearing

      for (int i = 0; i < V3; ++i) {
        for (int dir = 0; dir < 3; ++dir) {
          Eigen::Matrix3cd outer_staple = Eigen::Matrix3cd::Zero(); //Holding all smearing matrices for one link

          //filling each element of smearer using 3 links
          //debugging link
          Eigen::Matrix3cd outer_staple_test;
          for (int not_dir = 0; not_dir < 3; ++not_dir) {
            if (dir != not_dir) {

              //calculate plane in which was smeared
              int plane = ( (dir+1) ^ (not_dir+1) ) - 1;
              mu = iup[i][not_dir];
              nu = iup[mu][dir];
              eta = idown[nu][not_dir];

              //Staples in positive direction
              //replace directions by appropriate decor 
              int a,b;
              a = decor(not_dir,plane);
              b = decor(dir,plane);
              outer_staple += dec_timeslice[i][a]*
                (dec_timeslice[mu][b]*(dec_timeslice[eta][a].adjoint()));
              mu = idown[i][not_dir];
              nu = iup[mu][dir];
            
              //Staples in negative direction
              outer_staple += (dec_timeslice[mu][a].adjoint())*
                (dec_timeslice[mu][b]*dec_timeslice[nu][a]); //structure has to be a_dag, b, a
            }
          }
          eigen_timeslice_ts[i][dir] = (tslices.at(t)[i][dir] * (1.-alpha_1)) +
                                       (outer_staple * alpha_1/4.);
        }
      }
      for ( int i = 0; i < V3; ++i ) {
        for ( int mu = 0; mu < 3; ++mu) {
          tslices.at(t)[i][mu] = proj_to_su3_imp(eigen_timeslice_ts[i][mu]);
        }
      }
  }
  //clean up
}

///////////////////////////////////////////////////////////////////////////////
///Displacement routines, returning one Eigenvector////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//Derivative, toogle symmetrization via sym
Eigen::MatrixXcd GaugeField::disp(const Eigen::MatrixXcd& v,
                                     const size_t t, const size_t dir, bool sym ) {
  //parameter passing still to be improved
  const int LX = global_data -> get_Lx();
  const int LY = global_data -> get_Ly();
  const int LZ = global_data -> get_Lz();
  const int V3 = LX * LY * LZ;

  //Information on Matrix size
  const int dim_row = V3*3;
  const int dim_col = v.cols();
  //Loop over all eigenvectors in 
//    std::cout << "eigenvector: " << i << std::endl;
    //storing eigenvector
    Eigen::VectorXcd in(dim_row);
    Eigen::MatrixXcd out(dim_row, dim_col);
  for(int ev=0; ev < dim_col; ++ev){ 
    in = v.col(ev);

    //Displace eigenvector
    for (int spatial_ind = 0; spatial_ind < V3; ++spatial_ind) {
      //std::cout << "x: " << spatial_ind << std::endl;
      Eigen::Vector3cd tmp;
      Eigen::Vector3cd quark_up;
      Eigen::Vector3cd quark_down;

      //determine needed indices from lookup tables;
      int up_ind = iup[spatial_ind][dir];
      int down_ind = idown[spatial_ind][dir];

      quark_up = in.segment(3*up_ind,3);
      quark_down = in.segment(3*down_ind,3);
      if(sym) {
        tmp = 0.5*( ( (tslices.at(t))[spatial_ind][dir] * quark_up) - 
            ( ( (tslices.at(t))[down_ind][dir].adjoint() ) * quark_down) ); 
      }
      else { 
        Eigen::Vector3cd quark_point = in.segment(3*spatial_ind,3);
        tmp =  ( (tslices.at(t))[spatial_ind][dir] * quark_up) - quark_point ;
      }
      (out.col(ev)).segment(3*spatial_ind,3) = tmp;
    }//end spatial loop
  }//end eigenvector loop
  return out;
}

///////////////////////////////////////////////////////////////////////////////
///Gaugefield transformations//////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//Gauge-transform config and store transformed in gauge_config
void GaugeField::trafo_ts( Eigen::Matrix3cd* Omega) {

  int V3 = pars -> get_int("V3");
  /*
  //SU(3)-gauge field of same size as config
  Eigen::Matrix3cd* Omega = new Eigen::Matrix3cd [V3];
  build_gauge_array(Omega);
  write_gauge_matrices("ts_trafo_log.bin", Omega);
  */
  for ( int i = 0; i < V3; ++i ) {
    for ( int mu = 0; mu < 3; ++mu ) {
      int j = lookup -> get_up(i,mu);
      this -> set_gauge(i,mu, (Omega[i].adjoint())* ( (this -> get_gauge(i,mu) ) * Omega[j]) );
    //std::cout << "gauge config: \n" << gauge_config[i][mu] << "\n\n";
     // std::cout << "config: \n" << config[i][mu] << "\n\n\n";
    }
  }
//  delete[] Omega;
  /*for ( auto i = 0; i < V3; ++i ) {
    for ( auto mu = 0; mu < 3; ++mu ) {
      config[i][mu] = gauge_config[i][mu];
    }
  }*/

}
///////////////////////////////////////////////////////////////////////////////
///Data IO from and to files///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//Read in gauge field to vector of timeslices
void GaugeField::read_gauge_field(const size_t slice_i, const size_t slice_f){
  char filename[200];
  sprintf(filename, "/hiskp2/gauges/test4x4x4x4/conf.1000");
  const int vol4 = global_data -> get_V_for_lime();
  const int V_TS = global_data -> get_V_TS();
  double* configuration = new double[vol4];
  read_lime_gauge_field_doubleprec_timeslices(configuration, filename, slice_i, slice_f);
  for (auto t = slice_i; t <= slice_f; ++t) {
    double* timeslice = configuration + V_TS*t;
    map_timeslice_to_eigen(t, timeslice);
  }   
  delete[] configuration;
}

//Read in the gaugefield, deprecated tmLqcd routine
void GaugeField::read_lime_gauge_field_doubleprec_timeslices(double* gaugefield,
                                                             const char* filename,
                                                             const size_t slice_i,
                                                             const size_t slice_f) {

  try{
    const int Lt = global_data -> get_Lt();
    const int Lx = global_data -> get_Lx();
    const int Ly = global_data -> get_Ly();
    const int Lz = global_data -> get_Lz();
    const int verbose = global_data -> get_verbose();

    //const int slice_i = 0;
    //const int slice_f = Lt+1;

    FILE * ifs;
    int t, x, y, z, status;
    n_uint64_t bytes;
    char * header_type;
    LimeReader * limereader;
    double tmp[72], tmp2[72];
    int words_bigendian;

    if(verbose){
      printf("reading gauge fields from files:\n");
    }
    else{
      printf("\treading gauge fields:");
    }

    words_bigendian = big_endian();
    ifs = fopen(filename, "r");
    if(ifs == (FILE *)NULL) {
      fprintf(stderr, "Could not open file %s\n Aborting...\n", filename);
      exit(500);
    }
    limereader = limeCreateReader( ifs );
    if( limereader == (LimeReader *)NULL ) {
      fprintf(stderr, "Unable to open LimeReader\n");
      exit(500);
    }
    while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
      if(status != LIME_SUCCESS ) {
        fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n",
            status);
        status = LIME_EOF;
        break;
      }
      header_type = limeReaderType(limereader);
      if(strcmp("ildg-binary-data",header_type) == 0) break;
    }
    if(status == LIME_EOF) {
      fprintf(stderr, "no ildg-binary-data record found in file %s\n",filename);
      limeDestroyReader(limereader);
      fclose(ifs);
      exit(-2);
    }
    bytes = limeReaderBytes(limereader);
    if(bytes != (n_uint64_t)Lx*Ly*Lz*Lt*72*(n_uint64_t)sizeof(double)) {
      if(bytes != (n_uint64_t)Lx*Ly*Lz*Lt*72*(n_uint64_t)sizeof(double)/2) {
        fprintf(stderr, "Probably wrong lattice size or precision (bytes=%lu) in file %s expected %lu\n",
            (n_uint64_t)bytes, filename,
            (n_uint64_t)Lx*Ly*Lz*Lt*72*(n_uint64_t)sizeof(double));
        fprintf(stderr, "Aborting...!\n");
        fflush( stdout );
        exit(501);
      }
      else {
        fclose(ifs);
        fprintf(stderr, "single precision read!\n");

        //fprintf(stderr, "Not implemented!\n");
        exit(EXIT_FAILURE);
        //read_lime_gauge_field_singleprec(gaugefield, filename, Lt, Lx, Ly, Lz);
        return;
      }
    }

    bytes = (n_uint64_t)72*sizeof(double);

    for(t = 0; t < Lt; t++) {
      for(z = 0; z < Lz; z++) {
        for(y = 0; y < Ly; y++) {
          for(x = 0; x < Lx; x++) {

            // check for endianess and reading in data
            // the pointer limereader is internally increased by bytes
            // in the limeReaderReadData function
            if(!words_bigendian) {
              status = limeReaderReadData(tmp, &bytes, limereader);
              byte_swap_assign(tmp2, tmp, 72);
            }
            else
              status = limeReaderReadData(tmp2, &bytes, limereader);
            // check if reading was successfull
            if(status < 0 && status != LIME_EOR) {
              fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n",
                  status, filename);
              exit(500);
            }

            // we just want to read in data in the specific range of timeslices
            // must be here because the limereader pointer must be increased 
            // correctly
            // could be done with much more performance but might be tricky to 
            // do correctly
            if(t<slice_i || t>slice_f)
              continue;

            // copy of link variables from tmp2 into config
            // ILDG has mu-order: x,y,z,t so it is changed here to: t,x,y,z !
            const size_t p = (size_t) ( ((t-slice_i)*Lx*Lz*Lz + 
                  x*Ly*Lz + y*Lz + z) * 72); // position in config
            size_t k = 0;
            for(size_t mu = 1; mu <= 4; mu++) { // mu=4 is for the shift of U_t
              size_t index;
              if (mu != 4)
                index = p + mu*18; // for U_x, U_y and U_z
              else
                index = p; // U_t is copied into the beginning of
              // the (config+p) array

              for(size_t i = 0; i < 3; i++) {
                for(size_t j = 0; j < 3; j++) {
                  gaugefield[index+6*i+2*j] = tmp2[2*k];
                  gaugefield[index+6*i+2*j+1] = tmp2[2*k+1];
                  k++;
                }
              }

            } // loop over mu ends here

          } // loop over position space ends here
        }
      }
    }
    if(status < 0 && status != LIME_EOR) {
      fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n",
          status, filename);
      exit(500);
    }
    limeDestroyReader(limereader);
    fclose(ifs);
    return;

  }
  catch(std::exception& e){
    std::cout << e.what() << "in: ReadWrite::read_lime_gauge_field_doubleprec_timeslices\n";
    exit(0);
  }

}







