#include "Gaugefield.h"

static GlobalData * const global_data = GlobalData::Instance();

GaugeField::GaugeField(const size_t t0, const size_t tf, const size_t v3, const size_t ndir) : tslices() {
  tslices.resize(tf-t0+1);
  for(auto& t: tslices) t.resize(boost::extents[v3][ndir]);
}


void GaugeField::read_lime_gauge_field_doubleprec_timeslices(double* gaugefield, const char* filename, const size_t slice_i, const size_t slice_f) {

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









