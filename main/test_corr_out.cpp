/***************************************************************************/
//
// Short demonstration of gaugefield class
//
// Author: Christopher Helmes
//
/***************************************************************************/
#include <array>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <typeinfo>
#include <vector>

//#include "CorrelatorIo.h"
#include "BasicOperator.h"
#include "GlobalData.h"
#include "lime.h"
#include "typedefs.h"

// Function to set the Tag on correlator

 struct OpInfo{
    size_t id;
    std::array<int, 3> p;
    std::array<int, 3> dis_id;
    std::array<int, 4> gamma_id;

  };
  // This struct holds metainfo for lime files index 0: src, index 1: sink
 struct Tag{
    int mom[2];
    int dis[2];
    int gam[2];
  };
void set_tag(Tag& so_si);

void build_tag_name(char* name, Tag& tag);

int main(int ac, char* av[]) {

// dummy data
  // initialize one OpInfo struct
  OpInfo op_info;
  // one three momentum
  std::array<int,3> mom {{0,1,2}};
  // one displacement vector
  std::array<int,3> der {{3,4,5}};
  // gamma configuration 
  std::array<int,4> dirac {{5,4,4,4}};
  // initialize operator info
  op_info.p = mom;
  op_info.dis_id = der;
  op_info.gamma_id = dirac;
  // print everything as a check
  std::cout << "p:"; 
  for (auto mom_comp : op_info.p) 
    std::cout << " " << mom_comp;
    std::cout << std::endl;
  std::cout << "dis:"; 
  for (auto der : op_info.dis_id) 
    std::cout << " " << der;
    std::cout << std::endl;
  std::cout << "gamma:"; 
  for (auto dir : op_info.gamma_id) 
    std::cout << " " << dir;
    std::cout << std::endl;
  //generate random numbers for test
  float re[96];
  float im[96];
  rlxs_init(2,1337);
  ranlxs(re, 96);
  ranlxs(im, 96);

  const int T = 96;
  std::array<std::complex<double>,T > dummy_corr;
  for (int t = 0; t < T; ++t){
    dummy_corr.at(t) = std::complex<double>(re[t], im[t]);
  }


  // initialize MetaInfo
  Tag so_si;
  // fill MetaInfo
  set_tag(so_si);
  // Calculate momentum squared needed in header
  int p_sq = 0;
  for (auto p_comp : op_info.p) p_sq += p_comp*p_comp;  
  // Begin of Lime Stuff
  FILE *fp;
  LimeWriter *w;
  // Each message contains to records. Uneven records 
  LimeRecordHeader *id;
  LimeRecordHeader *corr;
  n_uint64_t bytes = T*2*sizeof(double);
  n_uint64_t tag_bytes = sizeof(so_si);
  // Access flags for record and message: MB (Message Begin): 1 if true, ME
  // (Message End): 1 if true
  int MB_flag, ME_flag;
  fp = fopen( "test_file", "w" );
  w = limeCreateWriter( fp );
  for(int msg = 0; msg < 7; ++msg){
  // write metainfo as struct of 3 2dim integer-arrays
  ME_flag = 0; MB_flag = 1;
  char headername[100];  
  sprintf(headername,"Tag of corr with p^2 = %d", p_sq);
  id = limeCreateHeader( MB_flag, ME_flag, headername , tag_bytes );
  limeWriteRecordHeader( id, w );
  limeDestroyHeader( id );
  limeWriteRecordData( &so_si, &tag_bytes, w );

  // Write the correlator to that metainfo
  ME_flag = 1; MB_flag = 0; 
  sprintf(headername,"Correlator with p^2 = %d", p_sq);
  corr = limeCreateHeader( MB_flag, ME_flag, headername, bytes );
  limeWriteRecordHeader( corr, w );
  limeDestroyHeader( corr );
  limeWriteRecordData( &dummy_corr, &bytes, w ); 
  }
  limeDestroyWriter( w );
  fclose( fp );

  // End of Lime Stuff
  return 0;

}

void set_tag(Tag& so_si){
  for (int i = 0; i < 2; ++i){
    so_si.mom[i] = i+1;
    so_si.dis[i] = i+2;
    so_si.gam[i] = i+3;
  }
}

