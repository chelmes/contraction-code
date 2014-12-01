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

  // This struct holds metainfo for lime files index 0: src, index 1: sink
 struct Tag{
    int mom[2];
    int dis[2];
    int gam[2];
  };
void set_tag(Tag& so_si);
void build_tag_name(char* name, Tag& tag);

int main(int ac, char* av[]) {


// goal: get one message of lime file to array dummy corr and according tag to
// tag structure
  // Correlator structure
  const int T = 96;
  std::array<std::complex<double>,T > dummy_corr;
  // initialize MetaInfo
  Tag so_si;
  
  // Begin of Lime Stuff
  FILE* fp;
  LimeReader* r;
  // Each message contains to records. Uneven records 
  LimeRecordHeader *id;
  LimeRecordHeader *corr;
  n_uint64_t bytes = T*2*sizeof(double);
  n_uint64_t tag_bytes = sizeof(so_si);
  // Access flags for record and message: MB (Message Begin): 1 if true, ME
  // (Message End): 1 if true
  int MB_flag, ME_flag;
  int status = 0; 
  fp = fopen( "test_file", "r" );
  r = limeCreateReader( fp );
  // read tag
  status = limeReaderNextRecord(r);
  MB_flag = limeReaderMBFlag(r);
  ME_flag = limeReaderMEFlag(r);
  if(MB_flag == 1 && MB_flag == 0)
    status = limeReaderReadData(&so_si, &tag_bytes, r);
  else if(MB_flag == 0 && MB_flag == 1) 
    status = limeReaderReadData(&dummy_corr, &bytes, r);
  std::cout << limeReaderMBFlag(r) << " " << limeReaderMEFlag(r) << std::endl; 
  std::cout << so_si.mom[0] << " " << so_si.mom[1] << std::endl;
  status = limeReaderNextRecord(r);
  for (auto& el : dummy_corr) std::cout << el << std::endl;

//  // Write the correlator to that metainfo
//  ME_flag = 1; MB_flag = 0; 
//  sprintf(headername,"Correlator with p^2 = %d", p_sq);
//  corr = limeCreateHeader( MB_flag, ME_flag, headername, bytes );
//  limeWriteRecordHeader( corr, w );
//  limeDestroyHeader( corr );
//  limeWriteRecordData( &dummy_corr, &bytes, w ); 
  limeDestroyReader( r );
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

