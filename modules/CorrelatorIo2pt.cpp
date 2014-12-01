#include "CorrelatorIo2pt.h"
///////////////////////////////////////////////////////////////////////////////
//Helper Routines//////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Calculate momentum squared needed in header
static int square_comp(const std::array<int, 3>& p1,
                       const std::array<int, 3>& p2){
  int square = 0;
  for (size_t i = 0; i < 3; ++i){
    square += p1[i]*p2[i];  
  }
  return square;
}

// Compare two tags of correlation functions
static bool compare_tags(const Tag& tag1, const Tag& tag2){
  bool flag = true;
  if (memcmp(tag1.mom, tag2.mom, sizeof(tag1.mom)) != 0) flag = false;
  if (memcmp(tag1.dis, tag2.dis, sizeof(tag1.dis)) != 0) flag = false;
  if (memcmp(tag1.gam, tag2.gam, sizeof(tag1.gam)) != 0) flag = false;
  return flag;
}

// Set the tag from two operator structures
static void set_tag(Tag& so_si, const pdg& op1, const pdg& op2){
  for (int i = 0; i < 2; ++i){
    so_si.mom[i] = i+1;
    so_si.dis[i] = i+2;
    so_si.gam[i] = i+3;
  }
}
///////////////////////////////////////////////////////////////////////////////
//Write 2pt correlator/////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void write_2pt_lime(const char* filename, const pdg& op_so, const pdg& op_si,
    std::vector<std::complex<double> >& corr){

  // first construct the tag from source and sink operator struct
  Tag so_si;
  set_tag(so_si, op_so, op_si);
  int p_sq = square_comp(op_so.p, op_si.p);
  // setup what is needed for output to lime
  FILE* fp;
  LimeWriter* w;
  // Each message contains to records. Uneven records 
  LimeRecordHeader *id;
  LimeRecordHeader *corr_hd;
  n_uint64_t tag_bytes = sizeof(so_si);
  n_uint64_t data_bytes = corr.size()*2*sizeof(double);
  std::cout << data_bytes << std::endl;
  char headername[100];  
  // Access flags for record and message: MB (Message Begin): 1 if true, ME
  // (Message End): 1 if true
  int MB_flag, ME_flag;
  fp = fopen( filename, "a" );
  w = limeCreateWriter( fp );

  // Record for Tag as struct of 3 2dim integer-arrays
  ME_flag = 0; MB_flag = 1;
  sprintf(headername,"Tag of Correlator with p^2 = %d", p_sq);
  id = limeCreateHeader( MB_flag, ME_flag, headername , tag_bytes );
  limeWriteRecordHeader( id, w );
  limeDestroyHeader( id );
  limeWriteRecordData( &so_si, &tag_bytes, w );

  // Record for correlator belonging to tag
  ME_flag = 1; MB_flag = 0; 
  sprintf(headername,"Correlator with p^2 = %d", p_sq);
  corr_hd = limeCreateHeader( MB_flag, ME_flag, headername, data_bytes );
  limeWriteRecordHeader( corr_hd, w );
  limeDestroyHeader( corr_hd );
  limeWriteRecordData( corr.data(), &data_bytes, w ); 
  limeDestroyWriter( w );
  fclose( fp );

}
///////////////////////////////////////////////////////////////////////////////
//Read 2pt correlator//////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


void read_2pt_lime(const char* filename, const Tag& tag,
    std::vector<std::complex<double> >& corr){

  // search for record in open lime file (compare record)
  Tag read_tag;

  // Begin of Lime Stuff
  FILE* fp;
  LimeReader* r;
  // Each message contains to records. Uneven records 
  n_uint64_t tag_bytes = sizeof(read_tag);
  n_uint64_t data_bytes = corr.size()*2*sizeof(double);
  // Access flags for record and message: -MB (Message Begin): 1 if true 
  //                                      -ME (Message End): 1 if true
  int MB_flag, ME_flag;
  int file_status = 0; 
  int act_status = 0; 
  fp = fopen( filename, "r" );
  r = limeCreateReader( fp );
  while(file_status != LIME_EOF){
    // read tag
    file_status = limeReaderNextRecord(r);
    MB_flag = limeReaderMBFlag(r);
    ME_flag = limeReaderMEFlag(r);
    std::cout << MB_flag << " " << ME_flag << std::endl;
    // get Metainfo of record in read tag
    if(MB_flag == 1 && ME_flag == 0){
      act_status = limeReaderReadData(&read_tag, &tag_bytes, r);
      // compare read_tag with tag
      if(compare_tags(read_tag, tag)){
        // if it is the same write next record to correlation function
        file_status = limeReaderNextRecord(r);
        MB_flag = limeReaderMBFlag(r);
        ME_flag = limeReaderMEFlag(r);
        std::cout << data_bytes << std::endl;
        act_status = limeReaderReadData(corr.data(), &data_bytes, r);
      }
    }
  }
}

