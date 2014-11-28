#include "CorrelatorIo2pt.h"
///////////////////////////////////////////////////////////////////////////////
//Write 2pt correlator/////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void write_2pt_lime(const char* filename, const pdg& op_so, const pdg& op_si,
                    const std::array<std::complex<double> >& corr){

  // first construct the tag from source and sink operator struct
  Tag so_si;
  build_tag(so_si, op_so, op_si);
  // setup what is needed for output to lime
  FILE* fp;
  LimeWriter* w;
  // Each message contains to records. Uneven records 
  LimeRecordHeader *id;
  LimeRecordHeader *corr;
  n_uint64_t tag_bytes = sizeof(so_si);
  n_uint64_t data_bytes = corr->size()*2*sizeof(double);
  char headername[100];  
  // Access flags for record and message: MB (Message Begin): 1 if true, ME
  // (Message End): 1 if true
  int MB_flag, ME_flag;
  fp = fopen( filename, "w" );
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
  corr = limeCreateHeader( MB_flag, ME_flag, headername, bytes );
  limeWriteRecordHeader( corr, w );
  limeDestroyHeader( corr );
  limeWriteRecordData( &dummy_corr, &bytes, w ); 
  limeDestroyWriter( w );
  fclose( fp );

}
///////////////////////////////////////////////////////////////////////////////
//Read 2pt correlator//////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


void read_2pt_lime(const char* filename, const Tag& tag,
                   std::array<std::complex<double> >& corr){

  // search for record in open lime file (compare record)
  Tag read_tag;

  // Begin of Lime Stuff
  FILE* fp;
  LimeReader* r;
  // Each message contains to records. Uneven records 
  LimeRecordHeader *id;
  LimeRecordHeader *corr;
  n_uint64_t tag_bytes = sizeof(read_tag);
  n_uint64_t data_bytes = corr->size()*2*sizeof(double);
  // Access flags for record and message: -MB (Message Begin): 1 if true 
  //                                      -ME (Message End): 1 if true
  int MB_flag, ME_flag;
  int file_status = 0; 
  int act_status = 0; 
  fp = fopen( "test_file", "r" );
  r = limeCreateReader( fp );
  while(file_status != LIME_EOF){
    // read tag
    file_status = limeReaderNextRecord(r);
    MB_flag = limeReaderMBFlag(r);
    ME_flag = limeReaderMEFlag(r);
    // get Metainfo of record in read tag
    if(MB_flag == 1 && MB_flag == 0)
      act_status = limeReaderReadData(&read_tag, &tag_bytes, r);
    // compare read_tag with tag
    if(read_tag == tag){
      // if it is the same write next record to correlation function
      file_status = limeReaderNextRecord(r);
      act_status = limeReaderReadData(corr, &_bytes, r);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//Helper Routines//////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void set_tag(Tag& so_si, const pdg& op1, const pdg& op2){

  for (int i = 0; i < 2; ++i){
    so_si.mom[i] = i+1;
    so_si.dis[i] = i+2;
    so_si.gam[i] = i+3;
  }
}
