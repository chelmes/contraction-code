#include "IoHelpers.h"

///////////////////////////////////////////////////////////////////////////////
// Helper Routines ////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// swap endianess of global data


//Tag handling/////////////////////////////////////////////////////////////////

// Calculate p^2
 int square_comp(const std::array<int, 3>& p1,
                       const std::array<int, 3>& p2){
  int square = 0;
  for (size_t i = 0; i < 3; ++i){
    square += p1[i]*p2[i];  
  }
  return square;
}

// Compare two tags of correlation functions
 bool compare_tags(const Tag& tag1, const Tag& tag2){
  bool flag = true;
  if (memcmp(tag1.mom, tag2.mom, sizeof(tag1.mom)) != 0) flag = false;
  if (memcmp(tag1.dis, tag2.dis, sizeof(tag1.dis)) != 0) flag = false;
  if (memcmp(tag1.gam, tag2.gam, sizeof(tag1.gam)) != 0) flag = false;
  return flag;
}

// Set the tag from two operator structures
 void set_tag(Tag& so_si, const pdg& op1, const pdg& op2){
  for (int i = 0; i < 2; ++i){
    so_si.mom[i] = i+1;
    so_si.dis[i] = i+2;
    so_si.gam[i] = i+3;
  }
}

// File handling ///////////////////////////////////////////////////////////////

// write global checksum, and runinfo in new file
void write_1st_msg(const char* filename, GlobalDat& dat,
                          size_t chksum){
  // open file
  FILE* fp;
  fp = fopen(filename, "w");
  LimeWriter* w;
  w = limeCreateWriter(fp);
  // first record is global checksum
  int MB_flag = 1;
  int ME_flag = 0;
  LimeRecordHeader* chk;
  n_uint64_t chk_bytes = sizeof(chksum);
  chk = limeCreateHeader(MB_flag, ME_flag, "Global Checksum", chk_bytes);
  limeWriteRecordHeader(chk, w);
  limeDestroyHeader(chk);
  limeWriteRecordData(&chksum, &chk_bytes, w);
  
  // second record is data on calculation
  MB_flag = 0; ME_flag = 1;
  
  LimeRecordHeader* run_id;
  n_uint64_t run_id_bytes = sizeof(dat);
  run_id = limeCreateHeader(MB_flag, ME_flag, "Runinfo", run_id_bytes);
  limeWriteRecordHeader(run_id,w);
  limeDestroyHeader(run_id);
  limeWriteRecordData(&dat, &run_id_bytes, w);
  limeDestroyWriter(w);
  fclose(fp);
} 

//append correlation functions to file
void append_msgs(const char* filename, std::vector<vec>& corr, std::vector<Tag>& tags,
              LimeWriter* w, bool be){
  // Each message contains three records:
  // 1st: Checksum for the Correlator
  // 2nd: Tag for the Correlator
  // 3rd: Correlationfunction itself
  LimeRecordHeader* corr_chk;
  LimeRecordHeader* id;
  LimeRecordHeader* corr_hd;

  boost::uint64_t corr_chksum;
  n_uint64_t tag_bytes = sizeof(tags[0]);
  n_uint64_t data_bytes = (corr[0]).size()*2*sizeof(double);
  char headername[100];  
  // Access flags for record and message: MB (Message Begin): 1 if true, ME
  // (Message End): 1 if true
  int MB_flag, ME_flag;
  for(size_t el = 0; el < corr.size(); ++el){

    // 1st record for checksum
    corr_chksum = checksum <vec> (corr[el], corr[el].size());
    if(be) corr_chksum = swap_endian<boost::uint64_t>(corr_chksum);
    n_uint64_t chk_bytes = sizeof(corr_chksum);
    ME_flag = 0; MB_flag = 1;
    corr_chk = limeCreateHeader(MB_flag, ME_flag,"Correlator checksum", chk_bytes);
    limeWriteRecordHeader( corr_chk, w );
    limeDestroyHeader( corr_chk );
    limeWriteRecordData( &corr_chksum, &chk_bytes, w );

    // 2nd record for Tag as struct of 3 2dim integer-arrays
    ME_flag = 0; MB_flag = 0;
    sprintf(headername,"Tag of Correlator with p^2 = %zd", el);
    id = limeCreateHeader( MB_flag, ME_flag, headername , tag_bytes );
    limeWriteRecordHeader( id, w );
    limeDestroyHeader( id );
    limeWriteRecordData( &tags[el], &tag_bytes, w );

    // 3rd record for correlator belonging to tag
    ME_flag = 1; MB_flag = 0; 
    sprintf(headername,"Correlator with p^2 = %zd", el);
    corr_hd = limeCreateHeader( MB_flag, ME_flag, headername, data_bytes );
    limeWriteRecordHeader( corr_hd, w );
    limeDestroyHeader( corr_hd );
    limeWriteRecordData( corr[el].data(), &data_bytes, w ); 
  }
}

// Check checksums
void file_check(const size_t glob_check,
                const std::vector<boost::uint64_t>& checksums,
                const std::vector<cmplx>& correlators){
  size_t tmp = 0;
  size_t bytecount = correlators.size()*sizeof(double);
  
  tmp = checksum<std::vector<cmplx> > (correlators, bytecount);
  if (big_endian()) swap_endian <size_t> (tmp);
  std::cout << bytecount << std::endl;
  if (tmp == glob_check) {
    std::cout << " File Checksum matches data. "<< std::endl;
  }
  else{
    std::cout << " Checksum broken, please check correlation functions ! "
       << std::endl;
    std::cout << "Read in: " << glob_check << " calculated: " << tmp << std::endl; 
  }

}
