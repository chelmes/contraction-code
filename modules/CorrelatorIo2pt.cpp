#include "CorrelatorIo2pt.h"
///////////////////////////////////////////////////////////////////////////////
// Write 2pt correlator ///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void write_2pt_lime(const char* filename, GlobalDat& dat, std::vector<Tag>& tags,
                    std::vector<vec>& corr){
  bool be = big_endian();
  if (be){
    swap_correlators(corr);
    swap_tag_vector(tags);
    dat = swap_glob_dat(dat);
  }
  // if the file should not exist initialize it with config info as first message
  if(!file_exist(filename)){
    // calculate checksum of all correlators
    size_t glob_bytes = corr.size() * (corr[0]).size();
    boost::uint64_t global_chksum = checksum <std::vector<vec> > (corr,
                                                                  glob_bytes);
    std::cout << "Global Checksum is: " << global_chksum << std::endl;
    if(be) global_chksum = swap_endian<boost::uint64_t>(global_chksum);
    write_1st_msg(filename, dat, global_chksum);
  }
  // setup what is needed for output to lime
  FILE* fp;
  LimeWriter* w;
  fp = fopen( filename, "a" );
  w = limeCreateWriter( fp );
  // writing the correlators to the end
  append_msgs(filename, corr, tags, w, be);
  limeDestroyWriter( w );
  fclose(fp);
}
///////////////////////////////////////////////////////////////////////////////
// Read 2pt correlator ////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


// IDEA: First read in everything and compare checksums, if not correct give
// warning and exit
// afterwards look for correlator with tag.

void read_2pt_lime(const char* filename, std::vector<Tag>& tags, std::vector<vec> correlators){
  if(tags.size() == correlators.size()){
    if(FILE* fp = fopen(filename, r)){
      std::vector<boost::uint64_t>checksums(tags.size(),0);
      int MB_flag = 0; int ME_flag = 0;
      int file_status = 0;
      int act_status = 0;
      limeReader r = limeCreateReader(fp);
      n_uint64_t tag_bytes = sizeof(Tag);
      n_uint64_t data_bytes = correlators[0].size()*2*sizeof(double);
      while (file_status != LIME_EOF){
      // read tag
      file_status = limeReaderNextRecord(r);
      MB_flag = limeReaderMBFlag(r);
      ME_flag = limeReaderMEFlag(r);
      std::cout << MB_flag << " " << ME_flag << std::endl;
      // get Metainfo of record in read tag
      if(MB_flag == 0 && ME_flag == 0){
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
    else std::cout << "File does not exist!" << std::endl;
  }
  else std::cout << "#elements for tags and correlators not equal" << std::endl;
}

void get_2pt_lime(const char* filename, const Tag& tag,
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
    if(MB_flag == 0 && ME_flag == 0){
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

