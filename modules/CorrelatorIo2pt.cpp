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

      //concatenate all correlation functions in one vector
      std::vector<cmplx> collect(glob_bytes);
      size_t length = glob_bytes*2*sizeof(double);
      std::cout << "lngth is: " << length << std::endl;
      for(auto& c : corr)
        for (auto& el : c) collect.push_back(el);
      size_t global_chksum = checksum <std::vector<cmplx> > (collect,
                                                                    length);
      std::cout << "Global Checksum is: " << global_chksum << std::endl;
      if(be) global_chksum = swap_endian<size_t>(global_chksum);
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

  void read_2pt_lime(const char* filename, std::vector<Tag>& tags, std::vector<vec>& correlators){
    bool bigend = big_endian();
    if(tags.size() == correlators.size()){
      if(FILE* fp = fopen(filename, "r")){
        size_t global_check = 0;
        std::vector<size_t> checksums(tags.size());
        
        int MB_flag = 0; int ME_flag = 0;
        int file_status = 0;
        int act_status = 0;
        LimeReader* r = limeCreateReader(fp);
        n_uint64_t check_bytes = sizeof(size_t);
        n_uint64_t tag_bytes = sizeof(Tag);
        n_uint64_t data_bytes = correlators[0].size()*2*sizeof(double);
        file_status = limeReaderNextRecord(r);
        // From first message read global checksum
        MB_flag = limeReaderMBFlag(r);
        ME_flag = limeReaderMEFlag(r);
        std::cout << MB_flag << " " << ME_flag << std::endl;
        if(MB_flag == 1 && ME_flag == 0){
          act_status = limeReaderReadData(&global_check, &check_bytes, r);
          std::cout << global_check << std::endl;
          if (bigend) global_check = swap_endian<size_t>(global_check);
          std::cout << global_check << std::endl;
      }
      //file_status = limeReaderNextRecord(r);
      // TODO: think about read in runinfo
      file_status = limeReaderNextRecord(r);
      // loop over all remaining messages
     int cnt = 0;
     do {
      //std::cout << "message: " << cnt/3 << std::endl;
      file_status = limeReaderNextRecord(r);
      MB_flag = limeReaderMBFlag(r);
      ME_flag = limeReaderMEFlag(r);
      //std::cout << MB_flag << " " << ME_flag << std::endl;
      // read correlator checksum
      if(MB_flag == 1 && ME_flag == 0){
        size_t check = 0;
        act_status = limeReaderReadData(&check, &check_bytes, r);
        if (bigend) check = swap_endian<size_t>(check);
        checksums.at(cnt/3) = check;
        //std::cout << check << std::endl;
      }
      // read correlator tag
      else if(MB_flag == 0 && ME_flag == 0){
        Tag read_tag;
        act_status = limeReaderReadData(&read_tag, &tag_bytes, r);
        if (bigend) read_tag = swap_single_tag(read_tag);
        tags.at(cnt/3) = read_tag;
      }
      // read correlator data
      else if(MB_flag == 0 && ME_flag == 1){
        vec corr;
        corr.resize(96);
        //std::cout << "initialized tmeporal correlator" << std::endl;
        act_status = limeReaderReadData(corr.data(), &data_bytes, r);
        //std::cout << "read_in tmeporal correlator" << std::endl;
        //std::cout << corr.size() << std::endl;
        if (bigend) corr = swap_single_corr(corr);
        correlators.at(cnt/3) = corr;
      }
      ++cnt;
      } while (file_status != LIME_EOF);
      // Check file integrity
     //std::cout << correlators.size();
     size_t glob_bytes = correlators.size()*correlators[0].size();
    //concatenate all correlation functions in one vector
    std::vector<cmplx> collect(glob_bytes);
    for(auto& c : correlators)
      for (auto& el : c) collect.push_back(el);
      file_check(global_check, checksums, collect);
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

