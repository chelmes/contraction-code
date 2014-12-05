#include "CorrelatorIo2pt.h"
///////////////////////////////////////////////////////////////////////////////
//Helper Routines//////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//Endianess////////////////////////////////////////////////////////////////////
// test endianess
static int big_endian () {
	union {
		int l;
		char c[sizeof(int)];
	} u;

	u.l = 1;
	return (u.c[sizeof(int) - 1] == 1);
}

// swapping endianess ----------------------------------------------------------

template <typename T> T swap_endian(T u) {
  union {
    T u;
    unsigned char u8[sizeof(T)];
  } source, dest;
  source.u = u;
  for (size_t k = 0; k < sizeof(T); k++)
    dest.u8[k] = source.u8[sizeof(T) - k - 1]; 
  return dest.u;
}

//
inline std::complex<double>  swap_complex(std::complex<double> val){

  return std::complex<double>(swap_endian<double>(std::real(val)),
                              swap_endian<double>(std::imag(val)));
}

// swap endianess of one correlation function
static std::vector<cmplx> swap_single_corr(const std::vector<cmplx>& corr){
  size_t ext = corr.size();
  // Temporary vector same size as input
  std::vector<cmplx> le(ext);
  // Swap endianness for every entry in corr
  for(size_t val = 0; val < ext; ++val) {
    le[val] = swap_complex(corr[val]);
  }
  return le;
}

// use for all correlation functions
inline void swap_correlators(std::vector<vec>& corr){
  for (auto& func : corr) func = swap_single_corr(func);
}

// swap endianess of one tag
static Tag swap_single_tag(const Tag& tag){
  Tag le_tag;
  for(size_t pos = 0; pos < 2; ++pos ){
    le_tag.mom[pos] = swap_endian<int>(tag.mom[pos]);
    le_tag.dis[pos] = swap_endian<int>(tag.dis[pos]);
    le_tag.gam[pos] = swap_endian<int>(tag.gam[pos]);
  }
  return le_tag;
}

// now for the whole vector
inline void swap_tag_vector(std::vector<Tag>& tags){
  for (auto& label : tags) label = swap_single_tag(label);
}

// swap endianess of global data

static GlobalDat swap_glob_dat(const GlobalDat& run_info){
  GlobalDat le_glob;
  for (auto seed = 0; seed < run_info.rnd_seeds.size(); ++seed)
    le_glob.rnd_seeds.push_back( swap_endian<size_t>(run_info.rnd_seeds[seed]) );
  le_glob.nb_rnd_vecs = swap_endian<size_t>(run_info.nb_rnd_vecs);
  le_glob.nb_perambs = swap_endian<size_t>(run_info.nb_perambs);
  return le_glob;
}

//Tag handling/////////////////////////////////////////////////////////////////

// Calculate p^2
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

//File handling////////////////////////////////////////////////////////////////
//check existence of a file
static bool file_exist(const char* name) {
  if (FILE* file = fopen(name, "r")) {
    fclose(file);
    return true;
  } 
  else {
    return false;
  }   
}
template <typename MyContainer> 
static boost::uint64_t checksum(const MyContainer& dat, size_t bytes){
  boost::crc_32_type chksum_agent;
  chksum_agent.process_bytes(&dat[0], bytes);
  return chksum_agent.checksum();
}

static void write_1st_msg(const char* filename, GlobalDat& dat,
                          boost::uint64_t chksum){
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

static void append_msgs(const char* filename, std::vector<vec>& corr, std::vector<Tag>& tags,
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
///////////////////////////////////////////////////////////////////////////////
//Write 2pt correlator/////////////////////////////////////////////////////////
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

