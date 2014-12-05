#ifndef IO_HELPERS_H_
#define IO_HELPERS_H_

#include <array>
#include <complex>
#include <cstdlib>
#include <vector>

#include "boost/crc.hpp"
#include "lime.h"
#include "typedefs.h"
///////////////////////////////////////////////////////////////////////////////
// Typedefs ///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

struct pdg{
  size_t id;
  std::array<int,3> p;
  std::array<int,3> dis;
  std::array<int,4> gamma;
};

struct Tag {
  int mom[2];
  int dis[2];
  int gam[2];
};

struct GlobalDat {
  std::vector<size_t> rnd_seeds;
  size_t nb_rnd_vecs;
  size_t nb_perambs;
};
 void write_1st_msg(const char* filename, GlobalDat& dat,
                          boost::uint64_t chksum);
 void append_msgs(const char* filename, std::vector<vec>& corr, std::vector<Tag>& tags,
              LimeWriter* w, bool be);
   
// Endianess ///////////////////////////////////////////////////////////////////
// test system endianess
inline int big_endian () {
	union {
		int l;
		char c[sizeof(int)];
	} u;

	u.l = 1;
	return (u.c[sizeof(int) - 1] == 1);
}

// swapping endianess 
template <typename T> inline T swap_endian(T u) {
  union {
    T u;
    unsigned char u8[sizeof(T)];
  } source, dest;
  source.u = u;
  for (size_t k = 0; k < sizeof(T); k++)
    dest.u8[k] = source.u8[sizeof(T) - k - 1]; 
  return dest.u;
}

inline std::complex<double>  swap_complex(std::complex<double> val){

  return std::complex<double>(swap_endian<double>(std::real(val)),
                              swap_endian<double>(std::imag(val)));
}

// swap endianess of one correlation function
inline std::vector<cmplx> swap_single_corr(const std::vector<cmplx>& corr){
  size_t ext = corr.size();
  // Temporary vector same size as input
  std::vector<cmplx> le(ext);
  // Swap endianness for every entry in corr
  for(size_t val = 0; val < ext; ++val) {
    le[val] = swap_complex(corr[val]);
  }
  return le;
}


// swap endianess of one tag
inline Tag swap_single_tag(const Tag& tag){
  Tag le_tag;
  for(size_t pos = 0; pos < 2; ++pos ){
    le_tag.mom[pos] = swap_endian<int>(tag.mom[pos]);
    le_tag.dis[pos] = swap_endian<int>(tag.dis[pos]);
    le_tag.gam[pos] = swap_endian<int>(tag.gam[pos]);
  }
  return le_tag;
}

// swap endaness of runinfo
inline GlobalDat swap_glob_dat(const GlobalDat& run_info){
  GlobalDat le_glob;
  for (auto seed = 0; seed < run_info.rnd_seeds.size(); ++seed)
    le_glob.rnd_seeds.push_back( swap_endian<size_t>(run_info.rnd_seeds[seed]) );
  le_glob.nb_rnd_vecs = swap_endian<size_t>(run_info.nb_rnd_vecs);
  le_glob.nb_perambs = swap_endian<size_t>(run_info.nb_perambs);
  return le_glob;
}

// type independent checksum of container
template <typename MyContainer> 
inline boost::uint64_t checksum(const MyContainer& dat, size_t bytes){
  boost::crc_32_type chksum_agent;
  chksum_agent.process_bytes(&dat[0], bytes);
  return chksum_agent.checksum();
}

// check existence of a file
inline bool file_exist(const char* name) {
  if (FILE* file = fopen(name, "r")) {
    fclose(file);
    return true;
  } 
  else {
    return false;
  }   
}
// Set the tag from two operator structures
 void set_tag(Tag& so_si, const pdg& op1, const pdg& op2);
// Compare two tags of correlation functions
 bool compare_tags(const Tag& tag1, const Tag& tag2);
// Calculate p^2
 int square_comp(const std::array<int, 3>& p1,
                       const std::array<int, 3>& p2);

// swap vector of all correlation functions
inline void swap_correlators(std::vector<vec>& corr){
  for (auto& func : corr) func = swap_single_corr(func);
}
// swap vector of all tags
inline void swap_tag_vector(std::vector<Tag>& tags){
  for (auto& label : tags) label = swap_single_tag(label);
}

#endif // IO_HELPERS_H_
