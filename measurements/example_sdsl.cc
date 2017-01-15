//g++ -std=c++11 -ggdb -I/home/nagym/DevSdsl/include -L/home/nagym/DevSdsl/lib ./example_sdsl.cc -o example -lsdsl -ldivsufsort -ldivsufsort64

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/uint128_t.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/suffix_trees.hpp>
#include <iostream>
#include <fstream>
//#include "immintrin.h" _lzcnt_
//#include "cstdint"

using namespace std;
using namespace sdsl;

typedef uint128_t number_type;

inline int clz_u128 (uint128_t u) {
  uint64_t hi = u>>64;
  uint64_t lo = u;
  int retval[3]={
    __builtin_clzll(hi),
    __builtin_clzll(lo)+64,
    128
  };
  int idx = !hi + ((!lo)&(!hi));
  return retval[idx];
}


int main() {
  //bit_vector b = bit_vector(80*(1<<20), 0);

  bit_vector b = bit_vector(64, 0); //starting from 0
  //b[0] = 1;
  //b[2] = 1;
  //b[3] = 1;
  b[5] = 1;
  b[10] = 1;
  //b[11] = 1;
  b[17] = 1;
  b[31] = 1;
  b[48] = 1;
  
  ef_pure<> ef(b);
  ef.printCompressedData();

}
