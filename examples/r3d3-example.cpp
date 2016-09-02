// usr/bin/c++ -std=c++11 -Wall -Wextra -DNDEBUG  -O3 -ffast-math -funroll-loops -msse4.2 \
// -DHAVE_CXA_DEMANGLE -I/home/nagym/GitHub/sdsl-lite/include -o r3d3-example r3d3-example.cpp \
// -lsdsl -ldivsufsort -ldivsufsort64

#include <sdsl/bit_vectors.hpp>
#include <iostream>

using namespace std;
using namespace sdsl;

int main()
{
  bit_vector b = bit_vector(63, 0); //starting from 0
  b[5] = 1;
  b[17] = 1;
  b[31] = 1;
  b[48] = 1;

  r3d3_vector<32> rb(b);
  r3d3_vector<32>::rank_1_type rank_rb(&rb);
  r3d3_vector<32>::select_1_type select_rb(&rb);  

  // access
  cout << rb[5] << endl;
  cout << rb[6] << endl;

  // rank
  cout << rank_rb(17) << endl;
  cout << rank_rb(18) << endl;

  // select
  cout << select_rb(2) <<endl;
  cout << select_rb(4) <<endl;

  r3d3i_vector<32> rbi(b);
  r3d3i_vector<32>::rank_1_type rank_rbi(&rbi);
  r3d3i_vector<32>::select_1_type select_rbi(&rbi);  

  // access
  cout << rbi[5] << endl;
  cout << rbi[6] << endl;

  // rank
  cout << rank_rbi(17) << endl;
  cout << rank_rbi(18) << endl;

  // select
  cout << select_rbi(2) <<endl;
  cout << select_rbi(4) <<endl;
}
