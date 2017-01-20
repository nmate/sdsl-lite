//g++ -std=c++11 -ggdb -I/home/nagym/DevSdsl/include -L/home/nagym/DevSdsl/lib ./generate_bitstring.cc -o generate_bitstring -lsdsl -ldivsufsort -ldivsufsort64

#include <iostream>
#include <map>
#include <random>
#include <iterator>
#include <stdio.h>
#include <fstream>
#include <sstream>  
#include <string>
#include <sdsl/bit_vectors.hpp>

#define N16      16
#define N32      32
#define N64      64
#define N128    128
#define N256    256
#define N1000  1000
#define N1M    100000000

#define PROB 0.1

using namespace std;
using namespace sdsl;

string int2str(const int& number){
  ostringstream convert;  // stream used for the conversion
  convert << number;  // insert the textual representation of 'Number'

  return convert.str();
}

string prob2str(const double &prob){
  ostringstream convert;
  convert << prob;
  string tmp = convert.str();
  size_t pointPos = tmp.find_last_of('.');
  pointPos++;
  size_t endPos = tmp.length();
  string toRet = tmp.substr(pointPos, endPos);

  return toRet;
}

//generates text file
void generateString_OwnWay(const int &n, const double &p, const string &dir){
  string fileName =  dir + "/bitmap_" + int2str(n) + "_0" + prob2str(p);
  ofstream myFile(fileName, ios::out);
  int numOfOnes = 0;

  std::random_device rd;
  std::mt19937 gen(rd());
  double weights[] =
          {1,     // generate '0' with base probability
           p};    // number 1 at p probability
  std::discrete_distribution<> d(std::begin(weights), std::end(weights));

  if ( myFile.is_open() ){
    for (int i=0; i<n; ++i){
      int value = d(gen);
      //cout<<value;
      myFile<<value;
      if (value == 1) numOfOnes++;
    }
    myFile<<endl;
    cout<<endl;
  }
  else {
    cerr<<"Cannot open the file!"<<endl;
  }

  cout<<"Num of '1's: "<<numOfOnes<<endl;
  cout<<"Num of '0's: "<<n-numOfOnes<<endl;
}

//generates binary file that can be read by sdsl
void generateString_SdslWay(const int &n, const double &p, const string &dir){
  //file to write
  string fileName =  dir + "/bitmap_" + int2str(n) + "_0" + prob2str(p);

  //create bitmap
  bit_vector b = bit_vector(n, 0);

  //fill it with appropriate probability
  int numOfOnes = 0;
  std::random_device rd;
  std::mt19937 gen(rd());
  double weights[] =
      {1-p,     // generate '0' with base probability
       p};    // number 1 at p probability
  std::discrete_distribution<> d(std::begin(weights), std::end(weights));

  for (int i=0; i<b.size(); ++i){
    int value = d(gen);
    if (value == 1) {
      b[i] = value;
      numOfOnes++;
    }
  }

  //store the bitvector to file
  store_to_file(b, fileName);

  cout<<"File: "<<fileName<<endl;
  cout<<"Num of '1's: "<<numOfOnes<<endl;
  cout<<"Num of '0's: "<<n-numOfOnes<<endl;

  numOfOnes = 0;
  for (int i=0; i< b.size(); ++i){
    if (b[i] == 1) numOfOnes++;
  }
  double pop = (double)numOfOnes / b.size();
  cout<<"Calculated population: "<<pop<<endl;
}

int main(int argc, char **argv){
  string prob;
  string dir;
  string size;
  if (argc == 4){
    size = argv[1];
    prob = argv[2];
    dir  = argv[3];
  }
  else{
    cout<<"Missing [size of bitmap] [probability] and [target dir] arguments!"<<endl;
    cout<<" Try this: generate_bitstring 1000 0.5 1"<<endl;
    exit(-1);
  }

  string::size_type sz;  
  double probD = stod(prob, &sz);
  int N = stoi(size, &sz);

  cout<<"Called with probability: "<<probD<<endl; 
  cout<<"Target directory: "<<dir<<endl;
  cout<<"Bitmap's size: "<<N<<endl;
  cout<<"Population: "<<probD<<endl;
  //generateString(N16,  probD, dir);
  //generateString(N32,  probD, dir);
  //generateString(N64,  probD, dir);
  //generateString(N128, probD, dir);
  //generateString(N256, probD, dir);
  //generateString_OwnWay(N, probD, dir);
  generateString_SdslWay(N, probD,dir);


}
