//g++ -D_BLOCK_SIZE_=16 -D_ARTIF_ -std=c++11 -ggdb -I/home/nagym/DevSdsl/include -L/home/nagym/DevSdsl/lib ./measurements.cc -o measurements -lsdsl -ldivsufsort -ldivsufsort64
//g++ -D_BLOCK_SIZE_ -D_ARTIF_ -std=c++11 -march=native -mavx -O3 -I/home/nagym/DevSdsl/include -L/home/nagym/DevSdsl/lib ./measurements.cc -o measurements -lsdsl -ldivsufsort -ldivsufsort64

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/bit_vectors.hpp>
#include <iostream>
#include <fstream>
#include <stdlib.h>

//#define BIG_PRIME 1294471
//#define BIG_PRIME 15485863
const uint64_t BIG_PRIME = 15485863; 

#define MAX_NUM_OF_RUNS 1
#define NUM_OF_QUERIES 1000

#define ACCESS  0x1
#define RANK    0x2
#define SELECT  0x3

#define OP_MODE SELECT

#define INPUT_PATH "/home/nagym/DevSdsl/measurements/input"
#define RESULT_PATH "/home/nagym/DevSdsl/measurements/results"

using namespace std;
using namespace sdsl;

// _BLOCK_SIZE_ must be given in compile time
// through g++ -DBLOCK_SIZE=X parameter !!!
// used 15, 16, 32, 48, 64, 96, 128, (160), 192, (224), 256

#ifndef _BLOCK_SIZE_
const int block_size = 3;//dummy
#else
const int block_size = _BLOCK_SIZE_;
#endif

// Input files are stored here
vector<string> inputFiles;
vector<string>::iterator fileIt;

__inline__ uint64_t rdtsc() {
  uint32_t low, high;
  __asm__ __volatile__ (
            "xorl %%eax,%%eax \n    cpuid"  //set 0 in eax and save info returned by cpuid
            ::: "%rax", "%rbx", "%rcx", "%rdx" );
  __asm__ __volatile__ (
            "rdtsc" : "=a" (low), "=d" (high));
  return (uint64_t)high << 32 | low;
}

string getOpModeString(){
  string opMode = "NOT_DEFINED";
  if (OP_MODE == ACCESS) opMode = "ACCESS";
  if (OP_MODE == RANK) opMode = "RANK";
  if (OP_MODE == SELECT) opMode = "SELECT";

  return opMode;
}

string stripPopulationFromFilename(const string& fileName){
    string inFile = fileName;
    size_t endpos   = fileName.size();
    size_t firstpos = fileName.find_last_of('_');
    firstpos++;

    string pop_str = inFile.substr(firstpos, endpos);
    pop_str.insert(1, ".");

    return pop_str;
}

double calculatePopulation(const bit_vector &bv){
    uint64_t numOfOnes = 0;
    for (int i=0; i< bv.size(); ++i){
      if (bv[i] == 1) numOfOnes++;
    }
    double pop = (double)numOfOnes / bv.size();
    return pop;
}

void generateRndmPositions(uint64_t posArr[], const int &length, const uint64_t &range){
  for (int i=0; i < length; ++i){
    posArr[i] = ((i+1)*BIG_PRIME) % range;
  }
}

void printArray(uint32_t arr[], const int &length){
  cout<<endl;
  cout<<" === Print Array ==="<<endl;
  for(int j=0; j < length; ++j){
    cout<<"["<<j<<"]: "<<arr[j]<<" ";
  }
  cout<<endl;
}

void measureCompressedBitmapSize(){
  string path;
  string resFile;
#ifdef _ARTIF_
  path    = (string)INPUT_PATH  + "/artificial/bitmaps/";
  resFile = (string)RESULT_PATH + "/artificial_av_SIZE/SIZE_artif.txt"; 
#elif _CORPUS_
  path = (string)INPUT_PATH + "/bitmap_corpus/";
  resFile = (string)RESULT_PATH + "/corpus_av_size/size_corpus.txt"; 
#endif
  vector<string>::iterator f_it;

  cout<<"========================================="<<endl;
  cout<<"      Starting SIZE measurement"<<endl;
  cout<<"========================================="<<endl;

  for (f_it=inputFiles.begin(); f_it < inputFiles.end(); ++f_it){
    string full_path_name = path + *f_it;
    ifstream file(full_path_name);
    if (!file.good()) {
      cout<<"File does not exist...exit."<<endl;
      exit(-1);
    }
    cout<<"Processing "<<full_path_name<<endl;

    bit_vector bv;
    load_from_file(bv, full_path_name);

    if (!bv.size()){
      cout<<"File seems to be empty...exit."<<endl;
      exit(-1);
    }

    cout<<"Size of file: "<<bv.bit_size()<<endl;

    rrr_vector<block_size>                      rrr(bv);
    rrri_vector<block_size>                     rrri(bv);
    r3d3_vector<block_size>                     r3d3(bv);
    r3d3i_vector<block_size>                    r3d3i(bv);
    ef_pure<>                                   ef_pure(bv);

    double size_rrr_Mb     = size_in_mega_bytes(rrr);
    double size_rrri_Mb    = size_in_mega_bytes(rrri);
    double size_r3d3_Mb    = size_in_mega_bytes(r3d3);
    double size_r3d3i_Mb   = size_in_mega_bytes(r3d3i);
    double size_ef_pure_Mb = size_in_mega_bytes(ef_pure);
    cout<<"  Size of RRR     in MB: " << size_rrr_Mb << endl;
    cout<<"  Size of RRRi    in MB: " << size_rrri_Mb << endl;
    cout<<"  Size of R3D3    in MB: " << size_r3d3_Mb << endl;
    cout<<"  Size of R3D3i   in MB: " << size_r3d3i_Mb << endl;
    cout<<"  Size of EF_PURE in MB: " << size_ef_pure_Mb << endl;

    //Obtain population (and entropy optionally)
#ifdef _ARTIF_
    string pop = stripPopulationFromFilename(*f_it);
    int N = bv.size();
    //Convert N to MByte
    double N_Mb = (double)N/(8*1024.0*1024.0);
    double popD = std::stod(pop); //C++11
#elif _CORPUS_
    double popD = calculatePopulation(bv);
    int N = bv.size();
    double N_Mb = (double)N/(8*1024.0*1024.0);
#endif
    // Entropy value
    double entropy = N_Mb*(popD*log2((double)1/popD) + (1-popD)*log2((double)1/(1-popD)));

    //rrr.printSizes();
    //rrri.printSizes();
    //r3d3.printSizes();
    //r3d3i.printSizes();

    // Writing result into a csv file (for R)
    ofstream outFile;
    outFile.open(resFile, ofstream::out | ofstream::app);

    outFile<<block_size<<" "<<popD<<" rrr "<<size_rrr_Mb<<" rrri "<<size_rrri_Mb<<
      " r3d3 "<<size_r3d3_Mb<<" r3d3i "<<size_r3d3i_Mb<<" ef_pure "<<size_ef_pure_Mb<<
      " en "<<entropy<<" "<<*f_it<<endl;
    //" en "<<entropy<<endl;

    outFile.close();
  }
}

void measureWTSize(){
    char* user;
    user = getenv ("USER");

    if (user == NULL)
       printf ("Error: USER variable is invalid: %s", user);

    string path = "/home/nagym";
    //path.append(user);
    path += "/ip_routing_src/EliasFano/src/input/text_corpus/";

    //path = "/tmp/";
    vector<string> files;
    vector<string>::iterator f_it;

    //files.push_back("html_source");
    //files.push_back("progc_source");
    files.push_back("shakespeare");
    files.push_back("scifi_book");
    files.push_back("bible");
    files.push_back("chr22_genome");
    files.push_back("chr7_genome");
    files.push_back("coli_genome");
    files.push_back("e.txt");
    files.push_back("pi_1m.txt");
    files.push_back("pi_10m.txt");

    //Set filenames according to operation mode
    const char operation[] = "results/size_wt";
    char fileName[100]; //to hold the result
    strcpy(fileName, operation);
    strcat(fileName, ".txt");

    cout<<"========================================="<<endl;
    cout<<"       Starting WT size measurement      "<<endl;
    cout<<"========================================="<<endl;

    for (f_it=files.begin(); f_it < files.end(); ++f_it){
      string full_path_name = path + *f_it;
      ifstream file(full_path_name);
      if (!file.good()) {
        cout<<"Input file does not exist...exit."<<endl;
        exit(-1);
      }
      cout<<"Processing "<<full_path_name<<endl;

      wt_huff< rrr_vector<block_size> >      wt_rrr;
      wt_huff< rrri_vector<block_size> >     wt_rrri;
      wt_huff< r3d3_vector<block_size> >     wt_r3d3;
      wt_huff< r3d3i_vector<block_size> >    wt_r3d3i;
      construct(wt_rrr, full_path_name, 1);
      construct(wt_rrri, full_path_name, 1);
      construct(wt_r3d3, full_path_name, 1);
      construct(wt_r3d3i, full_path_name, 1);

      cout<<"Size of loaded chars: "<<wt_rrr.size()<<endl;

      double size_wt_rrr_Mb   = size_in_mega_bytes(wt_rrr);
      double size_wt_rrri_Mb  = size_in_mega_bytes(wt_rrri);
      double size_wt_r3d3_Mb  = size_in_mega_bytes(wt_r3d3);
      double size_wt_r3d3i_Mb = size_in_mega_bytes(wt_r3d3i);
      cout<<"  Size of WT(RRR)	 in MB: " << size_wt_rrr_Mb << endl;
      cout<<"  Size of WT(RRRi)  in MB: " << size_wt_rrri_Mb << endl;
      cout<<"  Size of WT(R3D3)	 in MB: " << size_wt_r3d3_Mb << endl;
      cout<<"  Size of WT(R3D3i) in MB: " << size_wt_r3d3i_Mb << endl;

      // Entropy value
      //string pop = stripPopulationFromFilename(*f_it);
      //int N = bv.size();
      //Convert N to MByte
      //double N_Mb = (double)N/(8*1024.0*1024.0);
      //double popD = std::stod(pop); //C++11
      //double entropy = N_Mb*(popD*log2((double)1/popD) + (1-popD)*log2((double)1/(1-popD)));

      // Writing result into a txt file
      ofstream outFile;
      outFile.open(fileName, ofstream::out | ofstream::app);

      //method,t,block_size, population
      //replacing ',' with space
      //outFile<<block_size<<" "<<pop<<" or "<<size_or_Mb<<" or3 "<<size_or_3_Mb<<
      //        " ef "<<size_ef_Mb<<" ef3 "<<size_ef_3_Mb<<" en "<<entropy<<endl;

      outFile<<block_size<<" "<<" wt_rrr "<<size_wt_rrr_Mb<<" wt_rrri "<<size_wt_rrri_Mb<<
                 " wt_r3d3 "<<size_wt_r3d3_Mb<<" wt_r3d3i "<<size_wt_r3d3i_Mb<<" "<<*f_it<<endl;

      outFile.close();
    }
}

void measureBitmapOperation(){
  string path;
  string resFile;

  string opMode = getOpModeString();
#ifdef _ARTIF_
  path    = (string)INPUT_PATH  + "/artificial/bitmaps/";
  resFile = (string)RESULT_PATH + "/artificial_av_" + opMode + "/" + opMode + "_artif.txt"; 
#else
  path = (string)INPUT_PATH + "/bitmap_corpus/";
  resFile = (string)RESULT_PATH + "/corpus_av_" + opMode + "/" + opMode + "_corpus.txt"; 
#endif

  cout<<"========================================="<<endl;
  cout<<"   Starting "<<opMode<<" measurement  "<<endl;
  cout<<"   Number of queries: "<<NUM_OF_QUERIES<<endl;
  cout<<"   Output is written to: "<< resFile <<endl;
  cout<<"========================================="<<endl;

  for (fileIt=inputFiles.begin(); fileIt < inputFiles.end(); ++fileIt){
    string full_path_name = path + *fileIt;
    ifstream file(full_path_name);
    if (!file.good()) {
      cout<<"File does not exist...exit."<<endl;
      exit(-1);
    }
    cout<<"Processing "<<full_path_name<<endl;

    // Load input and calculate population (%)
    bit_vector bv;
    double popD;
#ifdef _ARTIF_
    load_from_file(bv, full_path_name);

    string pop = stripPopulationFromFilename(*fileIt);
    popD = std::stod(pop); //C++11
    unsigned int popcnt_bv = bv.size()*popD;
#else
    load_vector_from_file(bv, full_path_name, 0);

    popD = calculatePopulation(bv);
    unsigned int popcnt_bv = bv.size()*popD;
#endif
    //cout<<"# of bits: "<<bv.size()<<endl;
    //cout<<"# of bits: "<<bv.bit_size()<<endl;

    //Just to double check
    //cout<<" Rank/Select query bound: "<<popcnt_bv<<endl;

    // Building up Data Structures
    cout<<" Building up Data Structures ... ";
    rrr_vector<block_size>                    rrr(bv);
    rrri_vector<block_size>                   rrri(bv);
    r3d3_vector<block_size>                   r3d3(bv);
    r3d3i_vector<block_size>                  r3d3i(bv);
    rrr_vector<block_size>::rank_1_type       rrr_rank(&rrr);
    rrri_vector<block_size>::rank_1_type      rrri_rank(&rrri);
    r3d3_vector<block_size>::rank_1_type      r3d3_rank(&r3d3);
    r3d3i_vector<block_size>::rank_1_type     r3d3i_rank(&r3d3i);
    rrr_vector<block_size>::select_1_type     rrr_sel(&rrr);
    rrri_vector<block_size>::select_1_type    rrri_sel(&rrri);
    r3d3_vector<block_size>::select_1_type    r3d3_sel(&r3d3);
    r3d3i_vector<block_size>::select_1_type   r3d3i_sel(&r3d3i);
    cout << " Done. " <<endl;

    // Timers
    double t1_cpu_start, t1_cpu_stop;
    double t2_cpu_start, t2_cpu_stop;
    double t3_cpu_start, t3_cpu_stop;
    double t4_cpu_start, t4_cpu_stop;
    int k = 0;
    int i = 0;
    int pos = 0;
   
    // Generate random bit positions for acc/rank query
    cout<<" Generating random positions ... ";
    uint64_t rndmPosArrAccRank[NUM_OF_QUERIES];
    generateRndmPositions(rndmPosArrAccRank, NUM_OF_QUERIES, bv.size());
    //printArray(rndmPosArrAccRank, NUM_OF_QUERIES);

    // Generate random bit positions for select query
    uint64_t rndmPosArrSelect[NUM_OF_QUERIES];
    generateRndmPositions(rndmPosArrSelect, NUM_OF_QUERIES, popcnt_bv);
    //printArray(rndmPosArrSelect, NUM_OF_QUERIES);
    cout << " Done. " <<endl;

    // Measurements
    cout<<" Start measurements ... ";
    /******** RRR ********/
    uint64_t tmp=0;//for avoiding optimization
    t1_cpu_start = rdtsc();
    for (k=0; k < NUM_OF_QUERIES; ++k){
#if (OP_MODE == ACCESS)
      pos = rndmPosArrAccRank[k];
      tmp += rrr[pos];
#elif (OP_MODE == RANK)
      pos = rndmPosArrAccRank[k];
      tmp += rrr_rank.rank(pos);
#else
      pos = rndmPosArrSelect[k];
      tmp += rrr_sel.select(pos);
#endif
    }
    t1_cpu_stop = rdtsc();
    cout<<" "<<tmp;

    /******** RRRi ********/
    tmp=0; i = 0;
    t2_cpu_start = rdtsc();
    for (k=0; k < NUM_OF_QUERIES; ++k){
#if (OP_MODE == ACCESS)
      pos = rndmPosArrAccRank[k];
      tmp += rrri[pos];
#elif (OP_MODE == RANK)
      pos = rndmPosArrAccRank[k];
      tmp += rrri_rank.rank(pos);
#else
      pos = rndmPosArrSelect[k];
      tmp += rrri_sel.select(pos);
#endif
    }
    t2_cpu_stop = rdtsc();
    cout<<" "<<tmp;

    /******** R3D3 ********/
    tmp = 0; i = 0;
    t3_cpu_start = rdtsc();
    for (k=0; k < NUM_OF_QUERIES; ++k){
#if (OP_MODE == ACCESS)
      pos = rndmPosArrAccRank[k];
      tmp += r3d3[pos];    
#elif (OP_MODE == RANK)
      pos = rndmPosArrAccRank[k];
      tmp += r3d3_rank.rank(pos);
#else
      pos = rndmPosArrSelect[k];
      tmp += r3d3_sel.select(pos);    
#endif
    }
    t3_cpu_stop = rdtsc();
    cout<<" "<<tmp;

    /******** R3D3i ********/
    tmp = 0; i = 0;
    t4_cpu_start = rdtsc();
    for (k=0; k < NUM_OF_QUERIES; ++k){
#if (OP_MODE == ACCESS)
      pos = rndmPosArrAccRank[k];
      tmp += r3d3i[pos];
#elif (OP_MODE == RANK)
      pos = rndmPosArrAccRank[k];
      tmp += r3d3i_rank.rank(pos);
#else
      pos = rndmPosArrSelect[k];
      tmp += r3d3i_sel.select(pos);
#endif
    }
    t4_cpu_stop = rdtsc();
    cout<<" "<<tmp;
    cout<<" Done."<<endl;

    double delta_t1 = t1_cpu_stop - t1_cpu_start;
    double delta_t2 = t2_cpu_stop - t2_cpu_start;
    double delta_t3 = t3_cpu_stop - t3_cpu_start;
    double delta_t4 = t4_cpu_stop - t4_cpu_start;

    /*cout<<"File: "<<*f_it<<endl;
      cout<<"File size: "<<bv.size()<<" ("<<size_in_mega_bytes(bv)<<"MB)"<<" ; Block size: "<<block_size<<endl;
      cout<<"  OR  CPU time: "<<delta_t1<<endl;
      cout<<"  EF  CPU time: "<<delta_t2<<endl;
      cout<<"  EF2 CPU time: "<<delta_t3<<endl;
      cout<<endl;*/

    // Average operation time
    delta_t1 /= NUM_OF_QUERIES; 
    delta_t2 /= NUM_OF_QUERIES; 
    delta_t3 /= NUM_OF_QUERIES; 
    delta_t4 /= NUM_OF_QUERIES;

    // Writing result into file
    ofstream outFile;
    outFile.open(resFile, ofstream::out | ofstream::app);

    outFile<<block_size<<" "<<popD<<" rrr "<<delta_t1<<" rrri "<<delta_t2<<" r3d3 "<<
      delta_t3<<" r3d3i "<<delta_t4<<" "<<*fileIt<<" "<<rrr_rank.rank(bv.size())<<endl;

    outFile.close();
  }

}

//Todo refactor as above
void measureWTOperation(){
    //files to measure with different populations
    char* user;
    user = getenv ("USER");
    const int STEP = 1;

    if (user == NULL)
      printf ("Error: USER variable is invalid: %s", user);

    string path = "/home/nagym";
    path += "/ip_routing_src/EliasFano/src/input/text_corpus/";
    vector<string> files;
    vector<string>::iterator f_it;

    files.push_back("shakespeare");
    files.push_back("scifi_book");
    files.push_back("bible");
    files.push_back("chr22_genome");
    files.push_back("chr7_genome");
    files.push_back("coli_genome");
    files.push_back("e.txt");
    files.push_back("pi_1m.txt");
    files.push_back("pi_10m.txt");

    //Set filenames according to operation mode
    const char operation[] = "results/rank_wt";
    char fileName[100]; //to hold the result
    strcpy(fileName, operation);
    strcat(fileName, ".txt");//originally it was csv..

    cout<<"==================================================="<<endl;
    cout<<"   Starting "<<operation<<" measurement on WT      "<<endl;
    cout<<"==================================================="<<endl;

    for (f_it=files.begin(); f_it < files.end(); ++f_it){
        string full_path_name = path + *f_it;
        ifstream file(full_path_name);
        if (!file.good()) {
           cout<<"Input file does not exist...exit."<<endl;
           exit(-1);
        }
        cout<<"Processing "<<full_path_name<<endl;

        wt_huff< rrr_vector<block_size> >       wt_rrr;
        wt_huff< rrri_vector<block_size> >      wt_rrri;
        wt_huff< r3d3_vector<block_size> >      wt_r3d3;
        wt_huff< r3d3i_vector<block_size> >     wt_r3d3i;
        construct(wt_rrr,   full_path_name, 1);
        construct(wt_rrri,  full_path_name, 1);
        construct(wt_r3d3,  full_path_name, 1);
        construct(wt_r3d3i, full_path_name, 1);

        cout<<"# of chars: "<<wt_rrr.size()<<endl;

        double t1_cpu_start, t1_cpu_stop;
        double t2_cpu_start, t2_cpu_stop;
        double t3_cpu_start, t3_cpu_stop;
        double t4_cpu_start, t4_cpu_stop;
        int k = 0;

        /* Original RRR */
        //cout<<" WT - Original RRR "<<endl;
        char tmp = '\0';
        char symbol;
        if ( f_it->compare("e.txt") == 0 ||
             f_it->compare("pi_1m.txt") == 0 ||
             f_it->compare("pi_10m.txt") == 0 ){
          symbol = '5';
          cout<<"Symbol : 5"<<endl;
	}
        else{
	  symbol = 'a';
          cout<<"Symbol : a"<<endl;
        }
        cout<<"Quick test - rank(100): "<<wt_rrri.rank(100, symbol)<<endl;
        cout<<"Quick test - rank(100): "<<wt_r3d3i.rank(100, symbol)<<endl;
        t1_cpu_start = rdtsc();
        // NOTE wt_or.size() returns the length of the
        // original vector, eg. "aaa" -> 3
        for (k=0; k < wt_rrr.size(); k += STEP){
	  //tmp += wt_or[k];
          tmp += wt_rrr.rank(k, symbol);
        }
        t1_cpu_stop = rdtsc();
        cout<<"tmp: "<<tmp<<endl;

        /* RRR with relative indexes */
        //cout<<" WT - Indexed original RRR "<<endl;
        tmp = '\0';
        t2_cpu_start = rdtsc();
        for (k=0; k < wt_rrri.size() ; k += STEP){
	  //tmp += wt_or_3[k];
          tmp += wt_rrri.rank(k, symbol);
        }
        t2_cpu_stop = rdtsc();
        cout<<"tmp: "<<tmp<<endl;

        /* EliasFano'd RRR */
        //cout<<" WT - EliasFano RRR "<<endl;
        tmp = '\0';
        t3_cpu_start = rdtsc();
        for (k=0; k < wt_r3d3.size() ; k += STEP){
          //tmp += wt_ef[k];
          tmp += wt_r3d3.rank(k, symbol);
        }
        t3_cpu_stop = rdtsc();
        cout<<"tmp: "<<tmp<<endl;

        /* EliasFano'd RRR with relative indexes */
        //cout<<" WT - Indexed EliasFano RRR "<<endl;
        tmp = '\0';
        t4_cpu_start = rdtsc();
        for (k=0; k < wt_r3d3i.size() ; k += STEP){
          //tmp += wt_ef_3[k];
          tmp += wt_r3d3i.rank(k, symbol);
        }
        t4_cpu_stop = rdtsc();
        cout<<"tmp: "<<tmp<<endl;

        double delta_t1 = t1_cpu_stop - t1_cpu_start;
        double delta_t2 = t2_cpu_stop - t2_cpu_start;
        double delta_t3 = t3_cpu_stop - t3_cpu_start;
        double delta_t4 = t4_cpu_stop - t4_cpu_start;

        /*cout<<"File: "<<*f_it<<endl;
        cout<<"File size: "<<bv.size()<<" ("<<size_in_mega_bytes(bv)<<"MB)"<<" ; Block size: "<<block_size<<endl;
        cout<<"  OR  CPU time: "<<delta_t1<<endl;
        cout<<"  EF  CPU time: "<<delta_t2<<endl;
        cout<<"  EF2 CPU time: "<<delta_t3<<endl;
        cout<<endl;*/

        // Average operation time
        delta_t1 /= wt_rrr.size();
        delta_t2 /= wt_rrr.size();
        delta_t3 /= wt_rrr.size();
        delta_t4 /= wt_rrr.size();//miert nem bv.size()??
        //cout<<"wt_or size: "<<wt_or.size()<<endl;
        //cout<<"wt_ef size: "<<wt_ef.size()<<endl;
        //cout<<"wt_or_3 size: "<<wt_or_3.size()<<endl;
        //cout<<"wt_ef_3 size: "<<wt_ef_3.size()<<endl;

        // Writing result into a csv file (for R)
        ofstream outFile;
        outFile.open(fileName, ofstream::out | ofstream::app);

        outFile<<block_size<<" "<<" wt_rrr "<<delta_t1<<" wt_rrri "<<delta_t2<<" wt_r3d3 "<<
                delta_t3<<" wt_r3d3i "<<delta_t4<<" "<<*f_it<<endl;

        outFile.close();
      }
}

void measureBlockDecoderAccessSpeed(){
    //timers
    double t1_cpu_start, t1_cpu_stop;
    double t2_cpu_start, t2_cpu_stop;
    double delta_t1, delta_t2;

    string path = (string)INPUT_PATH + "/artificial/bitmaps/";
    string fileName = "bitmap_1000000_01"; 

    //vector<string> files;
    //vector<string>::iterator f_it;
    //files.push_back("bitmap_1000000_01");

    bit_vector bv;
    string full_path_name = path + fileName;
    ifstream file(full_path_name);
    if (!file.good()) {
      cout<<"File does not exist...exit."<<endl;
      exit(-1);
    }

    // Prepare rndm positions for query
    uint64_t rndmPosArrAcc[block_size];
    generateRndmPositions(rndmPosArrAcc, block_size, block_size);
    //printArray(rndmPosArrAcc, block_size);

    //load_from_file(bv, full_path_name);
    load_vector_from_file(bv, full_path_name, 0);
    string pop = stripPopulationFromFilename(fileName);

    //helper typedefs
    typedef rrr_helper<block_size> rrr_helper_type;
    typedef r3d3_helper<block_size> r3d3_helper_type;
    typedef rrr_helper_type::number_type number_type;

    uint64_t sum_t1 = 0;
    uint64_t sum_t2 = 0;
    uint64_t numOfBlocks  = 0;
    uint64_t sum_size_rrr  = 0;
    uint64_t sum_size_ef = 0;
    for (int pos = 0; pos < bv.size(); pos+=block_size){
      numOfBlocks++;
      //cout<<"pos: "<<pos<<endl;
      //*** OBTAIN BLOCK SIZE LONG CHUNK FROM BITVECTOR *** //
      // -> BIN
      number_type bin = rrr_helper_type::trait::get_int(bv, pos, block_size);
      int bt = rrr_helper_type::trait::popcount(bin);

      //*** ENCODE BIN TO RRR offset ***//
      number_type nr = rrr_helper_type::bin_to_nr(bin);
      uint16_t space_for_bt_rrr = rrr_helper_type::space_for_bt(bt);
      sum_size_rrr += space_for_bt_rrr;

      //*** ENCODE BIN TO EF offset ***//
      uint16_t msb_bit_pos = r3d3_helper_type::hi(bin);
      uint16_t space_for_bt_r3d3 = r3d3_helper_type::space_for_bt_i(msb_bit_pos, bt);
      bit_vector m_btnr(space_for_bt_r3d3);
      uint16_t l = (bt > 0 ? log2(msb_bit_pos/bt) : 0);
      r3d3_helper_type::compress_ef(m_btnr, 0, bin, l);
      sum_size_ef += space_for_bt_r3d3;

      //*** MEASUREMENT ***//
        int num = 0;
        int bit_pos = 0;
        //RRR access each bit position
        t1_cpu_start = rdtsc();
        for (int i=0; i < block_size; ++i){
          bit_pos = rndmPosArrAcc[i];
          num += rrr_helper_type::decode_bit(bt, nr, bit_pos);
        }
        t1_cpu_stop = rdtsc();
        delta_t1 = t1_cpu_stop-t1_cpu_start;
        delta_t1 /= block_size;
        //cout<<"CPU time: "<<delta_t1<<endl;
        sum_t1 += delta_t1;
        //cout<<"Popcnt of block: "<<num<<endl;

        //EliasFano access each bit position
        num = 0;
        t2_cpu_start = rdtsc();
        for (int i=0; i < block_size; ++i){
          bit_pos = rndmPosArrAcc[i];
          num += r3d3_helper_type::decode_bit_ef(m_btnr, 0, space_for_bt_r3d3, l, bt, bit_pos);
        }
        t2_cpu_stop = rdtsc();
        delta_t2 = t2_cpu_stop-t2_cpu_start;
        delta_t2 /= block_size;
        //cout<<"CPU time: "<<delta_t2<<endl;
        sum_t2 += delta_t2;
        //cout<<"Popcnt of block: "<<num<<endl;

    }
    //*** WRITE RESULTS TO FILE ***
    //cout<<"sum RRR access: "<<sum_t1<<endl;
    //cout<<"sum EF  access: "<<sum_t2<<endl;

    cout<<"Num of blocks: "<<numOfBlocks<<endl;
    double av_access_rrr  = (double)sum_t1/(double)numOfBlocks;
    double av_access_ef   = (double)sum_t2/(double)numOfBlocks;
    cout<<"average RRR access: "<<av_access_rrr<<endl;
    cout<<"average EF  access: "<<av_access_ef<<endl;

    double av_size_rrr  = (double)sum_size_rrr/(double)numOfBlocks;
    double av_size_ef   = (double)sum_size_ef/(double)numOfBlocks;
    cout<<"average RRR  size : "<<av_size_rrr<<endl;
    cout<<"average R3D3 size : "<<av_size_ef<<endl;

    ofstream outFile;
    outFile.open("results/artificial_block_decoder_ACCESS_compare/block_decoder_compare.txt", 
        ofstream::out | ofstream::app);

    outFile<<block_size<<" "<<pop<<" or "<<av_access_rrr<<" ef "<<av_access_ef
            <<" or_size "<<av_size_rrr<<" ef_size "<<av_size_ef<<endl;
    outFile.close();
}

void measureRRRSizeAndOperationPerBlock(){
  const int block_size = 256;
  const char operation[] = "rank";
  char fileName[20]; //to hold the result
  strcpy(fileName, operation);
  strcat(fileName, ".data");

  //create 256 bit long bitvector
  bit_vector bv = bit_vector(block_size, 0); //starting from 0
  double t1_cpu_start, t1_cpu_stop;
  double t2_cpu_start, t2_cpu_stop;
  int k = 0;

  //start to fill it with 1s, and do it in the worst
  //scenario from EF's point of view -> starting from
  //bit position 255
  for (int i = 0; i <= 64; ++i){
  //for (int i = 255; i >= 192; i--){
    bv[i] = 1;

    rrri_vector<block_size>                   rrri(bv);
    rrri_vector<block_size>::rank_1_type      rrri_rank(&rrri);
    rrri_vector<block_size>::select_1_type    rrri_sel(&rrri);
    r3d3i_vector<block_size>                  r3d3i(bv);
    r3d3i_vector<block_size>::rank_1_type     r3d3i_rank(&r3d3i);
    r3d3i_vector<block_size>::select_1_type   r3d3i_sel(&r3d3i);

    // MEASURE COMPRESSED SIZE
    double size_rrri_Mb    = size_in_mega_bytes(rrri);
    double size_r3d3i_Mb   = size_in_mega_bytes(r3d3i);
    double size_rrri_Kb    = size_rrri_Mb*1024;
    double size_r3d3i_Kb   = size_r3d3i_Mb*1024;

    // MEASURE OPERATION SPEED
    // - Original RRR
    t1_cpu_start = rdtsc();
    for (k=0; k < block_size; ++k){
       //int tmp = rrr_or_3[k];
       int tmp = rrri_rank.rank(k);
       //int tmp = rrr_or_3_sel.select(k);
    }
    t1_cpu_stop = rdtsc();

    // - EliasFano RRR
    t2_cpu_start = rdtsc();
    for (k=0; k < block_size; ++k){
       //int tmp = rrr_ef_3[k];
       int tmp = r3d3i_rank.rank(k);
       //int tmp = rrr_ef_3_sel.select(k);
    }
    t2_cpu_stop = rdtsc();

    double delta_t1 = t1_cpu_stop - t1_cpu_start;
    double delta_t2 = t2_cpu_stop - t2_cpu_start;

    // cutting down e+06
    delta_t1 /= 1000000;//00;
    delta_t2 /= 1000000;//00;

    // COLLECTING THE RESULTS
    cout<<" Bitvector length: "<<block_size<<", # of 1s: "<<i<<endl;
    cout<<" /Size/"<<endl;
    cout<<"   OR3: "<<size_rrri_Kb<<endl;
    cout<<"   EF3: "<<size_r3d3i_Kb<<endl;
    cout<<" /Speed/"<<endl;
    cout<<"   OR3 CPU time: "<<delta_t1<<endl;
    cout<<"   EF3 CPU time: "<<delta_t2<<endl;
    rrri.printSizes();
    r3d3i.printSizes();

    // WRITING INTO FILE
    ofstream dataFile;
    dataFile.open(fileName, ofstream::out | ofstream::app);

    dataFile<<"bs "<<block_size<<" pop "<<i<<" or3_speed "<<delta_t1<<" ef3_speed "<<delta_t2<<
           " or3_size "<<size_rrri_Kb<<" ef3_size "<<size_r3d3i_Kb<<endl;
    dataFile.close();
  }
}

void measureFIBs(){
   string path = "/home/nagym/ip_routing_src/EliasFano/src/input/fibs_bin/fibs/";

   //Same names in all directories
   string s_int = "S_int.bmp.bin_h";
   string s_alpha = "S_alpha2.txt";

   vector<string> dirs;
   vector<string>::iterator d_it;

   dirs.push_back("/hbone_szeged_2014_06_01_00_15_31.anon.txt/");
   dirs.push_back("/SFR_HMS_internet.txt/");
   dirs.push_back("/SFR_LNS_local.txt/");
   dirs.push_back("/orange.txt/");
   dirs.push_back("/hbone_vh1_2014_06_01_00_23_36.anon.txt/");

   for (d_it = dirs.begin(); d_it != dirs.end(); ++d_it){
     //Check if files exist
     string full_name_bm = path + *d_it + s_int;
     string full_name_txt = path + *d_it + s_alpha;
     ifstream fileBm(full_name_bm);
     ifstream fileTxt(full_name_txt);
     if ( !fileBm.good() || !fileTxt.good() ) {
       cout<<"File "<<full_name_bm<<" or "<<endl;
       cout<<"File "<<full_name_txt<<" does not exist...exit."<<endl;
       exit(-1);
     }
     else {
       cout<<"file s_int: "<<*d_it + s_int<<endl;
       cout<<"file s_alp: "<<*d_it + s_alpha<<endl;
     }

     //read S_int
     bit_vector bv;
     load_from_file(bv, full_name_bm);
     rrri_vector<block_size>                   rrri(bv);
     r3d3i_vector<block_size>                  r3d3i(bv);
     rrri_vector<block_size>::rank_1_type      rrri_rank(&rrri);
     r3d3i_vector<block_size>::rank_1_type     r3d3i_rank(&r3d3i);

     //measure s_int compressed size
     double size_s_int_rrri_Mb = size_in_mega_bytes(rrri);
     double size_s_int_r3d3i_Mb = size_in_mega_bytes(r3d3i);

     //read S_alpha
     wt_huff< rrri_vector<block_size> >        wt_rrri;
     wt_huff< r3d3i_vector<block_size> >       wt_r3d3i;
     construct(wt_rrri, full_name_txt, 1);
     construct(wt_r3d3i, full_name_txt, 1);

     //measure S_alpha compressed size
     double size_s_alpha_rrri_Mb = size_in_mega_bytes(wt_rrri);
     double size_s_alpha_r3d3i_Mb = size_in_mega_bytes(wt_r3d3i);

     /****************************************************************
      * X ACCESS on S_int:
      *
      */
     double t1_cpu_start, t1_cpu_stop;
     double t2_cpu_start, t2_cpu_stop;
     int tmp, i = 0;

     //generate rndm indexes offline
     // index[0-4]: for access
     // index[5] : for rank
     int x = 5;
     int index[x+1];
     for (int i = 0; i < x + 1; ++i){
         index[i] = ((i+1)*BIG_PRIME) % bv.size();
         //cout<<"index["<<i<<"]: "<<index[i]<<" ; total size: "<<bv.size()<<endl;
     }

     //RRRi
     t1_cpu_start = rdtsc();
     for (int k = 0; k < x; ++k){
        tmp += rrri[index[k]];
     }
     t1_cpu_stop = rdtsc();
     cout<<"tmp: "<<tmp<<endl;

     //R3D3i
     tmp = 0;
     t2_cpu_start = rdtsc();
     for (int k = 0; k < x; ++k){
        tmp += r3d3i[index[k]];
     }
     t2_cpu_stop = rdtsc();
     cout<<"tmp: "<<tmp<<endl;

     double delta_t1 = t1_cpu_stop - t1_cpu_start;
     double delta_t2 = t2_cpu_stop - t2_cpu_start;
     //delta_t1 /= x;
     //delta_t2 /= x;
//     cout<<"OR3 access time: "<<delta_t1<<endl;
//     cout<<"EF3 access time: "<<delta_t2<<endl;

     /****************************************************************
      * 1 RANK on S_int:
      *
      */
     double t3_cpu_start, t3_cpu_stop;
     double t4_cpu_start, t4_cpu_stop;
     //RRR
     tmp = 0;
     t3_cpu_start = rdtsc();
     tmp += rrri_rank.rank(index[5]);
     t3_cpu_stop = rdtsc();
     cout<<"tmp: "<<tmp<<endl;

     //R3D3i
     tmp = 0;
     t4_cpu_start = rdtsc();
     tmp+= r3d3i_rank.rank(index[5]);
     t4_cpu_stop = rdtsc();
     cout<<"tmp: "<<tmp<<endl;

     double delta_t3 = t3_cpu_stop - t3_cpu_start;
     double delta_t4 = t4_cpu_stop - t4_cpu_start;
//     cout<<"OR3 rank time: "<<delta_t3<<endl;
//     cout<<"EF3 rank time: "<<delta_t4<<endl;

     /****************************************************************
      * 1 WT access on S_alpha:
      *
      */
     int index_alpha = 3*BIG_PRIME % wt_rrri.size();
     //cout<<" index for WT access: "<<index_alpha<<" ; WT size: "<<wt_or_3.size()<<endl;
     tmp = 0;
     double t5_cpu_start, t5_cpu_stop;
     double t6_cpu_start, t6_cpu_stop;

     cout<<"WT RRRi  size: "<<wt_rrri.size()<<endl;
     cout<<"WT R3D3i size: "<<wt_r3d3i.size()<<endl;

     t5_cpu_start = rdtsc();
     tmp += wt_rrri[index_alpha];
     t5_cpu_stop = rdtsc();
     cout<<"tmp: "<<tmp<<endl;

     tmp = 0;
     t6_cpu_start = rdtsc();
     tmp += wt_r3d3i[index_alpha];
     t6_cpu_stop = rdtsc();
     cout<<"tmp: "<<tmp<<endl;

     double delta_t5 = t5_cpu_stop - t5_cpu_start;
     double delta_t6 = t6_cpu_stop - t6_cpu_start;
//     cout<<"WT OR3 access time: "<<delta_t5<<endl;
//     cout<<"WT EF3 access time: "<<delta_t6<<endl;

//     cout<<"========================================"<<endl;
//     cout<<"block_size: "<<block_size<<endl;
//     cout<<" Size OR3 S_int: "<<size_s_int_or_Mb<<endl;
//     cout<<" Size OR3 S_alpha: "<<size_s_alpha_or_Mb<<endl;
//     cout<<" Size EF3 S_int : "<<size_s_int_ef_Mb<<endl;
//     cout<<" Size EF3 S_alpha : "<<size_s_alpha_ef_Mb<<endl;
//     cout<<endl;
//     cout<<" Sum OR3: "<<delta_t1+delta_t3+delta_t5<<endl;
//     cout<<" Sum EF3: "<<delta_t2+delta_t4+delta_t6<<endl;
//     cout<<"========================================"<<endl;

     ofstream outFile;
     outFile.open("results/fib.txt", ofstream::out | ofstream::app);

     //method,t,block_size, population
     //using spaces instead of ','
     outFile<<block_size<<" "<<" or3_s_i "<<size_s_int_rrri_Mb<<" or3_s_a "<<
             size_s_alpha_rrri_Mb<<" ef3_s_i "<<size_s_int_r3d3i_Mb<<" ef_3_s_a "<<
             size_s_alpha_r3d3i_Mb<<" or3 "<<delta_t1+delta_t3+delta_t5<<" ef3 "<<
             delta_t2+delta_t4+delta_t6<<" "<<*d_it<<endl;

     outFile.close();

   }
}

int main() {

  if (block_size == 3){
    cout<<" When compiling measurements.cc, give the block size"<<endl;
    cout<<" in the GCC parameters like this: -DBLOCK_SIZE=32 ."<<endl;
    cout<<" Exiting..."<<endl;
    exit(-1);
  }

#ifdef _ARTIF_
  inputFiles.push_back("bitmap_500000000_000001");
  inputFiles.push_back("bitmap_500000000_000005");
  inputFiles.push_back("bitmap_500000000_00001");
  inputFiles.push_back("bitmap_500000000_00005");
  inputFiles.push_back("bitmap_500000000_0001");
  inputFiles.push_back("bitmap_500000000_0005");
  inputFiles.push_back("bitmap_500000000_001");
  inputFiles.push_back("bitmap_500000000_005");
  inputFiles.push_back("bitmap_500000000_01");
  inputFiles.push_back("bitmap_500000000_015");
  inputFiles.push_back("bitmap_500000000_025");
  inputFiles.push_back("bitmap_500000000_035");
  inputFiles.push_back("bitmap_500000000_05");
  //inputFiles.push_back("bitmap_1000000_01");
  //inputFiles.push_back("bitmap_500000000_01");

  /* Measurements on artificial bitmaps */

  // Fig 5.-6.
  //measureBlockDecoderAccessSpeed();

  //measureCompressedBitmapSize();

  measureBitmapOperation();

#elif _CORPUS_
  inputFiles.push_back("fax_calg.bin_h");
  inputFiles.push_back("reg-1.bmp.bin_h");
  inputFiles.push_back("phd.bmp.bin_h");
  inputFiles.push_back("zip.bmp.bin_h");
  inputFiles.push_back("caida_4.bmp.bin_h");
  inputFiles.push_back("caida_8.bmp.bin_h");
  inputFiles.push_back("caida_16.bmp.bin_h");
#elif _TEXT_
  inputFiles.push_back("shakespeare");
  inputFiles.push_back("scifi_book");
  inputFiles.push_back("bible");
  inputFiles.push_back("chr22_genome");
  inputFiles.push_back("chr7_genome");
  inputFiles.push_back("coli_genome");
  inputFiles.push_back("e.txt");
  inputFiles.push_back("pi_1m.txt");
  inputFiles.push_back("pi_10m.txt");
#else
  cout<<" There is no mode defined (_ARTIF_, _CORPUS_, _TEXT_), exit..."<<endl;
  exit(-1);
#endif


  /* Measure size */
  //measureCompressedBitmapSize();
  //measureBitmapOperation();

  /* Measurements on Wavelet Trees */
  //measureWTSize();
  //measureWTOperation();

  //fix the block size and check performance against population
  //measureSizeAndOperationPerBlock();

  // Fig 5.-6.
  //measureBlockDecoderAccessSpeed();

  //measureFIBs();
}
