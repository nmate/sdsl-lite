R3D3 - Fast Compressed Bitvectors for the Succinct Data Structure
=========

What is it?
-----------

R3D3 is a succinct data structure prototype that attains several times
space reduction beyond known compression techniques while it supports
operations on the compressed data at comparable speed. Our prototype is
built on top of the Succinct Data Structure Library (SDSL) [sdsl-lite][SDSLLIT].

Why R3D3?
--------
The volume of datasets waiting for being processed increases drastically nowadays.
Data mining, machine learning, pattern analysis and networking are the main
areas facing this challenge. The succinct data structures have very appealing 
properties, namely they not just compress information but they also make possible
to perform queries on the compressed data. These are the two key elements that
are needed for fast data processing. Compression makes possible to better use the
cache and memory, waiving the painful cost of disk accesses while enhanced operations
assure the efficient data processing on the other hand. 

A good example for such an opportunistic compression approach is the RRR compressed 
bitvector scheme due to Raman, Raman, and Rao attaining nHo bits on the data while
supporting access, rank and select queries in optimal O(1) time. A major shortcoming 
of compressed information processing is, however, that the storage size of the index 
can significantly outweigh (up to and beyond 8 times) that of the data.

To address this limitation R3D3 is made as a doubly opportunistic data structure
which, as opposed to conventional opportunistic schemes that compress only the data 
component, achieve information-theoretically minimal entropy-constrained space both 
on the data and the index at the same time. 

R3D3 is implemented in two versions. First, r3d3_vector is constructed with limited 
indexes, making it more compact at the price of sacrificing query performance. In 
opposite r3d3i_vector is `fully` indexed, enabling fast operation speed but with
an increased memory footprint.

Installation
------------

Since R3D3 is an extension to [sdsl-lite][SDSLLIT], please refer to it for installation
details.

Getting Started
------------

Using R3D3 is just like using a regular sdsl object. The following snippet
shows how to initialize R3D3 with a simple `b` bitvector and how to execute
access, rank and select operations in practice.

```cpp
#include <sdsl/bit_vectors.hpp>
#include <iostream>

using namespace sdsl;
int main(){
 bit_vector b = bit_vector(63, 1);
 r3d3_vector<32> rb(b);
 r3d3_vector<32>::rank_1_type rank_rb(&rb);
 r3d3_vector<32>::select_1_type select_rb(&rb);  
 std::cout << rb[5] << std::endl;
 std::cout << rank_rb(17) << std::endl;
 cout << select_rb(2) <<endl;
}
```
To compile the program using `g++` run:
```sh
usr/bin/c++ -std=c++11 -Wall -DNDEBUG -O3 -I{SDSL_DIR}/include -L{SDSL_DIR}/lib -o program program.cpp -lsdsl -ldivsufsort -ldivsufsort64
```

Test
----

R3D3 tests are inserted to [sdsl-lite][SDSLLIT]'s unit tests. All R3D3 related 
source code and testcase can be found with `r3d3` prefix. Running only R3D3 testcases
can be done by 

```sh
cd sdsl-lite/build
make r3d3-bit-vector-test
make r3d3-rank-support-test
make r3d3-select-support-test
```

Bug Reporting
------------
Mail <will be published later>

Licensing
---------
<will be detailed later>

Authors
-------
<will be detailed later>

[SDSLLIT]: https://github.com/simongog/sdsl-lite/wiki/Literature "Succinct Data Structure Literature"