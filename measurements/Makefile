# ############################################ #
#  Make binaries for all block sizes           #
# ############################################ #

CC=g++

#CXXFLAGS=-std=c++0x -g -ggdb -O0 -Wall -DDEBUG_LOG -I./

SDSLDIR = $(HOME)/DevSdsl
CXXFLAGS = -std=c++0x -march=native -mavx -O3 -I$(SDSLDIR)/include
LINKOBJ = -lsdsl -ldivsufsort -ldivsufsort64
TARGET = measure_16 measure_32 measure_64 measure_128 measure_256

SRCS=measurements.cc

# three type of measurements: 
#    - artificial, 
#    - bitmap corpus, 
#    - texts
ifndef TS
 $(error TESTSET is not set. Use TS=[_ARTIF_][_CORPUS_][_TEXT_])
endif

ifeq ($(TS),_ARTIF_)
 TESTFLAG = -D_ARTIF_
else ifeq ($(TS),_CORPUS_)
 TESTFLAG = -D_CORPUS_
else ifeq ($(TS),_TEXT_)
 TESTFLAG = -D_TEXT_
else 
 $(error Wrong TESTSET is given)
endif

all: $(TARGET)

measure_16: measurements.cc
	$(CC) -D_BLOCK_SIZE_=16 $(TESTFLAG) $(CXXFLAGS) measurements.cc -o $@ $(LINKOBJ)

measure_32: measurements.cc
	$(CC) -D_BLOCK_SIZE_=32 $(TESTFLAG) $(CXXFLAGS) measurements.cc -o $@ $(LINKOBJ)

measure_64: measurements.cc
	$(CC) -D_BLOCK_SIZE_=64 $(TESTFLAG) $(CXXFLAGS) measurements.cc -o $@ $(LINKOBJ)

measure_128: measurements.cc
	$(CC) -D_BLOCK_SIZE_=128 $(TESTFLAG) $(CXXFLAGS) measurements.cc -o $@ $(LINKOBJ)

measure_256: measurements.cc
	$(CC) -D_BLOCK_SIZE_=256 $(TESTFLAG) $(CXXFLAGS) measurements.cc -o $@ $(LINKOBJ)


clean:
	rm $(TARGET)
