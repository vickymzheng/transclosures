CXX=g++
CXXFLAGS=-std=c++14 -O2 -Wall -lboost_program_options

all: stream_cc

clean:
	rm -rf stream_cc