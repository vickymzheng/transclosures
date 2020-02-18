all: stream_cc nano-tools-stream

stream_cc: stream_cc.cpp
	g++ -std=c++14 -O2 -Wall -lboost_program_options stream_cc.cpp -o stream_cc

nano-tools-stream: nano-tools-stream.cpp
	g++ -std=c++14 -O2 -Wall nano-tools-stream.cpp -o nano-tools-stream

clean:
	rm -rf stream_cc
	rm -rf nano-tools-stream 