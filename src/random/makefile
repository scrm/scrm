all: speedtest
	./speedtest

speedtest:
	g++ -Ofast -std=c++11 -o $@ speedtest.cc fastfunc.cc mersenne_twister.cc random_generator.cc
