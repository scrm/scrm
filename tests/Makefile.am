TESTS = tests

check_PROGRAMS = tests

tests_SOURCES = test_forest.cc test_model.cc test_node.cc test_node_container.cc\
			   	test_runner.cc test_time_interval.cc test_fastfunc.cc test_param.cc\
				test_random_generator.cc test_seg.cc test_algorithm.cc\
			   	../src/forest.h ../src/model.h ../src/node.h ../src/node_container.h\
			   	../src/param.h ../src/param.cc ../src/seg.h ../src/seg.cc\
			   	../src/time_interval.h ../src/tree_point.h ../src/forest.cc ../src/model.cc\
			   	../src/node.cc ../src/node_container.cc ../src/time_interval.cc ../src/tree_point.cc\
				../src/forest-debug.cc\
			   	../src/random/constant_generator.cc ../src/random/constant_generator.h\
			   	../src/random/random_generator.cc ../src/random/random_generator.h\
			   	../src/random/mersenne_twister.cc ../src/random/mersenne_twister.h\
			   	../src/random/fastfunc.cc ../src/random/fastfunc.h  

#AM_CPPFLAGS= -I/opt/local/include # include the library for boost, given that boost was installed through macports
tests_CXXFLAGS = -DUNITTEST -DNDEBUG -std=c++0x
#tests_CXXFLAGS = -DUNITTEST -std=c++0x
#tests_LDADD= -lcppunit
#tests_LDFLAGS=-L/opt/local/lib
tests_LDADD= -lcppunit -ldl #link the cppunit unittest library in mac, cppunit was installed via macports

include_HEADERS = ../src/forest.h ../src/model.h ../src/node.h ../src/node_container.h ../src/time_interval.h\
				  ../src/tree_point.h ../src/param.h ../src/seg.h
