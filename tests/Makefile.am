TESTS = tests

check_PROGRAMS = tests

tests_SOURCES = test_forest.cc test_model.cc test_node.cc test_node_container.cc\
			   	test_runner.cc test_time_interval.cc\
			   	../src/forest.h ../src/model.h ../src/node.h ../src/node_container.h\
			   	../src/time_interval.h ../src/tree_point.h ../src/forest.cc ../src/model.cc\
			   	../src/node.cc ../src/node_container.cc ../src/time_interval.cc ../src/tree_point.cc\
				../src/forest-debug.cc\
			   	../src/random/constant_generator.cc ../src/random/constant_generator.h\
			   	../src/random/random_generator.cc ../src/random/random_generator.h\
			   	../src/random/mersenne_twister.cc ../src/random/mersenne_twister.h\
			   	../src/random/fastfunc.cc ../src/random/fastfunc.h  

#AM_CPPFLAGS= -I/opt/local/include # include the library for boost, given that boost was installed through macports
tests_CXXFLAGS = -DUNITTEST -DNDEBUG
#tests_LDADD= -lcppunit
#tests_LDFLAGS=-L/opt/local/lib
tests_LDADD= -lcppunit -ldl #link the cppunit unittest library in mac, cppunit was installed via macports

include_HEADERS = ../src/forest.h ../src/model.h ../src/node.h ../src/node_container.h ../src/time_interval.h ../src/tree_point.h
