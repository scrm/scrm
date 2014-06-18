include aminclude.am 
EXTRA_DIST = doc/manual.html m4/ac_doxygen.m4
bin_PROGRAMS = scrm 

TESTS = unit_tests algorithm_tests 
check_PROGRAMS = unit_tests algorithm_tests scrm_dbg scrm_prof
PROG = SCRM

dist-hook:
	export VERSION=${VERSION}; cd doc; make

unittests: unit_tests
	./unit_tests 

algorithmtests: algorithm_tests
	./algorithm_tests 

travis-test: unittests algorithmtests scrm scrm_dbg
	./tests/test_binaries.sh

nrml_src = src/param.cc src/forest.cc src/node.cc src/node_container.cc src/time_interval.cc \
		   src/model.cc src/tree_point.cc \
		   src/param.h src/forest.h src/node.h src/node_container.h src/time_interval.h \
		   src/model.h src/tree_point.h src/event.h

random_src = src/random/random_generator.cc src/random/mersenne_twister.cc \
			 src/random/fastfunc.cc \
			 src/random/random_generator.h src/random/mersenne_twister.h \
			 src/random/fastfunc.h 

sumstat_src = src/summary_statistics/tmrca.cc \
			  src/summary_statistics/seg_sites.cc \
			  src/summary_statistics/frequency_spectrum.cc \
			  src/summary_statistics/newick_tree.cc \
			  src/summary_statistics/tmrca.h \
			  src/summary_statistics/seg_sites.h \
			  src/summary_statistics/frequency_spectrum.h \
			  src/summary_statistics/newick_tree.h \
			  src/summary_statistics/summary_statistic.h \
			  src/summary_statistics/first_last_tmrca.cc \
			  src/summary_statistics/first_last_tmrca.h
 
scrm_src = $(nrml_src) $(random_src) $(sumstat_src)

debug_src = src/random/constant_generator.cc src/random/constant_generator.h \
			src/forest-debug.cc src/random/constant_generator.h

unit_test_src = tests/test_forest.cc tests/test_model.cc tests/test_node.cc tests/test_node_container.cc\
			   	tests/test_runner.cc tests/test_time_interval.cc tests/test_fastfunc.cc tests/test_param.cc\
				tests/test_random_generator.cc tests/test_summary_statistics.cc 

alg_test_src = tests/test_runner.cc tests/test_algorithm.cc 


scrm_SOURCES = $(scrm_src) src/scrm.cc
scrm_dbg_SOURCES = $(scrm_src) $(debug_src) src/scrm.cc
scrm_prof_SOURCES = $(scrm_src) src/scrm.cc
unit_tests_SOURCES = $(scrm_src) $(debug_src) $(unit_test_src)
algorithm_tests_SOURCES = $(scrm_src) $(alg_test_src)
 
scrm_CXXFLAGS= -DNDEBUG -std=c++0x
scrm_dbg_CXXFLAGS= -std=c++0x
scrm_prof_CXXFLAGS= -pg -DNDEBUG -std=c++0x
unit_tests_CXXFLAGS = -DUNITTEST -DNDEBUG -std=c++0x 
algorithm_tests_CXXFLAGS = -DUNITTEST -DNDEBUG -std=c++0x 
unit_tests_LDADD= -L/opt/local/lib -lcppunit -ldl #link the cppunit unittest library in mac, cppunit was installed via macports
algorithm_tests_LDADD= -L/opt/local/lib -lcppunit -ldl  #link the cppunit unittest library in mac, cppunit was installed via macports
