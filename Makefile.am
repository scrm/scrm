EXTRA_DIST = doc/create-manual.sh doc/knitr.css doc/scrm.1
bin_PROGRAMS = scrm
man_MANS = doc/scrm.1

TESTS = unit_tests algorithm_tests
check_PROGRAMS = unit_tests algorithm_tests scrm_dbg scrm_asan scrm_prof
PROG = SCRM

dist-hook:
	chmod u+w $(distdir)/doc
	cd $(distdir); ./doc/create-manual.sh $(VERSION)

test: unittests

unittests: unit_tests
	./unit_tests

algorithmtest: algorithm_tests
	./algorithm_tests

nrml_src = src/param.cc src/forest.cc src/node.cc src/node_container.cc src/time_interval.cc \
		   src/model.cc src/tree_point.cc \
		   src/param.h src/forest.h src/node.h src/node_container.h src/time_interval.h \
		   src/model.h src/tree_point.h src/event.h src/contemporaries_container.h \
		   src/macros.h

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
			  src/summary_statistics/oriented_forest.cc \
			  src/summary_statistics/oriented_forest.h

scrm_src = $(nrml_src) $(random_src) $(sumstat_src)

debug_src = src/random/constant_generator.cc src/random/constant_generator.h \
			src/forest-debug.cc src/random/constant_generator.h

unit_test_src = tests/unittests/test_forest.cc tests/unittests/test_model.cc\
				tests/unittests/test_node.cc tests/unittests/test_node_container.cc\
				tests/cppunit/test_runner.cc tests/unittests/test_time_interval.cc\
				tests/unittests/test_fastfunc.cc tests/unittests/test_param.cc\
				tests/unittests/test_random_generator.cc tests/unittests/test_summary_statistics.cc\
				tests/unittests/test_contemporaries_container.cc

alg_test_src = tests/cppunit/test_runner.cc tests/algorithmtest/test_algorithm.cc


scrm_SOURCES = $(scrm_src) src/scrm.cc
scrm_dbg_SOURCES = $(scrm_src) $(debug_src) src/scrm.cc
scrm_prof_SOURCES = $(scrm_src) src/scrm.cc
scrm_asan_SOURCES = $(scrm_src) src/scrm.cc
unit_tests_SOURCES = $(scrm_src) $(debug_src) $(unit_test_src)
algorithm_tests_SOURCES = $(scrm_src) $(alg_test_src)

scrm_CXXFLAGS= -DNDEBUG @OPT_CXXFLAGS@
scrm_dbg_CXXFLAGS= -g
scrm_prof_CXXFLAGS= -pg -DNDEBUG
scrm_asan_CXXFLAGS= -g -DNDEBUG -fsanitize=undefined,address -fno-sanitize-recover
unit_tests_CXXFLAGS = -g -DUNITTEST -DNDEBUG @TEST_CXXFLAGS@ 
algorithm_tests_CXXFLAGS = -g -DNDEBUG
unit_tests_LDADD= -L/opt/local/lib -lcppunit -ldl #link the cppunit unittest library in mac, cppunit was installed via macports
algorithm_tests_LDADD= -L/opt/local/lib -lcppunit -ldl  #link the cppunit unittest library in mac, cppunit was installed via macports
