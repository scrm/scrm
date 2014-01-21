SUBDIRS = src doc tests 
include aminclude.am 
EXTRA_DIST = oldMakefile bootstrap doc/manual.pdf acinclude.m4 tests/*sh tests/*r

#DISTFILES = $(DIST_COMMON) $(DIST_SOURCES) $(TEXINFOS) $(EXTRA_DIST)
distdir = $(PACKAGE)-$(VERSION)
CXXFLAGS= -I/opt/local/include
LDFLAGS= -L/opt/local/lib

unittests:
	cd tests; make unit_tests
	./tests/unit_tests 

algorithmtests:
	cd tests; make algorithm_tests
	./tests/algorithm_tests 

travis-test: unittests algorithm_tests
	cd src; make scrm scrm_dbg
	./tests/test_binaries.sh
