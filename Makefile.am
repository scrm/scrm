SUBDIRS = src doc tests 
include aminclude.am 
EXTRA_DIST = oldMakefile bootstrap doc/manual.pdf acinclude.m4 tests/*sh tests/*r
#DISTFILES = $(DIST_COMMON) $(DIST_SOURCES) $(TEXINFOS) $(EXTRA_DIST)
distdir = $(PACKAGE)-$(VERSION)
CXXFLAGS= -I/opt/local/include
LDFLAGS= -L/opt/local/lib

travis-test: 
	cd src; make scrm
	cd tests; make check 
	./tests/test_binaries.sh
