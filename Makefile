CC=g++
#CFLAGS=-c -Wall -O3 -pg -DNDEBUG 	# For gprof
#CFLAGS=-c -Wall -O3 -DNDEBUG 		# For speed testing
CFLAGS=-c -O3 -Wall -pg 			# For debugging
LDFLAGS=-pg
SOURCES=$(shell find src | grep .cc)
OBJECTS=$(SOURCES:.cc=.o)
TEST_SOURCES=$(shell find tests | grep .cc)
TEST_OBJECTS=$(TEST_SOURCES:.cc=.to) $(SOURCES:.cc=.to)
REBUILD=src/scrm
EXECUTABLE=scrm

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) Makefile
	  $(CC) $(LDFLAGS) $(OBJECTS) -o $@

run: $(EXECUTABLE)
	./scrm

test: unittest
	./unittest

unittest: $(TEST_OBJECTS) Makefile
	  $(CC)  $(TEST_OBJECTS) -o $@ $(LDFLAGS) `cppunit-config --libs`

%.o:  %.cc
	  $(CC) $(CFLAGS) $< -o $@

%.to:  %.cc
	  $(CC) $(CFLAGS) `cppunit-config --cflags` -DNDEBUG -DUNITTEST $< -o $@
	  #$(CC) $(CFLAGS) `cppunit-config --cflags` -DUNITTEST $< -o $@

%.cc: %.h
	  touch $@

clean: 
	rm $(OBJECTS) $(TEST_OBJECTS)
