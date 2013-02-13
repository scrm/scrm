CC=g++
#CFLAGS=-c -Wall -O3 -pg -DNDEBUG 	# For gprof
#CFLAGS=-c -Wall -O3 -DNDEBUG 		# For speed testing
CFLAGS=-c -O3 -Wall -pg 			# For debugging
LDFLAGS=-pg
SOURCES=$(shell find src | grep .cc)
HEADERS=$(shell find src | grep .h)
OBJECTS=$(SOURCES:.cc=.o)
TEST_SOURCES=$(shell find tests | grep .cc)
TEST_OBJECTS=$(TEST_SOURCES:.cc=.to) $(SOURCES:.cc=.to)
REBUILD=src/scrm
EXECUTABLE=scrm

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	  $(CC) $(LDFLAGS) $(OBJECTS) -o $@

run: $(EXECUTABLE)
	./scrm

test: unittest
	./unittest

unittest: $(TEST_OBJECTS) 
	  $(CC)  $(TEST_OBJECTS) -o $@ $(LDFLAGS) `cppunit-config --libs`

%.o:  %.cc $(HEADERS)
	  $(CC) $(CFLAGS) $< -o $@

%.to:  %.cc $(HEADERS)
	  $(CC) $(CFLAGS) `cppunit-config --cflags` -DNDEBUG -DUNITTEST $< -o $@
	  #$(CC) $(CFLAGS) `cppunit-config --cflags` -DUNITTEST $< -o $@

clean: 
	- rm $(OBJECTS) $(TEST_OBJECTS) 2> /dev/null
