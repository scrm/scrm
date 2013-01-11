CC=g++
#CFLAGS=-c -Wall -O3 -pg -DNDEBUG
CFLAGS=-c -Wall -O3 -pg
LDFLAGS=-pg
SOURCES=$(shell find src | grep .cc)
OBJECTS=$(SOURCES:.cc=.o)
TEST_SOURCES=$(shell find tests | grep .cc)
TEST_OBJECTS=$(TEST_SOURCES:.cc=.to) $(SOURCES:.cc=.to)
REBUILD=src/scrm
EXECUTABLE=scrm

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	  $(CC) $(LDFLAGS) $(OBJECTS) -o $@

test: unittest
	./unittest

unittest: $(TEST_OBJECTS)
	  $(CC)  $(TEST_OBJECTS) -o $@ $(LDFLAGS) `cppunit-config --libs`

%.o:  %.cc
	  $(CC) $(CFLAGS) $< -o $@

%.to:  %.cc
	  #$(CC) $(CFLAGS) `cppunit-config --cflags` -DNDEBUG -DUNITTEST $< -o $@
	  $(CC) $(CFLAGS) `cppunit-config --cflags` -DUNITTEST $< -o $@

%.cc: %.h
	  touch $@

clean: 
	rm $(OBJECTS) $(TEST_OBJECTS)
