CC=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=$(shell find src | grep .cc)
OBJECTS=$(SOURCES:.cc=.o)
TEST_SOURCES=$(shell find tests | grep .cc)
TEST_OBJECTS=$(TEST_SOURCES:.cc=.to) $(SOURCES:.cc=.to)
REBUILD=src/scrm
EXECUTABLE=scrm

all: $(SOURCES) $(EXECUTABLE)
	  
$(EXECUTABLE): $(OBJECTS) 
	  $(CC) $(LDFLAGS) $(OBJECTS) -o $@

tmp:	
	echo $(TEST_OBJECTS)

test: unittest
	./unittest

unittest: $(TEST_OBJECTS)
	  $(CC)  $(TEST_OBJECTS) -o $@ $(LDFLAGS) `cppunit-config --libs`

%.o:  %.cc  %.h
	  $(CC) $(CFLAGS) $< -o $@

%.to:  %.cc %.h
	  $(CC) $(CFLAGS) `cppunit-config --cflags` -DUNITTEST $< -o $@
