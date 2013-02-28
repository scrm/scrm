bin_PROGRAMS = scrm
PROG = SCRM
#CXXFLAGS=-c -Wall -O3 -pg -DNDEBUG  	# For gprof
#CXXFLAGS= -Wall -O3 -DNDEBUG 		# For speed testing
CXXFLAGS=  -O3 -Wall -pg 	-DNDEBUG		# For debugging
#scrm_CXXFLAGS=  -O3 -Wall -pg 	-DNDEBUG		# For debugging
#scrm_LDFLAGS=-pg #-pg flag is for debug

scrm_SOURCES = \ 
forest.cc forest.h forest-debug.cc \ 
model.cc model.h \ 
node.cc node.h \ 
node_container.cc node_container.h \
time_interval.cc time_interval.h \ 
tree_point.cc tree_point.h \ 
scrm.cc \
random/constant_generator.cc random/mersenne_twister.h  \
random/random_generator.cc  random/constant_generator.h \
random/mersenne_twister.cc   random/random_generator.h



include_HEADERS = forest.h model.h node.h node_container.h time_interval.h tree_point.h



