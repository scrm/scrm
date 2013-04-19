bin_PROGRAMS = scrm scrm_dbg scrm_prof
PROG = SCRM
scrm_CXXFLAGS=  -O3 -Wall -DNDEBUG
scrm_dbg_CXXFLAGS= -O3 -Wall
scrm_prof_CXXFLAGS= -O3 -Wall -pg -DNDEBUG

scrm_SOURCES = \
param.cc param.h \
forest.cc forest.h forest-debug.cc \
model.cc model.h \
node.cc node.h \
node_container.cc node_container.h \
time_interval.cc time_interval.h \
tree_point.cc tree_point.h \
scrm.cc \
random/constant_generator.cc random/mersenne_twister.h  \
random/random_generator.cc   random/constant_generator.h \
random/mersenne_twister.cc   random/random_generator.h \
random/fastfunc.cc           random/fastfunc.h

scrm_dbg_SOURCES = \
param.cc param.h \
forest.cc forest.h forest-debug.cc \
model.cc model.h \
node.cc node.h \
node_container.cc node_container.h \
time_interval.cc time_interval.h \
tree_point.cc tree_point.h \
scrm.cc \
random/constant_generator.cc random/mersenne_twister.h  \
random/random_generator.cc   random/constant_generator.h \
random/mersenne_twister.cc   random/random_generator.h \
random/fastfunc.cc           random/fastfunc.h

scrm_prof_SOURCES = \
param.cc param.h \
forest.cc forest.h forest-debug.cc \
model.cc model.h \
node.cc node.h \
node_container.cc node_container.h \
time_interval.cc time_interval.h \
tree_point.cc tree_point.h \
scrm.cc \
random/constant_generator.cc random/mersenne_twister.h  \
random/random_generator.cc   random/constant_generator.h \
random/mersenne_twister.cc   random/random_generator.h \
random/fastfunc.cc           random/fastfunc.h

include_HEADERS = forest.h model.h node.h node_container.h time_interval.h tree_point.h \
	          random/constant_generator.h random/mersenne_twister.h random/random_generator.h
