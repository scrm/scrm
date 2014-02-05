bin_PROGRAMS = scrm scrm_dbg scrm_prof
PROG = SCRM

scrm_CXXFLAGS= -DNDEBUG -std=c++0x
scrm_dbg_CXXFLAGS= -std=c++0x
scrm_prof_CXXFLAGS= -pg -DNDEBUG -fno-inline -std=c++0x

normal_src = \
param.cc 					param.h \
forest.cc 					forest.h \
model.cc 					model.h \
node.cc 					node.h \
node_container.cc 			node_container.h \
time_interval.cc 			time_interval.h \
tree_point.cc 				tree_point.h \
seg.cc 						seg.h \
event.h \
scrm.cc \
random/random_generator.cc random/random_generator.h  \
random/mersenne_twister.cc random/mersenne_twister.h  \
random/fastfunc.cc         random/fastfunc.h \
summary_statistics/summary_statistic.h
 
debug_src = random/constant_generator.cc random/constant_generator.h \
			forest-debug.cc 

scrm_SOURCES = $(normal_src) 

scrm_dbg_SOURCES = $(normal_src) $(debug_src)

scrm_prof_SOURCES = $(normal_src)

include_HEADERS = seg.h forest.h model.h node.h node_container.h time_interval.h tree_point.h \
	          random/mersenne_twister.h random/random_generator.h \
			  summary_statistics/summary_statistic.h
