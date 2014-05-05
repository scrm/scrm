#!/usr/bin/env python
import ld_test as ld

import sys

if __name__ == "__main__":
    _use_param = ld.read_param_file ( sys.argv[1] )
    _use_param.printing()
    
    processed_data = []
    runtime = []
    for job in _use_param.jobs:
        print job
        f = open ( job+"tmrcaRho", "r" )        
        processed_data.append( [float(x) for x in f.readline().strip('[').strip(']\n').split(',') ] )
        f.close()
        f = open ( job+"time", "r" )
        runtime.append( float( f.readline().strip() ) )
        f.close()
    print runtime
    _legend = [ job[len(_use_param.case):-1] for job in _use_param.jobs ]
    _colors = [ "orange", "purple", "green",  "red",  "blue", "black",  "yellow", "cyan", "magenta"]

    
    ld.myfigures ( _use_param.small_delta, processed_data , _use_param.case+"tmrca", _legend, _colors)
    ld.time_figure ( ld.calculate_acurrcy_array(processed_data, _use_param.small_delta) , runtime , _use_param.case+"tmrca", _legend, _colors)

    
    print "Done"
