#!/usr/bin/env python

import sys

def extract_tmrca ( file_name ):
    tmrca_file = open( file_name, "r")
    tmrca = [ float(x) for x in tmrca_file.read().split() ]
    tmrca_file.close()
    return tmrca


def compute_tmrca_moment( prefix , out ):
    tmrca = extract_tmrca ( prefix+"Tmrca" )
    duration = extract_tmrca ( prefix+"TreeFreq" )
    seqlen = sum(duration)
    
    moment = []
    for power in range (1, 5):
        ith_moment = 0
        for i, tmrca_i in enumerate( tmrca ):
            ith_moment += tmrca_i**power * duration[i] / seqlen
            
        moment.append ( ith_moment )
    outfile = open(out, "aw" )
    outfile.write ( `moment`.strip("[").strip("]") + '\n' )
    outfile.close()
    
    
if __name__ == "__main__":
    _in =  sys.argv[1] 
    _out = sys.argv[2] 
    compute_tmrca_moment( _in, _out)

