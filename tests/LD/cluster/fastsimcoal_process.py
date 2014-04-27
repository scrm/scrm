#!/usr/bin/env python
import sys
if __name__ == "__main__":
    prefix = sys.argv[1]
    tmrca_in_name  = prefix + "Tmrcaraw"
    tmrca_out_name = prefix + "Tmrca"
    _4Ne = float(40000.0)
    tmrca = [ float(x)/_4Ne for x in open( tmrca_in_name, "r" )]
    
    tmrca_out = open( tmrca_out_name, "w")
    for tmrca_i in tmrca:
        tmrca_out.write(`tmrca_i`+"\n")
    tmrca_out.close()
    
    freq_in_name  = prefix + "change"
    freq_out_name = prefix + "TreeFreq"
    change = [ int(x) for x in open ( freq_in_name, "r" ) ]
    #print change
    
    freq_out = open( freq_out_name, "w")
    freq = [] 
    for i in range(1, len(change)):
        freq.append( change[i] - change[i-1]) 
        freq_out.write(`change[i] - change[i-1]` + "\n")
    freq_out.close()
    #print freq
