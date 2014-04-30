#!/usr/bin/env python
import sys
if __name__ == "__main__":
    prefix = sys.argv[1]
    seqlen = float(sys.argv[2])
    
    freq_in_name  = prefix + "change"
    freq_out_name = prefix + "TreeFreq"
    change = [ float(x)*seqlen for x in open ( freq_in_name, "r" ) ]
    change.append(seqlen)
    #print change
    
    freq_out = open( freq_out_name, "w")
    freq = [] 
    for i in range(1, len(change)):
        freq.append( change[i] - change[i-1]) 
        freq_out.write(`int(change[i] - change[i-1])` + "\n")
    freq_out.close()
    #print numpy.sum(freq)
